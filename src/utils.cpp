#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <omp.h>

double sum_vec(const std::vector<double>& vec)
{
	double res{0};
	for(auto i = 0u; i < vec.size(); i++)
	{
		res += vec[i];
	}
	return res;
}


double mean_vec(const std::vector<double>& vec)
{
	return sum_vec(vec) / vec.size();
}

double median(const tripplet& nums)
{
	if( nums[0] > nums[1] ) {
		if( nums[0] < nums[2] ) {
			return nums[0];
		}
		else if( nums[2] > nums[1] ) {
			return nums[2];
		}
		else {
			return nums[1];
		}
	}
	else {
		if( nums[1] < nums[2] ) {
			return nums[1];
		}
		else if( nums[2] > nums[0] ) {
			return nums[2];
		}
		else {
			return nums[0];
		}
	}
}

/**
 * Helper function that does calculates partial sum for the function G_n(t) (from the paper)
 * with the bounds [start, end). It's intent is to be used to calculate nested sums from
 * function G_n(t) in parallel and thus gaining significant performance boost.
 */
double median_sumation(const std::vector<double>& x, int n, int start, int end, double t)
{
	double partial_sum = 0;
	for(auto i = start; i < end; i++) {
		for(auto j = 0; j < n; j++) {
			for(auto k = j; k < n; k++) {
				for(auto l = k; l < n; l++) {
					double med = median(tripplet{x[j], x[k], x[l]});
					auto tmp = x[i] + med < t;
					if( j == k && k == l) {
						partial_sum += tmp;
					}
					else if( j != k && j != l && k != l ) {
						partial_sum += tmp*6;
					}
					else {
						partial_sum += tmp*3;
					}
				}
			}
		}
	}
	return partial_sum;
}

/**
 * Function G_n(t) from the paper. Uses OpenMP multithreading library. Currently,
 * it uses predefined number of threads (4) because this is still testing environment
 * and that is the number of threads that are available on the development machine.
 */
double median_sum(const std::vector<double>& x, double t)
{
	auto length = x.size();
	auto squared = length*length;
	double quad = squared * squared;

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;
	omp_set_num_threads(4);
	num_threads = omp_get_num_threads();
	std::vector<double> res;
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		res.push_back(median_sumation(x, length, start, end, t));
	}
	double sum = sum_vec(res) / quad;
#else
	auto res = median_sumation(x, length, 0, length-1, t);
	double sum = res / quad;
#endif

	return sum;
}

/**
 * Helper function that does calculates partial sum for the function K_n(t) (from the paper).
 * with the bounds [start, end). It's intent is to be used to calculate nested sums from
 * function K_n(t) in parallel and thus gaining significant performance boost.
 */
double max_sumation(const std::vector<double>& x, int n, int start, int end, double t)
{
	double partial_sum = 0;
	for(auto i = start; i < end; i++) {
		for(auto j = i; j < n; j++) {
			for(auto k = j; k < n; k++) {
				double max = std::max({x[j], x[k], x[i]});
				auto tmp = max < t;
				if( i == j && j == k ) {
					partial_sum += tmp;
				}
				else if( i != j && i != k && j != k ) {
					partial_sum += tmp*6;
				}
				else {
					partial_sum += tmp*3;
				}
			}
		}
	}
	return partial_sum;
}

/**
 * Function K_n(t) from the paper. Uses OpenMP multithreading library. Currently,
 * it uses predefined number of threads (4) because this is still testing environment
 * and that is the number of threads that are available on the development machine.
 */
double max_sum(const std::vector<double>& x, double t)
{
	auto length = x.size();
	double cube = length * length * length;

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;
	omp_set_num_threads(4);
	std::vector<double> res;
	num_threads = omp_get_num_threads();

	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		res.push_back(max_sumation(x, length, start, end, t));
	}

	double sum = sum_vec(res) / cube;
#else
	double sum = max_sumation(x, length, 0, length, t) / cube;
#endif

	return sum;
}
