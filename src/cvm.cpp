#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include "utils.hpp"

/**
 * Eps <= 10^-4
 */
double par_sum_cvm(const std::vector<double>& x, int start, int end)
{
	double sum{0};
	auto n = x.size();
	for(int i = start; i < end; i++) {
		sum += std::pow(x[i] - (2.0 * (i + 1) - 1) / (2 * n), 2);
	}
	return sum;
}

void make_sorted_vec(std::vector<double>& x)
{
	auto mean = mean_vec(x);
	auto n = x.size();
	for(int i = 0; i < n; i++)
		x[i] = 1 - std::exp(-1 * x[i] / mean);

	std::sort(x.begin(), x.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector cvm_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	Rcpp::NumericVector result{1};

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;

	std::vector<double> y_tmp;
	omp_set_num_threads(4);
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	make_sorted_vec(__x);

	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_cvm(__x, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = sum_vec(y_tmp);
#else
	auto y = par_sum_cvm(__x, 0, length);
#endif

	auto r = 1.0 / (12 * length) + y;
	result[0] = r;
	return result;
}
