#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

std::pair<double, double> par_sum_gnedenko(const std::vector<double>& x, int r, int start, int end)
{
	double sum1{0}, sum2{0};
	auto n = x.size();
	for(int i = start; i < end; i++) {
		if(i == 0){
			sum1 += x[i] * n;
		}
		else if(i < r) {
			sum1 += (x[i] - x[i - 1]) * (n - i);
		}
		else {
			sum2 += (x[i] - x[i - 1]) * (n - i);
		}
	}
	return std::make_pair(sum1, sum2);
}

// [[Rcpp::export]]
Rcpp::NumericVector gnedenko_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	Rcpp::NumericVector result{1};
	auto length = __x.size();
	auto r = length / 2;
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	std::sort(__x.begin(), __x.end());

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;

	std::vector<double> y_tmp;
	std::vector<double> z_tmp;
	omp_set_num_threads(4);

	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_gnedenko(__x, r, start, end);

		y_tmp.push_back(p_sum.first);
		z_tmp.push_back(p_sum.second);
	}

	auto y = sum_vec(y_tmp) / r;
	auto z = sum_vec(z_tmp) / (length - r);
#else
	auto p_sum = par_sum_gnedenko(__x, r, 0, length);
	auto y = p_sum.first/r;
	auto z = p_sum.second/(length-r);
#endif

	result[0] = y / z;
	return result;
}
