#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double par_sum_hegazy2(const std::vector<double>& x, int start, int end)
{
	double sum{0};
	auto n = x.size();
	for(int i = start; i < end; i++) {
		sum += std::pow(x[i] + std::log(1 - ((i + 1) / (n + 1.0))), 2);
	}
	return sum;
}

// [[Rcpp::export]]
Rcpp::NumericVector hegazy2_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	Rcpp::NumericVector result{1};
	auto length = __x.size();
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;

	std::vector<double> y_tmp;
	omp_set_num_threads(4);

	std::sort(__x.begin(), __x.end());
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_hegazy2(__x, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = sum_vec(y_tmp);
#else
	auto y = par_sum_hegazy2(__x, 0, length);
#endif

	auto r = y / length;
	result[0] = r;
	return result;
}
