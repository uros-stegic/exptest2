#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double par_sum_frozini(const std::vector<double>& x, double mean, int start, int end)
{
	double sum{0};
	auto n = x.size();
	for(int i = start; i < end; i++) {
		sum += std::abs(1 - std::exp(-(x[i] / mean)) - (i + 1 - 0.5) / n);
	}
	return sum;
}

// [[Rcpp::export]]
Rcpp::NumericVector frozini_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	Rcpp::NumericVector result{1};
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	auto mean = mean_vec(__x);
	std::sort(__x.begin(), __x.end());

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;

	std::vector<double> y_tmp;
	omp_set_num_threads(4);
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_frozini(__x, mean, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = sum_vec(y_tmp);

#else
	auto y = par_sum_frozini(__x, mean, 0, length);
#endif

	auto r = 1.0 / std::sqrt(length) * y;
	result[0] = r;
	return result;
}
