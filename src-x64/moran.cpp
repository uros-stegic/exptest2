#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include "utils.hpp"

void make_log_vec(std::vector<double>& x, double mean, int start, int end)
{
	for(int i = start; i < end; i++) {
		x[i] = std::log(x[i] / mean);
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector moran_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	Rcpp::NumericVector result{1};
	auto length = __x.size();
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	auto mean = mean_vec(__x);

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

		make_log_vec(__x, mean, start, end);
	}

#else
	make_log_vec(__x, mean, 0, length);
#endif

	auto r = - boost::math::digamma(1) + mean_vec(__x);
	result[0] = r;
	return result;
}
