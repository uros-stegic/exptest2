#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

void par_sum_kimber(std::vector<double>& x, double mean, int start, int end)
{
	double pi = 2 / 3.1415;
	auto n = x.size();
	for(int i = start; i < end; i++) {
		x[i] = std::abs(pi * std::asin(std::sqrt((i + 0.5) / n))
                      - pi * std::asin(std::sqrt(1 - (std::exp(-(x[i] / mean))))));
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector kimber_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);


	Rcpp::NumericVector result{1};
	auto length = __x.size();
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	auto mean = mean_vec(__x);
	std::sort(__x.begin(), __x.end());

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;
	omp_set_num_threads(4);
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		par_sum_kimber(__x, mean, start, end);
	}

#else
	par_sum_kimber(__x, mean, 0, length);
#endif

	auto r = std::max_element(__x.begin(), __x.end());
	result[0] = *r;
	return result;
}
