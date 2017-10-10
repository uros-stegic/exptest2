#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double par_sum_we(const std::vector<double>& x, double mean, int start, int end)
{
	double sum{0};
	for(int i = start; i < end; i++) {
		sum += std::pow((x[i] - mean), 2);
	}
	return sum;
}

// [[Rcpp::export]]
Rcpp::NumericVector we_test(Rcpp::NumericVector x) {
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

		auto p_sum = par_sum_we(__x, mean, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = (1.0 / (length - 1)) * sum_vec(y_tmp);

#else
	auto y = (1.0 / (length - 1)) * par_sum_we(__x, mean, 0, length);
#endif

	auto r = (length - 1) * y / (std::pow(length, 2)*std::pow(mean, 2));
	result[0] = r;
	return result;
}
