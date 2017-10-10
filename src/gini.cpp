#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

std::pair<double, double> par_sum_gini(const std::vector<double>& x, int start, int end)
{
	double sum{0}, sumx{0};
	auto n = x.size();
	for(int i = start; i < end; i++) {
		if(i != 0){
			sum += (i * (n - i)) * (x[i] - x[i-1]);
		}
		sumx += x[i];
	}
	return std::make_pair(sum, sumx);
}

// [[Rcpp::export]]
Rcpp::NumericVector gini_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	Rcpp::NumericVector result{1};
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

		auto p_sum = par_sum_gini(__x, start, end);

		y_tmp.push_back(p_sum.first);
		z_tmp.push_back(p_sum.second);
	}

	auto y = sum_vec(y_tmp);
	auto z = sum_vec(z_tmp);

#else
	auto p_sum = par_sum_gini(__x, 0, length);
	auto y = p_sum.first;
	auto z = p_sum.second;
#endif

	auto r = y / ((length - 1) * z);
	result[0] = r;
	return result;
}
