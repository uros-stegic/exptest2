#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

std::pair<double, double>  par_sum_harris(const std::vector<double>& x, int r, int start, int end)
{
	double sum1{0}, sum2{0};
	auto n = x.size();
	int m = n - r;
	for(int i = start; i < end; i++) {
		if(i == 0) {
			sum1 += x[i] * n;
		}
		else if(i < r || i >= m) {
			sum1 += (x[i] - x[i - 1]) * (n - i);
		}
		else {
			sum2 += (x[i] - x[i - 1]) * (n - i);
		}
	}
	return std::make_pair(sum1, sum2);
}

// [[Rcpp::export]]
Rcpp::NumericVector harris_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	auto r = length / 4;
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
	std::vector<double> w_tmp;
	omp_set_num_threads(4);
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_harris(__x, r, start, end);

		y_tmp.push_back(p_sum.first);
		z_tmp.push_back(p_sum.second);
	}

	auto y = sum_vec(y_tmp);
	auto z = sum_vec(z_tmp);

#else
	auto p_sum = par_sum_harris(__x, r, 0, length);
	auto y = p_sum.first;
	auto z = p_sum.second;
#endif

	result[0] = (y / (2.0 * r)) / (z / (length - 2.0 * r));
	return result;
}
