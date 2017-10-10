#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double par_sum_deshpande(const std::vector<double>& x, double b, int start, int end)
{
	int sum{0}, n = x.size();
	for(int i = start; i < end; i++) {
		for(int j = i + 1; j < n; j++) {
			if(x[i] > b * x[j])
				sum++;
			if(x[j] > b * x[i])
				sum++;
		}
	}
	return sum;
}

// [[Rcpp::export]]
Rcpp::NumericVector deshpande_test(Rcpp::NumericVector x, double b = 0.44) {
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

	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_deshpande(__x, b, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = sum_vec(y_tmp);

#else
	auto y = par_sum_deshpande(__x, b, 0, length);
#endif

	auto r = 1.0 / (length * (length - 1)) * y;
	result[0] = r;
	return result;
}
