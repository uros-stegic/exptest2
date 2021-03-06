#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double par_sum_hollander(const std::vector<double>& x, int start, int end)
{
	int sum = 0;
	for(int i = start; i < end; i++) {
		for(int j = i + 1; j < end; j++) {
			for(int k = j + 1; k < end; k++) {
				if(i != j && i != k && j < k && (x[i] > (x[j] + x[k]))) {
					sum++;
				}
				else if(j != i && j != k && i < k && (x[j] > (x[i] + x[k]))) {
					sum++;
				}
				else if(k != i && k != j && i < j && (x[k] > (x[i] + x[j]))) {
					sum++;
				}
			}
		}
	}
	return sum;
}

// [[Rcpp::export]]
Rcpp::NumericVector hollander_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	Rcpp::NumericVector result{1};
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

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

		auto p_sum = par_sum_hollander(__x, start, end);

		y_tmp.push_back(p_sum);
	}

	auto y = sum_vec(y_tmp);
#else
	auto y = par_sum_hollander(__x, 0, length);
#endif

	auto r = (2.0 / (length * (length - 1) * (length - 2))) * y;
	result[0] = r;
	return result;
}
