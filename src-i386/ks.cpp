#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

double max_ks(const std::vector<double>& x,  int start, int end)
{
	double max_el, a, b;
	auto n = x.size();
	a = (start + 1.0) / n - x[start];
	b = x[start] - (start / (n * 1.0));
	max_el = (a > b) ? a : b;
	for(int i = start + 1; i < end; i++) {
		a = (i + 1.0) / n - x[i];
		b = x[i] - (i / (n * 1.0));
        if(a > b && a > max_el) {
        	max_el = a;
       	}
      	else if(b > a && b > max_el) {
      		max_el = b;
      	}
	}
	return max_el;
}

void make_sorted_ks(std::vector<double>& x, double mean) {
	auto n = x.size();
	for(unsigned i = 0; i < n; i++) {
		x[i] = 1 - std::exp(-(x[i] / mean));
	}
	std::sort(x.begin(), x.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector ks_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);

	Rcpp::NumericVector result{1};
	std::vector<double> y_tmp;
	auto length = __x.size();
	if(length == 0)
	{
		result[0] = 0;
		return result;
	}

	auto mean = mean_vec(__x);
	make_sorted_ks(__x, mean);

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;
	omp_set_num_threads(4);
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = max_ks(__x, start, end);

		y_tmp.push_back(p_sum);
	}

	auto r = *(std::max_element(y_tmp.begin(), y_tmp.end()));
#else
	auto r = max_ks(__x, 0, length);
#endif

	result[0] = r;
	return result;
}
