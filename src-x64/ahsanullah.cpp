#include <Rcpp.h>
#include <omp.h>
#include "utils.hpp"

double par_sum_ahsanullah(const std::vector<double>& x, int start, int end)
{
	double result{0};
	auto n = x.size();
	for(int t = start; t < end; t++) {
		for(unsigned i = 0; i < n; i++) {
			for(unsigned j = 0; j < n; j++) {
				result += (std::fabs(x[i] - x[j]) < x[t]) - (2 * std::min(x[i], x[j]) < x[t]);
			}
		}
	}
	return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector ahsanullah_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();
	double cube = length * length * length;

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;
	omp_set_num_threads(4);
	std::vector<double> res;
	num_threads = omp_get_num_threads();

	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		res.push_back(par_sum_ahsanullah(__x, start, end));
	}

	auto fin = sum_vec(res);
#else
	auto fin = par_sum_ahsanullah(__x, 0, length);
#endif

	Rcpp::NumericVector result{1};
	result[0] = fin / cube;
 	return result;
}
