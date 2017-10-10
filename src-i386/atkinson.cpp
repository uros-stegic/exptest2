#include <Rcpp.h>
#include <omp.h>
#include <boost/math/special_functions/gamma.hpp>
#include "utils.hpp"

std::pair<double, double> par_sum_atkinson(const std::vector<double>& x, double p, int start, int end)
{
	double real_mean{0};
	double power_mean{0};
	for(int i = start; i < end; i++) {
		real_mean += x[i];
		power_mean += std::pow(x[i], p);
	}
	return std::make_pair(real_mean, power_mean);
}

// [[Rcpp::export]]
Rcpp::NumericVector atkinson_test(Rcpp::NumericVector x, double p = 0.99) {
	auto __x = Rcpp::as<std::vector<double>>(x);
	auto length = __x.size();

#ifdef USE_OPENMP
	int thread_num, num_threads, start, end;

	omp_set_num_threads(4);
	std::vector<double> y_tmp;
	std::vector<double> m_tmp;
	num_threads = omp_get_num_threads();
	#pragma omp parallel private(thread_num,start,end)
	{
		thread_num = omp_get_thread_num();
		start = thread_num * length / num_threads;
		end = (thread_num + 1) * length / num_threads;

		auto p_sum = par_sum_atkinson(__x, p, start, end);

		y_tmp.push_back(p_sum.first);
		m_tmp.push_back(p_sum.second);
	}

	auto y = sum_vec(y_tmp)/length;
	auto m = sum_vec(m_tmp)/length;
#else
	auto p_sum = par_sum_atkinson(__x, p, 0, length);
	auto y = p_sum.first/length;
	auto m = p_sum.second/length;
#endif

	Rcpp::NumericVector result{1};
	auto r = std::sqrt(length)*std::abs(std::pow(m, 1/p)/y - std::pow(boost::math::tgamma(1+p), 1/p));
	result[0] = r;
	return result;
}
