#include <Rcpp.h>
#include <omp.h>
#include <boost/math/special_functions/factorials.hpp>
#include "utils.hpp"

std::pair<double, double> par_sum_rossberg(const std::vector<double>& x, int start, int end)
{
	double sumh{0}, sumg{0};
	int h, g, n = x.size();
	double fac, fac2, fac3, min, max;
	fac = boost::math::factorial<double>(n);
	fac2 = boost::math::factorial<double>(n - 2);
	fac3 = boost::math::factorial<double>(n - 3);

	for(int m = start; m < end; m++) {
		h = 0;
		g = 0;
		for(int i = 0; i < n; i++) {
			for(int j = i+1; j < n; j++) {
				if(std::min(x[i], x[j]) < x[m]) {
					g++;
				}

				for(int k = j+1; k < n; k++) {
					min = std::min({x[i], x[j], x[k]});
					max = std::max({x[i], x[j], x[k]});

					if((x[i] + x[j] + x[k] - 2 * min - max ) < x[m]) {
						h++;
					}
				}
			}
		}
		sumh += (6 * fac3 / fac) * h;
		sumg += (2 * fac2 / fac) * g;
	}
	return std::make_pair(sumh, sumg);
}

// [[Rcpp::export]]
Rcpp::NumericVector rossberg_test(Rcpp::NumericVector x) {
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

		auto p_sum = par_sum_rossberg(__x, start, end);

		y_tmp.push_back(p_sum.first);
		m_tmp.push_back(p_sum.second);
	}

	auto y = sum_vec(y_tmp);
	auto m = sum_vec(m_tmp);

#else
	auto p_sum = par_sum_rossberg(__x, 0, length);
	auto y = p_sum.first;
	auto m = p_sum.second;
#endif

	Rcpp::NumericVector result{1};
	auto r = (y - m) / length;
 	result[0] = r;
	return result;
}
