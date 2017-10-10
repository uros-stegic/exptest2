#include <Rcpp.h>
#include <omp.h>
#include <cmath>
#include "utils.hpp"

std::pair<double, double> max_min(const std::vector<double>& x) {
	auto max = x[0], min = x[0];
	auto n = x.size();
	for(unsigned i = 1; i < n; i++) {
		if(x[i] > max) {
			max = x[i];
		}
		if(x[i] < min) {
			min = x[i];
		}
	}
	return std::make_pair(max, min);
}
// [[Rcpp::export]]
Rcpp::NumericVector ww_test(Rcpp::NumericVector x) {
	auto __x = Rcpp::as<std::vector<double>>(x);

	Rcpp::NumericVector result{1};
	if(__x.size() == 0) {
		result[0] = 0;
	}
	else {
		auto res = max_min(__x);
		if(res.second != 0) {
			result[0] = res.first / res.second;
		}
		else {
			result[0] = 0;
		}
	}
	return result;
}
