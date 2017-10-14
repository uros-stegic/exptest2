#include <Rcpp.h>
#include <cmath>
#include "utils.hpp"

/**
 * Kolmogorov type goodness-of-fit exponentiality test from the paper.
 */
// [[Rcpp::export]]
Rcpp::NumericVector kolmogorov_test(Rcpp::NumericVector x)
{
	double tmp{-1};
	auto _x = Rcpp::as<std::vector<double>>(x);
	for(unsigned int i = 1; i <= x.size(); i++) {
		tmp = std::max(tmp, std::fabs(max_sum(_x, i) - median_sum(_x, i)));
	}
	Rcpp::NumericVector result{1};
	result[0] = tmp;
	return result;
}
