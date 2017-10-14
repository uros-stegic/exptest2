#include <Rcpp.h>
#include "utils.hpp"

/**
 * Integral type goodness-of-fit exponentiality test from the paper.
 */
// [[Rcpp::export]]
Rcpp::NumericVector integral_test(Rcpp::NumericVector x)
{
	double tmp{0.0};
	auto _x = Rcpp::as<std::vector<double>>(x);
	for(unsigned int i = 1; i <= x.size(); i++) {
		tmp += (max_sum(_x, i) - median_sum(_x, i));
	}
	Rcpp::NumericVector result{1};
	result[0] = tmp / x.size();
	return result;
}
