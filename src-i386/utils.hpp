#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <array>
#include <vector>

/* Creates alias for vector of three elements of type */
using tripplet = std::array<double, 3>;

/**
 * Helper function that sums a vector of type T.
 */
double sum_vec(const std::vector<double>&);

/**
 * Helper function that calculates vectors mean.
 */
double mean_vec(const std::vector<double>&);

#endif // __UTILS_HPP__

