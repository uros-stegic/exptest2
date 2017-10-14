#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <array>
#include <vector>

/* Creates alias for vector of three elements of type */
using tripplet = std::array<double, 3>;

/*template<>
double std::max<tripplet>(const tripplet& t) {
    return std::max(std::max(t[0], t[1]), t[2]);
}*/

/**
 * Helper function that sums a vector of type T.
 */
double sum_vec(const std::vector<double>&);

/**
 * Helper function that calculates vectors mean.
 */
double mean_vec(const std::vector<double>&);

/**
 * Calculates median of three elements.
 */
double median(const tripplet& nums);

/**
 * Helper function that does calculates partial sum for the function G_n(t) (from the paper)
 * with the bounds [start, end). It's intent is to be used to calculate nested sums from
 * function G_n(t) in parallel and thus gaining significant performance boost.
 */
double median_sumation(const std::vector<double>& x, int n, int start, int end, double t);

/**
 * Function G_n(t) from the paper. Uses OpenMP multithreading library. Currently,
 * it uses predefined number of threads (4) because this is still testing environment
 * and that is the number of threads that are available on the development machine.
 */
double median_sum(const std::vector<double>& x, double t);

/**
 * Helper function that does calculates partial sum for the function K_n(t) (from the paper).
 * with the bounds [start, end). It's intent is to be used to calculate nested sums from
 * function K_n(t) in parallel and thus gaining significant performance boost.
 */
double max_sumation(const std::vector<double>& x, int n, int start, int end, double t);

/**
 * Function K_n(t) from the paper. Uses OpenMP multithreading library. Currently,
 * it uses predefined number of threads (4) because this is still testing environment
 * and that is the number of threads that are available on the development machine.
 */
double max_sum(const std::vector<double>& x, double t);

#endif // __UTILS_HPP__

