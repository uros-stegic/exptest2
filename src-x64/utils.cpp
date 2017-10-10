#include "utils.hpp"

double sum_vec(const std::vector<double>& vec)
{
	double res{0};
	for(auto i = 0u; i < vec.size(); i++)
	{
		res += vec[i];
	}
	return res;
}


double mean_vec(const std::vector<double>& vec)
{
	return sum_vec(vec) / vec.size();
}

