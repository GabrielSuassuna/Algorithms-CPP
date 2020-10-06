#include "../linearAlgebra.h"

std::vector<double> linalg::vectorSubtraction(std::vector<double> x, std::vector<double> y)
{
  std::vector<double> result;

  for (size_t i = 0; i < x.size(); i++)
  {
    result.push_back(x[i] - y[i]);
  }

  return result;
}