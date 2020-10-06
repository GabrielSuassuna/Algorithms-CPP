#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::vectorOuterProduct(std::vector<double> x, std::vector<double> y)
{
  std::vector<std::vector<double>> result(x.size(), std::vector<double>(y.size(), 0.0));

  for (size_t i = 0; i < x.size(); i++)
  {
    for (size_t j = 0; j < y.size(); j++)
    {
      result[i][j] = x[i] * y[j];
    }
  }

  return result;
}