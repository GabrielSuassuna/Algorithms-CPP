#include "../linearAlgebra.h"

double linalg::dotProduct(std::vector<double> x, std::vector<double> y)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error(std::string("Vectors with diferent sizes"));
  }

  double sum = 0;

  for (size_t i = 0; i < x.size(); i++)
  {
    sum = sum + x[i] * y[i];
  }

  return sum;
}