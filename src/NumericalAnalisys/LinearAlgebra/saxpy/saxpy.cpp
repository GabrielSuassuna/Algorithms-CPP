#include "../linearAlgebra.h"

std::vector<double> linalg::saxpy(std::vector<double> vector, double alpha)
{
  std::vector<double> result;

  for (size_t i = 0; i < vector.size(); i++)
  {
    result.push_back(alpha * vector[i]);
  }

  return result;
}