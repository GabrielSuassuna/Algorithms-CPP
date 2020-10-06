#include "../linearAlgebra.h"

std::vector<double> linalg::gaxpy(std::vector<std::vector<double>> matrix, std::vector<double> vector)
{
  std::vector<double> result;
  double sum;

  for (size_t i = 0; i < matrix.size(); i++)
  {
    sum = 0;

    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      sum = sum + matrix[i][j] * vector[i];
    }

    result.push_back(sum);
  }

  return result;
}