#include "../linearAlgebra.h"

std::vector<double> linalg::gaxpy(std::vector<std::vector<double>> matrix, std::vector<double> vector)
{
  std::vector<double> result(matrix.size(), 0.0);

  for (size_t i = 0; i < matrix.size(); i++)
  {
    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      result[i] += matrix[i][j] * vector[j];
    }
  }

  return result;
}