#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::scalarMatrixMultiplication(std::vector<std::vector<double>> matrix, double alpha)
{
  std::vector<std::vector<double>> result = matrix;

  for (size_t i = 0; i < matrix.size(); i++)
  {
    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      result[i][j] = alpha * matrix[i][j];
    }
  }

  return result;
}