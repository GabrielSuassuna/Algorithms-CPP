#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::transpose(std::vector<std::vector<double>> matrix)
{
  std::vector<std::vector<double>> result = matrix;

  for (size_t i = 0; i < matrix.size(); i++)
  {
    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      result[i][j] = matrix[j][i];
    }
  }

  return result;
}