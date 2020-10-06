#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::matrixSubtraction(std::vector<std::vector<double>> matrixA, std::vector<std::vector<double>> matrixB)
{
  std::vector<std::vector<double>> result = matrixA;

  for (size_t i = 0; i < matrixA.size(); i++)
  {
    for (size_t j = 0; j < matrixA[i].size(); j++)
    {
      result[i][j] = matrixA[i][j] - matrixB[i][j];
    }
  }

  return result;
}