#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::matrixMultiplication(std::vector<std::vector<double>> matrixA, std::vector<std::vector<double>> matrixB)
{
  std::vector<std::vector<double>> result(matrixA.size(), std::vector<double>(matrixB[0].size(), 0.0));

  for (size_t i = 0; i < matrixA.size(); i++)
  {
    for (size_t j = 0; j < matrixB[0].size(); j++)
    {
      for (size_t k = 0; k < matrixA[i].size(); k++)
      {
        result[i][j] += matrixA[i][k] * matrixB[k][j];
      }
    }
  }

  return result;
}