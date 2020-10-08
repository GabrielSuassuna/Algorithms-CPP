#include "../linearAlgebra.h"

double linalg::sumSquaresTermsBelowDiagonal(std::vector<std::vector<double>> matrixA)
{
  double sum = 0;
  for (size_t i = 1; i < matrixA.size(); i++)
  {
    for (size_t j = 0; j < i; j++)
    {
      sum += pow(matrixA[i][j], 2);
    }
  }

  return sum;
}