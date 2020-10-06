#include "../linearAlgebra.h"

std::vector<double> linalg::normalize(std::vector<double> vector)
{
  std::vector<double> resultVector = vector;
  double norm = std::sqrt(linalg::dotProduct(vector, vector));

  for (size_t i = 0; i < resultVector.size(); i++)
  {
    resultVector[i] = resultVector[i] / norm;
  }

  return resultVector;
}