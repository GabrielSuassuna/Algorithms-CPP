#include "../linearAlgebra.h"

std::vector<std::vector<double>> linalg::identityMatrix(int n)
{
  std::vector<std::vector<double>> result;

  for (int i = 0; i < n; i++)
  {
    result.push_back({});
    for (int j = 0; j < n; j++)
    {
      if (i == j)
      {
        result[i].push_back(1);
      }
      else
      {
        result[i].push_back(0);
      }
    }
  }

  return result;
}