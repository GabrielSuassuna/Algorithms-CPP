#include "../linearAlgebra.h"

void linalg::printVector(std::vector<double> vector)
{
  std::cout << "[ ";
  for (size_t i = 0; i < vector.size(); i++)
  {
    if (i != vector.size() - 1)
    {
      std::cout << vector[i] << ", ";
    }
    else
    {
      std::cout << vector[i];
    }
  }
  std::cout << " ]" << std::endl;
}
