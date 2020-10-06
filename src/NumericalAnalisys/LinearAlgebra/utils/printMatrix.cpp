#include "../linearAlgebra.h"

void linalg::printMatrix(std::vector<std::vector<double>> matrix)
{
  for (size_t i = 0; i < matrix.size(); i++)
  {
    std::cout << "| ";
    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      if (j != matrix[i].size() - 1)
      {
        std::cout << std::setprecision(6) << matrix[i][j] << ", ";
      }
      else
      {
        std::cout << std::setprecision(6) << matrix[i][j];
      }
    }
    std::cout << " |" << std::endl;
  }
}
