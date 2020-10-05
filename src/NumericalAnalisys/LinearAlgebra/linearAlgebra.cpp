#include "linearAlgebra.h"

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

double linalg::dotProduct(std::vector<double> x, std::vector<double> y)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error(std::string("Vectors with diferent sizes"));
  }

  double sum = 0;

  for (size_t i = 0; i < x.size(); i++)
  {
    sum = sum + x[i] * y[i];
  }

  return sum;
}

// void linalg::saxpy(double a, std::vector<double> x, std::vector<double> y)
// {
//   if (x.size() != y.size())
//   {
//     throw std::runtime_error(std::string("Vectors with diferent sizes"));
//   }

//   for (size_t i = 0; i < x.size(); i++)
//   {
//     y[i] = y[i] + a * x[i];
//   }
// }

std::vector<double> linalg::gaxpy(std::vector<std::vector<double>> matrix, std::vector<double> vector)
{
  std::vector<double> y;
  double sum;

  for (size_t i = 0; i < matrix.size(); i++)
  {
    sum = 0;

    for (size_t j = 0; j < matrix[i].size(); j++)
    {
      sum = sum + matrix[i][j] * vector[i];
    }

    y.push_back(sum);
  }

  return y;
}

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