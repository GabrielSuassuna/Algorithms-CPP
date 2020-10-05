#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

namespace linalg
{
  void printVector(std::vector<double> vector);
  double dotProduct(std::vector<double> x, std::vector<double> y);
  std::vector<double> normalize(std::vector<double> vector);
  std::vector<double> gaxpy(std::vector<std::vector<double>> matrix, std::vector<double> vector);
}