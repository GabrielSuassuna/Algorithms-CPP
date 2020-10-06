#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

namespace linalg
{
  void printVector(std::vector<double> vector);
  double dotProduct(std::vector<double> x, std::vector<double> y);
  std::vector<double> normalize(std::vector<double> vector);
  std::vector<double> saxpy(std::vector<double> vector, double alpha);
  std::vector<double> gaxpy(std::vector<std::vector<double>> matrix, std::vector<double> vector);
  std::vector<std::vector<double>> matrixSubtraction(std::vector<std::vector<double>> matrixA, std::vector<std::vector<double>> matrixB);
  std::vector<std::vector<double>> scalarMatrixMultiplication(std::vector<std::vector<double>> matrix, double alpha);
  std::vector<std::vector<double>> identityMatrix(int n);
} // namespace linalg