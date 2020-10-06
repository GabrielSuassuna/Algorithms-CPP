#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> holseholderMethod(std::vector<std::vector<double>> matrix, int i)
{
  std::vector<double> w(matrix.size(), 0.0);
  std::vector<double> wL(matrix.size(), 0.0);

  for (size_t position = i + 1; position < matrix.size(); position++)
  {
    w[position] = matrix[position][i];
  }

  double wSize = std::sqrt(linalg::dotProduct(w, w));

  wL[i + 1] = wSize;

  std::vector<double> n = linalg::normalize(linalg::vectorSubtraction(w, wL));

  return linalg::matrixSubtraction(linalg::identityMatrix(matrix.size()), linalg::scalarMatrixMultiplication(linalg::vectorOuterProduct(n, n), 2));
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> matrixEigenvalue::householder(std::vector<std::vector<double>> matrix)
{
  std::vector<std::vector<double>> matrixHAux, matrixABar;
  std::vector<std::vector<double>> matrixH = linalg::identityMatrix(matrix.size());
  std::vector<std::vector<double>> matrixA = matrix;

  for (size_t i = 0; i < matrix.size() - 2; i++)
  {
    matrixHAux = holseholderMethod(matrixA, i);
    matrixABar = linalg::matrixMultiplication(linalg::matrixMultiplication(linalg::transpose(matrixHAux), matrixA), matrixHAux);
    matrixA = matrixABar;
    matrixH = linalg::matrixMultiplication(matrixH, matrixHAux);
  }

  return std::make_tuple(matrixABar, matrixH);
}