#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> jacobiMatrixBasedOldMatrix(std::vector<std::vector<double>> matrixA, int i, int j)
{
  std::vector<std::vector<double>> matrixJij = linalg::identityMatrix(matrixA.size());
  double theta, error = 0.000001;

  if (abs(matrixA[i][j]) <= error)
  {
    return matrixJij;
  }
  if (abs(matrixA[i][i] - matrixA[j][j]) <= error)
  {
    // PI = atan(1.0) * 4
    theta = atan(1.0);
  }
  else
  {
    theta = atan((-2 * matrixA[i][j]) / (matrixA[i][i] - matrixA[j][j])) / 2;
  }
  matrixJij[i][i] = cos(theta);
  matrixJij[j][j] = cos(theta);
  matrixJij[i][j] = sin(theta);
  matrixJij[j][i] = -sin(theta);

  return matrixJij;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jacobiScan(std::vector<std::vector<double>> matrixA)
{
  std::vector<std::vector<double>> matrixJ = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixJij;

  for (size_t j = 0; j < matrixA[0].size() - 1; j++)
  {
    for (size_t i = j + 1; i < matrixA.size(); i++)
    {
      matrixJij = jacobiMatrixBasedOldMatrix(matrixAOld, i, j);
      matrixANew = linalg::matrixMultiplication(
          linalg::matrixMultiplication(
              linalg::transpose(matrixJij), matrixAOld),
          matrixJij);
      matrixAOld = matrixANew;
      matrixJ = linalg::matrixMultiplication(matrixJ, matrixJij);
    }
  }

  return std::make_tuple(matrixANew, matrixJ);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::jacobiMethod(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixJ;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;
  do
  {
    std::tie(matrixANew, matrixJ) = jacobiScan(matrixAOld);
    matrixAOld = matrixANew;
    matrixP = linalg::matrixMultiplication(matrixP, matrixJ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
  } while (val > toleranceError);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::jacobiHouseholder(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld, matrixH;
  std::tie(matrixAOld, matrixH) = matrixEigenvalue::householder(matrixA);
  std::vector<std::vector<double>> matrixANew, matrixJ;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;

  do
  {
    std::tie(matrixANew, matrixJ) = jacobiScan(matrixAOld);
    matrixAOld = matrixANew;
    matrixP = linalg::matrixMultiplication(matrixP, matrixJ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
  } while (val > toleranceError);

  matrixP = linalg::matrixMultiplication(matrixH, matrixP);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}