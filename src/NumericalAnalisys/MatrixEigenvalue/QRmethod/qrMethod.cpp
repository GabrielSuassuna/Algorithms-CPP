#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> jacobiMatrixBasedOldRMatrix(std::vector<std::vector<double>> matrixA, int i, int j)
{
  std::vector<std::vector<double>> matrixJij = linalg::identityMatrix(matrixA.size());
  double theta, error = 0.000001;

  if (abs(matrixA[i][j]) <= error)
  {
    return matrixJij;
  }
  if (abs(matrixA[j][j]) <= error)
  {
    if (matrixA[i][j] < 0)
    {
      // PI/2 = atan(1.0) * 2
      theta = atan(1.0) * 2;
    }
    else
    {

      // PI/2 = atan(1.0) * 2
      theta = -atan(1.0) * 2;
    }
  }
  else
  {
    theta = atan(-matrixA[i][j] / matrixA[j][j]);
  }
  matrixJij[i][i] = cos(theta);
  matrixJij[j][j] = cos(theta);
  matrixJij[i][j] = sin(theta);
  matrixJij[j][i] = -sin(theta);

  return matrixJij;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> decompositionQR(std::vector<std::vector<double>> matrixA)
{
  std::vector<std::vector<double>> matrixQT = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixROld = matrixA;
  std::vector<std::vector<double>> matrixRNew, matrixJij, matrixQ, matrixR;

  for (size_t j = 0; j < matrixA[0].size() - 1; j++)
  {
    for (size_t i = j + 1; i < matrixA.size(); i++)
    {
      matrixJij = jacobiMatrixBasedOldRMatrix(matrixROld, i, j);
      matrixRNew = linalg::matrixMultiplication(matrixJij, matrixROld);
      matrixROld = matrixRNew;
      matrixQT = linalg::matrixMultiplication(matrixJij, matrixQT);
    }
  }

  matrixQ = linalg::transpose(matrixQT);
  matrixR = matrixRNew;

  return std::make_tuple(matrixQ, matrixR);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::qrMethod(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixQ, matrixR;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;
  do
  {
    std::tie(matrixQ, matrixR) = decompositionQR(matrixAOld);
    matrixANew = linalg::matrixMultiplication(matrixR, matrixQ);
    matrixAOld = matrixANew;
    matrixP = linalg::matrixMultiplication(matrixP, matrixQ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
  } while (val > toleranceError);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::qrHouseholder(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixANew, matrixQ, matrixR, matrixAOld, matrixH;
  std::tie(matrixAOld, matrixH) = matrixEigenvalue::householder(matrixA);
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;
  do
  {
    std::tie(matrixQ, matrixR) = decompositionQR(matrixAOld);
    matrixANew = linalg::matrixMultiplication(matrixR, matrixQ);
    matrixAOld = matrixANew;
    matrixP = linalg::matrixMultiplication(matrixP, matrixQ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
  } while (val > toleranceError);

  matrixP = linalg::matrixMultiplication(matrixH, matrixP);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}