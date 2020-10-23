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

  std::cout << "Matriz Q:" << std::endl;
  linalg::printMatrix(matrixQ);
  std::cout << "Matriz R:" << std::endl;
  linalg::printMatrix(matrixR);

  return std::make_tuple(matrixQ, matrixR);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::qrMethod(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixQ, matrixR;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;

  std::cout << "iii. Matriz que sai de cada varredura de Jacobi:" << std::endl;
  int count = 1;

  do
  {
    std::cout << std::endl
              << "Passo " << count << ":" << std::endl
              << std::endl;
    std::tie(matrixQ, matrixR) = decompositionQR(matrixAOld);
    matrixANew = linalg::matrixMultiplication(matrixR, matrixQ);
    matrixAOld = matrixANew;
    std::cout << "Matriz A:" << std::endl;
    linalg::printMatrix(matrixANew);
    matrixP = linalg::matrixMultiplication(matrixP, matrixQ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
    count++;
  } while (val > toleranceError);

  std::cout << std::endl
            << "i. Matriz diagonal:" << std::endl;
  linalg::printMatrix(matrixANew);

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

  std::cout << std::endl
            << "i. Imprima matriz Anova que sai de cada iteração QR" << std::endl;
  int count = 1;

  do
  {
    std::cout << std::endl
              << "Passo " << count << ":" << std::endl
              << std::endl;
    std::tie(matrixQ, matrixR) = decompositionQR(matrixAOld);
    matrixANew = linalg::matrixMultiplication(matrixR, matrixQ);
    matrixAOld = matrixANew;
    std::cout
        << "Matriz Anova:" << std::endl;
    linalg::printMatrix(matrixANew);
    matrixP = linalg::matrixMultiplication(matrixP, matrixQ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
    count++;
  } while (val > toleranceError);

  std::cout << std::endl
            << "ii. As colunas P não são autovetores de A:" << std::endl;
  linalg::printMatrix(matrixP);

  matrixP = linalg::matrixMultiplication(matrixH, matrixP);

  std::cout << std::endl
            << "iii. P = HP; As colunas P são autovetores de A:" << std::endl;
  linalg::printMatrix(matrixP);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}