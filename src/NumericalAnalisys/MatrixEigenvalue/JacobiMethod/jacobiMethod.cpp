#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> jacobiMatrixBasedOldMatrix(std::vector<std::vector<double>> matrixA, int i, int j)
{
  // Matriz identidade com n x n elementos
  std::vector<std::vector<double>> matrixJij = linalg::identityMatrix(matrixA.size());
  double theta, error = 0.000001;

  // Considerar Aij = 0, retornar matriz identidade
  if (abs(matrixA[i][j]) <= error)
  {
    return matrixJij;
  }

  // Considerar Aii = Ajj retornar pi/4;
  if (abs(matrixA[i][i] - matrixA[j][j]) <= error)
  {
    // PI = atan(1.0) * 4
    theta = atan(1.0);
  }

  // Esta função já retorna um ângulo +/-
  // no primeiro quadrante sentido anti-horário (+)
  // no primeiro quadrante sentido horário (-)
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
  // Matriz que contém os produtos das matrizes ortogonais Jij para recuperar os
  // autovetores da matriz original.
  std::vector<std::vector<double>> matrixJ = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixJij;

  // Para transformar a matriz original, precisamos zerar todos
  // os elementos abaixo do elemento da diagonal principal em cada coluna exceto na
  // última. Portanto, o loop j das colunas vai até a penúltima coluna. O loop
  // das linhas para uma dada coluna j, percorre todos os elementos abaixo da linha
  // da diagonal, ou seja, da linha (j+1) até a última linha (linha n). Como a matriz
  // é simétrica, o que fizermos para as colunas acontecerá igualmente nas linhas

  // Laço das colunas
  for (size_t j = 0; j < matrixA[0].size() - 1; j++)
  {
    // Laço das linhas
    for (size_t i = j + 1; i < matrixA.size(); i++)
    {
      // Construção da matriz de Jacobi Jij. A matriz de Jacobi é uma matriz de rotação.
      matrixJij = jacobiMatrixBasedOldMatrix(matrixAOld, i, j);

      // Transformação de similaridade do passo ij
      matrixANew = linalg::matrixMultiplication(
          linalg::matrixMultiplication(
              linalg::transpose(matrixJij), matrixAOld),
          matrixJij);

      // Salvar para o próximo passo.
      matrixAOld = matrixANew;

      // Acumular o produto das matrizes de Jacobi
      matrixJ = linalg::matrixMultiplication(matrixJ, matrixJij);
    }
  }

  // No final do loop externo, o formato da matriz A já está mais próximo
  // do formato de uma matriz diagonal.

  std::cout << "Matriz A:" << std::endl;
  linalg::printMatrix(matrixANew);
  std::cout << "Matriz J:" << std::endl;
  linalg::printMatrix(matrixJ);

  return std::make_tuple(matrixANew, matrixJ);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::jacobiMethod(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld = matrixA;
  std::vector<std::vector<double>> matrixANew, matrixJ;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;

  std::cout << "iii. Matriz que sai de cada varredura de Jacobi:" << std::endl;
  int count = 1;

  // Laço das varreduras de diagonalização
  do
  {
    std::cout << std::endl
              << "Passo " << count << ":" << std::endl
              << std::endl;
    // Varredura de Jacobi (devolve uma matriz que deve
    // se aproximar de uma matriz diagonal)
    std::tie(matrixANew, matrixJ) = jacobiScan(matrixAOld);

    // Salvar para o próximo passo.
    matrixAOld = matrixANew;

    // Acumular o produto das matrizes de Jacobi
    matrixP = linalg::matrixMultiplication(matrixP, matrixJ);

    // Verificar se a matriz nova já é diagonal
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
    count += 1;
  } while (val > toleranceError);

  // Ao sair do loop, o formato da matriz nova já está suficientemente próximo do formato de
  // uma matriz diagonal. Assim, os elementos da diagonal são os autovalores da matriz original
  // de entrada e as colunas da matriz P são os autovetores correspondente.
  std::cout << std::endl
            << "i. Matriz diagonal:" << std::endl;
  linalg::printMatrix(matrixANew);

  // Copia os elementos da diagonal da matriz no vetor Lamb
  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jacobiScanHouseholder(std::vector<std::vector<double>> matrixA)
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
      std::cout << std::endl
                << "Matriz A logo após a linha Anova <- JitT*Avelha*Jit:" << std::endl
                << std::endl;
      linalg::printMatrix(matrixANew);
      matrixAOld = matrixANew;
      matrixJ = linalg::matrixMultiplication(matrixJ, matrixJij);
    }
  }

  std::cout << std::endl
            << "Matriz A ao final da varredura:" << std::endl;
  linalg::printMatrix(matrixANew);
  std::cout << std::endl
            << "Matriz J ao final da varredura:" << std::endl;
  linalg::printMatrix(matrixJ);

  return std::make_tuple(matrixANew, matrixJ);
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> matrixEigenvalue::jacobiHouseholder(std::vector<std::vector<double>> matrixA, double toleranceError)
{
  std::vector<std::vector<double>> matrixP = linalg::identityMatrix(matrixA.size());
  std::vector<std::vector<double>> matrixAOld, matrixH;
  std::tie(matrixAOld, matrixH) = matrixEigenvalue::householder(matrixA);
  std::vector<std::vector<double>> matrixANew, matrixJ;
  std::vector<double> lamb(matrixA.size(), 0.0);
  double val;

  std::cout << std::endl
            << "ii. Imprima matriz Anova logo após a linha Anova <- JitT*Avelha*Jit:" << std::endl;

  do
  {
    std::tie(matrixANew, matrixJ) = jacobiScanHouseholder(matrixAOld);
    matrixAOld = matrixANew;
    matrixP = linalg::matrixMultiplication(matrixP, matrixJ);
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
  } while (val > toleranceError);

  std::cout << std::endl
            << "iii. As colunas P não são autovetores de A:" << std::endl;
  linalg::printMatrix(matrixP);

  matrixP = linalg::matrixMultiplication(matrixH, matrixP);

  std::cout << std::endl
            << "iv. P = HP; As colunas P são autovetores de A:" << std::endl;
  linalg::printMatrix(matrixP);

  for (size_t i = 0; i < matrixANew.size(); i++)
  {
    lamb[i] = matrixANew[i][i];
  }

  return std::make_tuple(matrixP, lamb);
}