#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> jacobiMatrixBasedOldRMatrix(std::vector<std::vector<double>> matrixA, int i, int j)
{
  // Matriz identidade com n x n elementos
  std::vector<std::vector<double>> matrixJij = linalg::identityMatrix(matrixA.size());
  double theta, error = 0.000001;

  // Considerar Aij = 0, retornar matriz identidade
  if (abs(matrixA[i][j]) <= error)
  {
    return matrixJij;
  }

  // Considerar Ajj = 0, então
  if (abs(matrixA[j][j]) <= error)
  {
    // O numerador será positivo e assumimos tangente tende a +Inf
    if (matrixA[i][j] < 0)
    {
      // PI/2 = atan(1.0) * 2
      theta = atan(1.0) * 2;
    }
    // O numerador será negativo e assumimos tangente tende a -Inf
    else
    {

      // PI/2 = atan(1.0) * 2
      theta = -atan(1.0) * 2;
    }
  }

  // Esta função já retorna um ângulo +/-
  // no primeiro quadrante sentido anti-horário (+)
  // no primeiro quadrante sentido horário (-)
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
  // Na inicialização, R velha não tem a estrutura de uma matriz triangular superior
  std::vector<std::vector<double>> matrixROld = matrixA;
  std::vector<std::vector<double>> matrixRNew, matrixJij, matrixQ, matrixR;

  // Para transformar a matriz original, precisamos zerar todos os elementos
  // abaixo do elemento da diagonal principal em cada coluna exceto na última.
  // Portanto, o loop j das colunas vai até a penúltima coluna. O loop das linhas
  // para uma dada coluna j, percorre todos os elementos abaixo da linha da
  // diagonal, ou seja, da linha (j+1) até a última linha (linha n). Como a matriz
  // é simétrica, o que fizermos para as colunas acontecerá igualmente nas linhas

  // Laço das colunas
  for (size_t j = 0; j < matrixA[0].size() - 1; j++)
  {
    // Laço das linhas
    for (size_t i = j + 1; i < matrixA.size(); i++)
    {
      // Construção da matriz de Jacobi Jij. A matriz de Jacobi é uma matriz de rotação.
      matrixJij = jacobiMatrixBasedOldRMatrix(matrixROld, i, j);
      // Matriz modificada com elemento (i,j) zerado
      matrixRNew = linalg::matrixMultiplication(matrixJij, matrixROld);
      // Salvar para o próximo passo.
      matrixROld = matrixRNew;
      // Acumular o produto das matrizes de Jacobi
      matrixQT = linalg::matrixMultiplication(matrixJij, matrixQT);
    }
  }

  matrixQ = linalg::transpose(matrixQT);

  // No final do loop externo, o formato da matriz R nova é triangular superior.
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
    // Decomposição QR (devolve as matrizes Q e R tais que A velha = QR
    // onde Q é ortogonal e R é triangular superior).
    std::tie(matrixQ, matrixR) = decompositionQR(matrixAOld);
    // Calcula a nova matriz como A nova = RQ (na ordem reversa)
    // Isso é equivalente a transformação de similaridade (Q^T A Q)
    matrixANew = linalg::matrixMultiplication(matrixR, matrixQ);
    // Salvar para a próxima iteração
    matrixAOld = matrixANew;
    std::cout << "Matriz A:" << std::endl;
    linalg::printMatrix(matrixANew);
    // Acumular o produto das matrizes Q
    matrixP = linalg::matrixMultiplication(matrixP, matrixQ);
    // Verificação se a matriz A nova é diagonal
    val = linalg::sumSquaresTermsBelowDiagonal(matrixANew);
    count++;
  } while (val > toleranceError);

  std::cout << std::endl
            << "i. Matriz diagonal:" << std::endl;
  linalg::printMatrix(matrixANew);

  // Ao sair do loop, o formato da matriz A nova já está suficientemente próximo do formato de
  // uma matriz diagonal. Assim, os elementos da diagonal são os autovalores da matriz original
  // de entrada e as colunas de P são os autovetores correspondente.

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