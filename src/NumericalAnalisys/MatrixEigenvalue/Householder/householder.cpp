#include "../matrixEigenvalue.h"

std::vector<std::vector<double>> holseholderMethod(std::vector<std::vector<double>> matrix, int i)
{
  // Inicializamos os vetores w e wl
  std::vector<double> w(matrix.size(), 0.0);
  std::vector<double> wL(matrix.size(), 0.0);

  // Copia os elementos abaixo da diagonal da coluna i da matriz A para as
  // respectvas posições no vetor w
  for (size_t position = i + 1; position < matrix.size(); position++)
  {
    w[position] = matrix[position][i];
  }

  // Calculamos os comprimento do vetor w
  double wSize = std::sqrt(linalg::dotProduct(w, w));

  // Copiamos o wSize para a posição i + 1 do vetor wL
  wL[i + 1] = wSize;

  // Calculamos um vetor N = w - wL e normalizo-o.
  std::vector<double> n = linalg::normalize(linalg::vectorSubtraction(w, wL));

  // Montamos a matriz de Householder
  // H = I -2nn^T
  return linalg::matrixSubtraction(
      linalg::identityMatrix(matrix.size()),
      linalg::scalarMatrixMultiplication(linalg::vectorOuterProduct(n, n), 2));
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> matrixEigenvalue::householder(std::vector<std::vector<double>> matrix)
{
  std::vector<std::vector<double>> matrixHAux, matrixABar;
  std::vector<std::vector<double>> matrixH = linalg::identityMatrix(matrix.size());
  std::vector<std::vector<double>> matrixA = matrix;

  for (size_t i = 0; i < matrix.size() - 2; i++)
  {
    // Construção da matriz de Householder do passo i
    // Ele é baseado na coluna i da matriz do passo anterior
    matrixHAux = holseholderMethod(matrixA, i);

    // Transformação de similaridade do passo i
    matrixABar = linalg::matrixMultiplication(
        linalg::matrixMultiplication(
            linalg::transpose(matrixHAux), matrixA),
        matrixHAux);

    // Salvar para o próximo passo.
    matrixA = matrixABar;

    // Acumular o produto das matrizes de Householder
    matrixH = linalg::matrixMultiplication(matrixH, matrixHAux);
  }

  // No final do loop, a matriz A_bar já está no formato de uma matriz tridiagonal
  return std::make_tuple(matrixABar, matrixH);
}