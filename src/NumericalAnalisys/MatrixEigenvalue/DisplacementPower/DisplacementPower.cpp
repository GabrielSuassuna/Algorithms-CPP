#include "../matrixEigenvalue.h"

std::tuple<double, std::vector<double>> matrixEigenvalue::displacementPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError, double displacement)
{
  // O método da potência com deslocamento acha o par (lambda_i, x_i) correspondente
  // ao lambda_i que estiver mais próximo de u e seja diferente de seus vizinhos.
  double eigenvalueHat, eigenvalue;
  std::vector<double> vectorHat;

  // Â = [A − uI]
  std::vector<std::vector<double>> matrixHat = linalg::matrixSubtraction(
      matrix, linalg::scalarMatrixMultiplication(linalg::identityMatrix(matrix.size()), displacement));

  // Aplicando o método da Potência Inverso sobre a matriz Â, vamos encontrar o
  // autovalor lambda_hat de menor valor absoluto, e isso significa que
  // |lambda_i − u| é o menor possível.
  std::tie(eigenvalueHat, vectorHat) = matrixEigenvalue::inversePower(matrixHat, initialGuess, toleranceError);

  // O autovalor lambda_i da matriz A é simplimente lambda_hat + u
  eigenvalue = eigenvalueHat + displacement;
  return std::make_tuple(eigenvalue, vectorHat);
}