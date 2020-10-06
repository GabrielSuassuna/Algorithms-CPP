#include "../matrixEigenvalue.h"

std::tuple<double, std::vector<double>> matrixEigenvalue::displacementPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError, double displacement)
{
  double eigenvalueHat, eigenvalue;
  std::vector<double> vectorHat;

  std::vector<std::vector<double>> matrixHat = linalg::matrixSubtraction(
      matrix, linalg::scalarMatrixMultiplication(linalg::identityMatrix(matrix.size()), displacement));
  std::tie(eigenvalueHat, vectorHat) = matrixEigenvalue::inversePower(matrixHat, initialGuess, toleranceError);

  eigenvalue = eigenvalueHat + displacement;
  return std::make_tuple(eigenvalue, vectorHat);
}