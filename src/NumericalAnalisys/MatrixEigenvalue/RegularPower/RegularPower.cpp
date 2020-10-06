#include "../matrixEigenvalue.h"

std::tuple<double, std::vector<double>> matrixEigenvalue::regularPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError)
{
  double eigenvalueNew = 0;
  double eigenvalueOld;
  std::vector<double> vectorNew = initialGuess;
  std::vector<double> vectorOld;

  do
  {
    eigenvalueOld = eigenvalueNew;
    vectorOld = vectorNew;
    vectorOld = linalg::normalize(vectorOld);
    vectorNew = linalg::gaxpy(matrix, vectorOld);
    eigenvalueNew = linalg::dotProduct(vectorOld, vectorNew);
  } while (abs((eigenvalueNew - eigenvalueOld) / eigenvalueNew) > toleranceError);

  return std::make_tuple(eigenvalueNew, vectorOld);
}