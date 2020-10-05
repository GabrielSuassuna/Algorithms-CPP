#include "../matrixEigenvalue.h"

void matrixEigenvalue::powerMethod(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError)
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

  std::cout << "Autovalor Dominante: " << eigenvalueNew << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(vectorOld);
  std::cout << std::endl;
}