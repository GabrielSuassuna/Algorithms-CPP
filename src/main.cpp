#include "./NumericalAnalisys/MatrixEigenvalue/matrixEigenvalue.h"

int main(int argc, char const *argv[])
{
  double eigenvalue;
  std::vector<double> result;
  std::vector<std::vector<double>> matrixABar, matrixH;
  std::vector<double> guess = {1, 1, 1, 1, 1};
  std::vector<std::vector<double>> A = {
      {40, 8, 4, 2, 1},
      {8, 30, 12, 6, 2},
      {4, 12, 20, 1, 2},
      {2, 6, 1, 25, 4},
      {1, 2, 2, 4, 5},
  };
  ;
  std::tie(matrixABar, matrixH) = matrixEigenvalue::householder(A);

  std::cout << "----------------------------------------" << std::endl
            << std::endl;
  std::cout << "Matriz tridiagonal: " << std::endl;
  linalg::printMatrix(matrixABar);
  std::cout << "Matriz acumulada: " << std::endl;
  linalg::printMatrix(matrixH);
  std::cout << "----------------------------------------" << std::endl
            << std::endl;

  std::tie(eigenvalue, result) = matrixEigenvalue::regularPower(matrixABar, guess, 0.000001);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(linalg::gaxpy(matrixH, result));
  std::cout << std::endl;

  std::tie(eigenvalue, result) = matrixEigenvalue::displacementPower(matrixABar, guess, 0.000001, 31);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(linalg::gaxpy(matrixH, result));
  std::cout << std::endl;

  std::tie(eigenvalue, result) = matrixEigenvalue::displacementPower(matrixABar, guess, 0.000001, 23);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(linalg::gaxpy(matrixH, result));
  std::cout << std::endl;

  std::tie(eigenvalue, result) = matrixEigenvalue::displacementPower(matrixABar, guess, 0.000001, 11);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(linalg::gaxpy(matrixH, result));
  std::cout << std::endl;

  std::tie(eigenvalue, result) = matrixEigenvalue::inversePower(matrixABar, guess, 0.000001);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(linalg::gaxpy(matrixH, result));
  std::cout << std::endl;

  return 0;
}
