#include "./NumericalAnalisys/MatrixEigenvalue/matrixEigenvalue.h"

int main(int argc, char const *argv[])
{
  std::vector<std::vector<double>> a = {{9, 4}, {4, 3}};
  std::vector<double> b = {1, 1};
  double eigenvalue;
  std::vector<double> c;
  std::tie(eigenvalue, c) = matrixEigenvalue::regularPower(a, b, 0.000001);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(c);
  std::cout << std::endl;

  std::tie(eigenvalue, c) = matrixEigenvalue::inversePower(a, b, 0.000001);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(c);
  std::cout << std::endl;

  std::tie(eigenvalue, c) = matrixEigenvalue::displacementPower(a, b, 0.00001, -5);

  std::cout << "Autovalor Dominante: " << eigenvalue << std::endl;
  std::cout << "Autovetor Correspondente: ";
  linalg::printVector(c);
  std::cout << std::endl;
  return 0;
}
