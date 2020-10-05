#include "./NumericalAnalisys/MatrixEigenvalue/matrixEigenvalue.h"

int main(int argc, char const *argv[])
{
  std::vector<std::vector<double>> a = {{9, 4}, {4, 3}};
  std::vector<double> b = {1, 1};
  matrixEigenvalue::powerMethod(a, b, 0.000001);
  return 0;
}
