#include "./NumericalMethods/BoundaryValue/boundaryValue.h"

int main(int argc, char const *argv[])
{
  std::cout << "Alunos: Gabriel Suassuna Almeida - 412715 e\n";
  std::cout << "        Pedro Victor de Oliveira Carvalho - 417338\n";

  boundaryValue::predictorCorrector(77907, {5, 200}, 0.0001, 0, 0.000001);

  return 0;
}
