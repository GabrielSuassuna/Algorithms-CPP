#include "../matrixEigenvalue.h"

std::tuple<double, std::vector<double>> matrixEigenvalue::inversePower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError)
{
  double eigenvalueNew = 0;
  double eigenvalueOld, eigenvalue;
  std::vector<double> vectorNew = initialGuess;
  std::vector<double> vectorOld;

  do
  {
    eigenvalueOld = eigenvalueNew;
    vectorOld = vectorNew;
    vectorOld = linalg::normalize(vectorOld);

    // Esse passo é equivalente ao produto entre a inversa da matriz A e o vetor
    // v_k reescalonado
    vectorNew = systemEquations::luPartialPivoting(matrix, vectorOld);

    eigenvalueNew = linalg::dotProduct(vectorOld, vectorNew);
  } while (abs((eigenvalueNew - eigenvalueOld) / eigenvalueNew) > toleranceError);

  // Utilizamos o método da potência regular a partir de pequenas modicações para
  // para gerarmos um passo equivalente a multiplicação da inversa, obtendo 1/lambda_n
  // que é o autovalor dominante de A^-1. No qual, lambda_n seria o menor autovalor
  eigenvalue = 1 / eigenvalueNew;

  return std::make_tuple(eigenvalue, vectorOld);
}