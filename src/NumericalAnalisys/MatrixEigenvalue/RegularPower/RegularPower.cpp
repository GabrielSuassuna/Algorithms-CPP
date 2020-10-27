#include "../matrixEigenvalue.h"

std::tuple<double, std::vector<double>> matrixEigenvalue::regularPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError)
{
  double eigenvalueNew = 0;
  double eigenvalueOld;

  // A variável recebe o vetor v_0 qualquer (chute inicial).
  std::vector<double> vectorNew = initialGuess;
  std::vector<double> vectorOld, outputVector;

  do
  {
    eigenvalueOld = eigenvalueNew;
    vectorOld = vectorNew;

    // Para evitar que o vetor v_k cresça ou diminua muito, e sabendo que tudo que
    // nos interessa é a direção do vetor, em cada passo é feito um reescalonamento
    // do tamanho do vetor
    outputVector = linalg::normalize(vectorOld);

    // Multiplicamos o vetor v_k-1 pela matriz A para obter v_k
    vectorNew = linalg::gaxpy(matrix, outputVector);

    // O autovalor é aproximadamente igual ao produto escalar entre o vetor v_k e
    // o vetor reescalonado
    eigenvalueNew = linalg::dotProduct(outputVector, vectorNew);
  } while (abs((eigenvalueNew - eigenvalueOld) / eigenvalueNew) > toleranceError);

  return std::make_tuple(eigenvalueNew, outputVector);
}