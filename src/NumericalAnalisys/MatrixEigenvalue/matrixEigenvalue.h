#include <vector>
#include <iostream>
#include <tuple>
#include <math.h>
#include "../LinearAlgebra/linearAlgebra.h"
#include "../SystemEquations/systemEquations.h"

namespace matrixEigenvalue
{
  std::tuple<double, std::vector<double>> regularPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError);
  std::tuple<double, std::vector<double>> inversePower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError);
  std::tuple<double, std::vector<double>> displacementPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError, double displacement);
  std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> householder(std::vector<std::vector<double>> matrix);
  std::tuple<std::vector<std::vector<double>>, std::vector<double>> jacobiMethod(std::vector<std::vector<double>> matrixA, double toleranceError);
}; // namespace matrixEigenvalue
