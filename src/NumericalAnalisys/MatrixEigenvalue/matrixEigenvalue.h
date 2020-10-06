#include <vector>
#include <iostream>
#include <tuple>
#include "../LinearAlgebra/linearAlgebra.h"
#include "../SystemEquations/systemEquations.h"

namespace matrixEigenvalue
{
  std::tuple<double, std::vector<double>> regularPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError);
  std::tuple<double, std::vector<double>> inversePower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError);
  std::tuple<double, std::vector<double>> displacementPower(std::vector<std::vector<double>> matrix, std::vector<double> initialGuess, double toleranceError, double displacement);
}; // namespace matrixEigenvalue
