#include <vector>
#include <iostream>
#include <math.h>
#include <tuple>

namespace boundaryValue
{
  void predictorCorrector(int n, std::vector<double> initialCondition, double deltaT, double t0, double error);
}