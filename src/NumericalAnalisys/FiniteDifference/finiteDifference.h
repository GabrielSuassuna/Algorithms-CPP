#include <vector>
#include <iostream>
#include <tuple>
#include <math.h>
#include "../LinearAlgebra/linearAlgebra.h"
#include "../SystemEquations/systemEquations.h"

namespace finiteDifference
{
  double foward(double x, double deltaX, double function(double));
  double backward(double x, double deltaX, double function(double));
  double centered(double x, double deltaX, double function(double));
};