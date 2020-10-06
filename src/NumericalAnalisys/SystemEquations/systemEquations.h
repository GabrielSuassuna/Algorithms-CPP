#include <vector>
#include <stdexcept>
#include <math.h>
#include <iostream>

namespace systemEquations
{
  void choosePivot(std::vector<std::vector<double>> A, int k, double &pv, int &r);
  void permute(std::vector<std::vector<double>> &A, std::vector<int> &p, int k, int r);
  void fakePermute(std::vector<std::vector<double>> &A, int k, int r);
  std::vector<double> retroativeIterations(std::vector<std::vector<double>> A, std::vector<double> b);
  std::vector<double> sucessiveIterations(std::vector<std::vector<double>> A, std::vector<double> b);
  std::vector<double> luPartialPivoting(std::vector<std::vector<double>> A, std::vector<double> b);
} // namespace systemEquations