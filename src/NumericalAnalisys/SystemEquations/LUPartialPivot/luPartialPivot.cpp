#include "../systemEquations.h"

void systemEquations::choosePivot(std::vector<std::vector<double>> A, int k, double &pv, int &r)
{
  double pivot = fabs(A[k][k]);
  pv = A[k][k];
  r = k;

  for (uint i = k + 1; i < A.size(); i++)
  {
    if (fabs(A[i][k]) > pivot)
    {
      pivot = fabs(A[i][k]);
      pv = A[i][k];
      r = i;
    }
  }
};

void systemEquations::fakePermute(std::vector<std::vector<double>> &A, int k, int r)
{
  for (uint j = 0; j < A.size(); j++)
  {
    std::swap(A[k][j], A[r][j]);
  }
}

void systemEquations::permute(std::vector<std::vector<double>> &A, std::vector<int> &p, int k, int r)
{
  std::swap(p[k], p[r]);
  for (uint j = 0; j < A.size(); j++)
    std::swap(A[k][j], A[r][j]);
}

std::vector<double> systemEquations::retroativeIterations(std::vector<std::vector<double>> A, std::vector<double> b)
{
  std::vector<double> x;
  for (uint i = 0; i < b.size(); i++)
  {
    x.push_back(0);
  }
  x[b.size() - 1] = b[A.size() - 1] / A[A.size() - 1][A.size() - 1];
  double sum;
  for (int i = A.size() - 2; i >= 0; i--)
  {
    sum = 0;
    for (uint j = i + 1; j < A.size(); j++)
    {
      sum += A[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / A[i][i];
  }

  return x;
}

std::vector<double> systemEquations::sucessiveIterations(std::vector<std::vector<double>> A, std::vector<double> b)
{
  std::vector<double> x;
  for (uint i = 0; i < b.size(); i++)
    x.push_back(0);

  double sum;

  for (uint i = 0; i < b.size(); i++)
  {
    sum = 0;
    for (uint j = 0; j < i; j++)
    {
      sum += A[i][j] * x[j];
    }
    x[i] = b[i] - sum;
  }

  return x;
}

std::vector<double> systemEquations::luPartialPivoting(std::vector<std::vector<double>> A, std::vector<double> b)
{
  std::vector<int> p;
  std::vector<std::vector<double>> L, U;
  std::vector<double> blin;
  int r;

  for (uint i = 0; i < A.size(); i++)
  {
    L.push_back({});
    U.push_back({});
    for (uint j = 0; j < A.size(); j++)
    {
      L[i].push_back(0);
      U[i].push_back(0);
    }
  }

  for (uint i = 0; i < b.size(); i++)
    p.push_back(i);

  for (int k = 0; k < int(b.size()); k++)
  {
    double pivot;
    int r;
    choosePivot(A, k, pivot, r);

    if (pivot == 0)
      return {};

    if (k != r)
    {
      permute(A, p, k, r);
      fakePermute(L, k, r);
    }
    double m;
    for (int i = 0; i < k + 1; i++)
    {
      U[i][k] = A[i][k];
    }

    for (uint i = k + 1; i < b.size(); i++)
    {
      m = A[i][k] / A[k][k];
      A[i][k] = m;
      L[i][k] = m;
      for (uint j = k + 1; j < b.size(); j++)
      {
        A[i][j] = A[i][j] - m * A[k][j];
      }
    }
  }

  for (uint i = 0; i < b.size(); i++)
  {
    r = p[i];
    blin.push_back(b[r]);
  }

  std::vector<double> y = sucessiveIterations(A, blin);
  std::vector<double> x = retroativeIterations(A, y);
  return x;
}