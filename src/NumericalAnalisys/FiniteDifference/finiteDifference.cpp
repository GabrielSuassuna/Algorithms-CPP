#include "./finiteDifference.h"

double finiteDifference::foward(double x, double deltaX, double function(double))
{
  return (function(x + deltaX) - function(x)) / deltaX;
}

double finiteDifference::backward(double x, double deltaX, double function(double))
{
  return (function(x) - function(x - deltaX)) / deltaX;
}

double finiteDifference::centered(double x, double deltaX, double function(double))
{
  return (function(x + deltaX) - function(x - deltaX)) / (2 * deltaX);
}