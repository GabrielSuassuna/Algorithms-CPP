#include "./numericalIntegration.h"

// Closed approach
double numericalIntegration::trapezoid(double initialPoint,
                                       double finalPoint,
                                       double function(double))
{
  double deltaX = finalPoint - initialPoint;

  return (deltaX / 2) * (function(initialPoint) +
                         function(initialPoint + deltaX));
}

double numericalIntegration::simpsonFormula1_3(double initialPoint,
                                               double finalPoint,
                                               double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 2;

  return (h / 3) * (function(initialPoint) +
                    4 * function(initialPoint + h) +
                    function(finalPoint));
}

double numericalIntegration::simpsonFormula3_8(double initialPoint,
                                               double finalPoint,
                                               double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 3;

  return ((3 * h) / 8) * (function(initialPoint) +
                          3 * function(initialPoint + h) +
                          3 * function(initialPoint + 2 * h) +
                          function(finalPoint));
}

double numericalIntegration::substituition4Closed(double initialPoint,
                                                  double finalPoint,
                                                  double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 4;

  return ((2 * h) / 45) * (7 * function(initialPoint) +
                           32 * function(initialPoint + h) +
                           12 * function(initialPoint + 2 * h) +
                           32 * function(initialPoint + 3 * h) +
                           7 * function(finalPoint));
}

// Open approach
double numericalIntegration::openTrapezoid(double initialPoint,
                                           double finalPoint,
                                           double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 3;

  return (deltaX / 2) * (function(initialPoint + h) +
                         function(initialPoint + 2 * h));
}

double numericalIntegration::milneRule(double initialPoint,
                                       double finalPoint,
                                       double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 4;

  return ((4 * h) / 3) * (2 * function(initialPoint + h) -
                          function(initialPoint + 2 * h) +
                          2 * function(initialPoint + 3 * h));
}

double numericalIntegration::substituition3Open(double initialPoint,
                                                double finalPoint,
                                                double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 5;

  return ((5 * h) / 24) * (11 * function(initialPoint + h) +
                           function(initialPoint + 2 * h) +
                           function(initialPoint + 3 * h) +
                           11 * function(initialPoint + 4 * h));
}

double numericalIntegration::substituition4Open(double initialPoint,
                                                double finalPoint,
                                                double function(double))
{
  double
      deltaX = (finalPoint - initialPoint),
      h = deltaX / 6;
  return ((3 * h) / 10) * (11 * function(initialPoint + h) -
                           14 * function(initialPoint + 2 * h) +
                           26 * function(initialPoint + 3 * h) -
                           14 * function(initialPoint + 4 * h) +
                           11 * function(initialPoint + 5 * h));
}