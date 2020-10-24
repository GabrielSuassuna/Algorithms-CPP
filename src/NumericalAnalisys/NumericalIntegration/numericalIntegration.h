namespace numericalIntegration
{
  double trapezoid(double initialPoint,
                   double finalPoint,
                   double function(double));
  double simpsonFormula1_3(double initialPoint,
                           double finalPoint,
                           double function(double));
  double simpsonFormula3_8(double initialPoint,
                           double finalPoint,
                           double function(double));
  double substituition4Closed(double initialPoint,
                              double finalPoint,
                              double function(double));
  double openTrapezoid(double initialPoint,
                       double finalPoint,
                       double function(double));
  double milneRule(double initialPoint,
                   double finalPoint,
                   double function(double));
  double substituition3Open(double initialPoint,
                            double finalPoint,
                            double function(double));
  double substituition4Open(double initialPoint,
                            double finalPoint,
                            double function(double));

} // namespace numericalIntegration