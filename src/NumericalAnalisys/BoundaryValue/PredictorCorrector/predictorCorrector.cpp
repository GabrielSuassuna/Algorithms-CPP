#include "../boundaryValue.h"

std::vector<double> function(std::vector<double> s)
{
  std::vector<double> result = {
      -10 - (0.25 / 2) * s[0],
      s[0],
  };
  return result;
}

void boundaryValue::predictorCorrector(int n, std::vector<double> initialCondition, double deltaT, double t0, double error)
{
  std::vector<double> f1, f2, f3, f4, sBarra(2, 0.0), fBarra(2, 0.0), sOld(2, 0.0), sNew(2, 0.0);
  std::vector<std::vector<double>> s(n + 1, std::vector<double>(2, 0.0));
  double
      maxY = initialCondition[1],
      maxT = initialCondition[0],
      maxI = 0;

  s[0] = initialCondition;
  // Usar Runge Kutta para inicialização
  for (int i = 0; i < 3; i++)
  {
    f1 = function(s[i]);

    f2 = function({
        s[i][0] + (deltaT / 2.0) * f1[0],
        s[i][1] + (deltaT / 2.0) * f1[1],
    });

    f3 = function({
        s[i][0] + (deltaT / 2.0) * f2[0],
        s[i][1] + (deltaT / 2.0) * f2[1],
    });

    f4 = function({
        s[i][0] + deltaT * f3[0],
        s[i][1] + deltaT * f3[1],
    });

    s[i + 1] = {s[i][0] + (deltaT / 6) * (f1[0] + 2 * f2[0] + 2 * f3[0] + f4[0]),
                s[i][1] + (deltaT / 6) * (f1[1] + 2 * f2[1] + 2 * f3[1] + f4[1])};

    if (std::abs(s[i + 1][1]) < 1)
    {

      std::cout << "Passos relevantes quanto à queda no mar: " << std::endl;
      std::cout << "t = " << t0 + i * deltaT + deltaT << std::endl;
      std::cout << "v_" << i + 1 << " : " << s[i + 1][0] << std::endl;
      std::cout << "y_" << i + 1 << " : " << s[i + 1][1] << std::endl
                << std::endl;
    }

    if (maxY < s[i + 1][1])
    {
      maxY = s[i + 1][1];
      maxT = t0 + i * deltaT + deltaT;
      maxI = i + 1;
    }
  }

  for (int i = 3; i < n; i++)
  {
    sBarra[0] = s[i][0] + (deltaT / 24.0) * (55 * function(s[i])[0] - 59 * function(s[i - 1])[0] + 37 * function(s[i - 2])[0] - 9 * function(s[i - 3])[0]);
    sBarra[1] = s[i][1] + (deltaT / 24.0) * (55 * function(s[i])[1] - 59 * function(s[i - 1])[1] + 37 * function(s[i - 2])[1] - 9 * function(s[i - 3])[1]);

    fBarra = function(sBarra);

    sNew[0] = s[i][0] + (deltaT / 24.0) * (9 * fBarra[0] + 19 * function(s[i])[0] - 5 * function(s[i - 1])[0] + function(s[i - 2])[0]);
    sNew[1] = s[i][1] + (deltaT / 24.0) * (9 * fBarra[1] + 19 * function(s[i])[1] - 5 * function(s[i - 1])[1] + function(s[i - 2])[1]);

    sOld = {-976, -976};

    do
    {
      if (sNew[0] < -976 || sNew[1] < -976)
        sOld = sNew;
      sNew[0] = s[i][0] + (deltaT / 24.0) * (9 * fBarra[0] + 19 * function(s[i])[0] - 5 * function(s[i - 1])[0] + function(s[i - 2])[0]);
      sNew[1] = s[i][1] + (deltaT / 24.0) * (9 * fBarra[1] + 19 * function(s[i])[1] - 5 * function(s[i - 1])[1] + function(s[i - 2])[1]);
    } while (abs(abs(sOld[0] - sNew[0]) / sNew[0]) < error && abs(abs(sOld[1] - sNew[1]) / sNew[1]) < error);

    s[i + 1] = sNew;
    if (std::abs(s[i + 1][1]) < 1)
    {

      std::cout << "Passos relevantes quanto à queda no mar: " << std::endl;
      std::cout << "t = " << t0 + i * deltaT + deltaT << std::endl;
      std::cout << "v_" << i + 1 << " : " << s[i + 1][0] << std::endl;
      std::cout << "y_" << i + 1 << " : " << s[i + 1][1] << std::endl
                << std::endl;
    }

    if (maxY < s[i + 1][1])
    {
      maxY = s[i + 1][1];
      maxT = t0 + i * deltaT + deltaT;
      maxI = i + 1;
    }
  }
  std::cout << "Passo da altura máxima na trajetória: " << std::endl;
  std::cout << "t = " << maxT << std::endl;
  std::cout << "y_" << maxI << " : " << maxY << std::endl
            << std::endl;
}
