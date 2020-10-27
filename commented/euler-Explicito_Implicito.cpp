#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

long double exactFormulaPVI1(long double y){
  double euler = 2.71828182845904523536;
  long double response = 2*pow(euler,(2.0/3.0)*y);
  return response;
}

long double myFormula(long double y){
  return ((2.0/3.0) * y);
}

void explicitEulerPVI1(long double timeLimit, long double timeStep ,long double initialValueY, long double initialValueStep){
  vector< long double >responses = { initialValueY };
  int responseIndexer = 0;
  long double steper = initialValueStep;

  while (steper <= timeLimit + timeStep){
    long double newResponse = responses[responseIndexer] + timeStep*myFormula(responses[responseIndexer]);

    responses.push_back(newResponse);
    responseIndexer++;

    long double exactResponse = exactFormulaPVI1(steper);
    long double relativeError = ( exactResponse - responses[responseIndexer]) / exactResponse;


    cout << "Resposta aproximada: " << newResponse << "\n";
    cout << "Resposta exata: " << exactResponse << "\n";
    cout << "Erro relativo: " << relativeError << "\n";
    cout << "---------------------------------------" << "\n";
    
    steper = steper + timeStep;
  }
  return ;
}

void implicitEulerPVI1(long double timeLimit, long double timeStep ,long double initialValueY, long double initialValueStep){
  vector< long double >responses = { initialValueY };
  int responseIndexer = 0;
  long double steper = initialValueStep;

  while( steper <= timeLimit + timeStep) {

    //formula no formato: y1 = y0 + deltaT*( CONSTANTE * y1 + RESTO )
    //abaixo deve ser escrito no formato: y1 = (y0 + deltaT * RESTO) / (1 - deltaT * CONSTANTE)
    long double newResponse = ( responses[responseIndexer] + timeStep*0 )/(1 - timeStep * (2.0/3.0));

    responses.push_back(newResponse);
    responseIndexer++;

    long double exactResponse = exactFormulaPVI1(steper);
    long double relativeError = ( exactResponse - responses[responseIndexer]) / exactResponse;


    cout << "Resposta aproximada: " << newResponse << "\n";
    cout << "Resposta exata: " << exactResponse << "\n";
    cout << "Erro relativo: " << relativeError << "\n";
    cout << "---------------------------------------" << "\n";
    
    steper = steper + timeStep;
  }

}

long double myFormulaForV(long double g, long double k, long double m, long double v){

  return (-g) - (k/m)*v;

}

long double myFormulaForY( long double v ){
  return v;
}

long double exactFormulaVPVI2(long double step){
  double euler = 2.71828182845904523536;

  long double resultV = -10 + 13*pow(euler,-step);
  
  return resultV;
}

long double exactFormulaYPVI2(long double step){
  double euler = 2.71828182845904523536;

  long double resultY = 150 - 10*(step) - 13*(pow(euler, -step) - 1);

  return resultY;
}

long double explicitEulerPVI2(long double timeLimit,
                              long double timeStep,
                              long double initialValueV,
                              long double initialValueY,
                              long double initialValueStep,
                              long double g,
                              long double k,
                              long double m){

  vector< long double >responsesY = { initialValueY };
  int responseIndexerForY = 0;

  vector< long double >responsesV = { initialValueV };
  int responseIndexerForV = 0;
  
  long double steper = initialValueStep;

  //long double passos = 1;
  //vector<long double> respostaAltura = { timeStep ,initialValueY, initialValueStep, passos };

  while( steper <= timeLimit + timeStep ) {

    //Si+1 = Si + dT*F(Si,Ti)
    long double newResponseForY = responsesY[responseIndexerForY] + timeStep*myFormulaForY(responsesV[responseIndexerForV]);

    responsesY.push_back(newResponseForY);
    responseIndexerForY++;

    /*if( newResponseForY >= respostaAltura[1] ){

      respostaAltura = {timeStep, newResponseForY ,steper, passos};

    }*/

    long double newResponseForV = responsesV[responseIndexerForV] + timeStep*myFormulaForV(g,k,m, responsesV[responseIndexerForV]);

    responsesV.push_back(newResponseForV);
    responseIndexerForV++;


    //Vamos calcular a resposta exata para V
    long double exactResponseV = exactFormulaVPVI2(steper);

    //E agora vamos calcular o erro relativo
    long double relativeErrorV = ( exactResponseV - responsesV[responseIndexerForV]) / exactResponseV;

    //Vamos calcular a resposta exata para Y
    long double exactResponseY = exactFormulaYPVI2(steper);

    //E agora vamos calcular o erro relativo
    long double relativeErrorY = ( exactResponseY - responsesY[responseIndexerForY]) / exactResponseY;

  
    cout << "Step: " << steper << "\n";
    cout << "Resposta aproximada de V: " << newResponseForV << "\n";
    cout << "Resposta aproximada de Y: " << newResponseForY << "\n";
    cout << "Resposta exata de V: " << exactResponseV << "\n";
    cout << "Resposta exata de Y: " << exactResponseY << "\n";
    cout << "Erro relativo de V: " << relativeErrorV << "\n";
    cout << "Erro relativo de Y: " << relativeErrorY << "\n";
    cout << "---------------------------------------" << "\n";
  
    steper = steper + timeStep;
    //passos++;
  }

  /*
  cout << "=================================" << "\n";
  for (int i = 0; i <= 3; i++){
    cout << respostaAltura[i] << "  |  ";
  }
  cout << "\n" << "=================================" << "\n";
  */
}

long double formulaForYInImplicitPVI2(long double m,
                                      long double k,
                                      long double timeStep,
                                      long double g,
                                      long double v,
                                      long double y) {
  
  return y + (( m*timeStep )/( m + k*timeStep)) * ( v - g*timeStep );

}

long double formulaForVInImplicitPVI2(long double m,
                                      long double k,
                                      long double timeStep,
                                      long double g,
                                      long double v){

  return ( (m)/(m + k*timeStep) ) * (v - g * timeStep);

}

long double implicitEulerPVI2(long double timeLimit,
                              long double timeStep,
                              long double initialValueV,
                              long double initialValueY,
                              long double initialValueStep,
                              long double g,
                              long double k,
                              long double m){
  
  vector< long double >responsesY = { initialValueY };
  int responseIndexerForY = 0;

  vector< long double >responsesV = { initialValueV };
  int responseIndexerForV = 0;
  
  long double steper = initialValueStep;

  long double passos = 1;

  vector<long double> respostaAltura = { timeStep ,initialValueY, initialValueStep, passos };


  while( steper <= timeLimit + timeStep ) {

    //y1 = y0 + (m*dT/(m+k*dT))*(v0 - g*dT)
    long double newResponseForY = formulaForYInImplicitPVI2(m, k, timeStep, g, responsesV[responseIndexerForV], responsesY[responseIndexerForY]);

    responsesY.push_back(newResponseForY);
    responseIndexerForY++;

    if( newResponseForY >= respostaAltura[1] ){

      respostaAltura = {timeStep, newResponseForY ,steper, passos};

    }

    //v1 = (m*dT/(m+k*dT))*(v0 - g*dT)
    long double newResponseForV = formulaForVInImplicitPVI2(m, k, timeStep, g, responsesV[responseIndexerForV]);

    responsesV.push_back(newResponseForV);
    responseIndexerForV++;

    //vamos agora calcular a resposta exata para V
    long double exactResponseV = exactFormulaVPVI2(steper);

    //E então seu erro relativo
    long double relativeErrorV = ( exactResponseV - responsesV[responseIndexerForV]) / exactResponseV;

    //vamos agora calcular a resposta exata para Y
    long double exactResponseY = exactFormulaYPVI2(steper);

    //E então seu erro relativo
    long double relativeErrorY = ( exactResponseY - responsesY[responseIndexerForY]) / exactResponseY;


    cout << "Step: " << steper << "\n";
    cout << "Resposta aproximada de V: " << newResponseForV << "\n";
    cout << "Resposta aproximada de Y: " << newResponseForY << "\n";
    cout << "Resposta exata de V: " << exactResponseV << "\n";
    cout << "Resposta exata de Y: " << exactResponseY << "\n";
    cout << "Erro relativo de V: " << relativeErrorV << "\n";
    cout << "Erro relativo de Y: " << relativeErrorY << "\n";
    cout << "---------------------------------------" << "\n";
  
    steper = steper + timeStep;
    passos++;
  }
 
  cout << "=================================" << "\n";
  for (int i = 0; i <= 3; i++){
    cout << respostaAltura[i] << "  |  ";
  }
  cout << "\n" << "=================================" << "\n";

} 

int main(){
  //implicitEulerPVI1(10, 0.5, 2,0.5);

  //As entradas são:
  //timeLimit, timeStep, initialValueV, initialValueY, initialValueStep, g, k, m
  
  //explicitEulerPVI2(16.3, 0.1, 5, 200, 0.1, 10, 0.25, 2.0);

  implicitEulerPVI2(7.8, 0.0001, 5, 200, 0.0001, 10, 0.25, 2.0);

  return 0;
}