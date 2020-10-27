#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

long double myFormulaForV(long double g, long double k, long double m, long double v){

  return (-g) - (k/m)*v;

}

long double myFormulaForY( long double v ){
  return v;
}

long double calculateWithF(long double F1, long double F2, long double F3){
    return ((1.0/6.0)*F1 + (4.0/6.0)*F2 + (1.0/6.0)*F3);
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

void rungeKuttaTercOrdemPVI2(long double timeLimit,
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

    while (steper <= timeLimit + timeStep){
        //Precisamos de alguns ingredientes para a fórmula final, S0, Si+1/2, barra Si+1
        //barrado Si+1/2 = S0 + dT/2*F(S0,t0)
        long double vBarra1 = responsesV[responseIndexerForV] + (timeStep/2.0)*myFormulaForV(g, k, m, responsesV[responseIndexerForV]);
        long double yBarra1 = responsesY[responseIndexerForY] + (timeStep/2.0)*myFormulaForY(responsesV[responseIndexerForV]);

        //barrado Si+1 = S0 + dT*F(S0,T0)
        long double vBarra2 = responsesV[responseIndexerForV] + timeStep*myFormulaForV(g, k, m, responsesV[responseIndexerForV]);
        long double yBarra2 = responsesY[responseIndexerForY] + timeStep*myFormulaForY(responsesV[responseIndexerForV]);

        //Si+1 = S0 + dT[1/6*F(S0,T0) + 4/6*F(Si1/2) + 1/6*F(S1,T1)]
        long double newResponseForV = responsesV[responseIndexerForV] + timeStep*((1.0/6.0)*myFormulaForV(g, k, m, responsesV[responseIndexerForV]) + (4.0/6.0)*myFormulaForV(g, k, m, vBarra1) + (1.0/6.0)*myFormulaForV(g, k, m, vBarra2));
        long double newResponseForY = responsesY[responseIndexerForY] + timeStep*((1.0/6.0)*myFormulaForY(responsesV[responseIndexerForV]) + (4.0/6.0)*myFormulaForY(vBarra1) + (1.0/6.0)*myFormulaForY(vBarra2));

        //Atualizamos a ultima resposta calculada
        responsesV.push_back(newResponseForV);
        responseIndexerForV++;
        
        responsesY.push_back(newResponseForY);
        responseIndexerForY++;

        //Vamos agora calcular a resposta exata para V
        long double exactResponseV = exactFormulaVPVI2(steper);

        //E então calcular nosso erro relativo para V
        long double relativeErrorV = ( exactResponseV - responsesV[responseIndexerForV]) / exactResponseV;

        //Vamos agora calcular a resposta exata para Y
        long double exactResponseY = exactFormulaYPVI2(steper);

        //E então calcular nosso erro relativo para Y
        long double relativeErrorY = ( exactResponseY - responsesY[responseIndexerForY]) / exactResponseY;

        //Para termos uma visualização melhor, vamos armazenar a maior altura
        if( newResponseForY >= respostaAltura[1] ){
            respostaAltura = {timeStep, newResponseForY ,steper, passos};
        }
        
        cout << "Step: " << steper << "\n";
        cout << "Resposta aproximada de V: " << newResponseForV << "\n";
        cout << "Resposta aproximada de Y: " << newResponseForY << "\n";
        cout << "Resposta exata de V: " << exactResponseV << "\n";
        cout << "Resposta exata de Y: " << exactResponseY << "\n";
        cout << "Erro relativo de V: " << relativeErrorV << "\n";
        cout << "Erro relativo de Y: " << relativeErrorY << "\n";
        cout << "---------------------------------------------" << "\n";
        
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
  //Esses são os valores a serem inseridos
  //timeLimit, timeStep, initialValueV, initialValueY, initialValueStep, g, k, m
  rungeKuttaTercOrdemPVI2(7.8, 0.0001, 5, 200, 0.0001, 10, 0.25, 2.0);

  return 0;
}