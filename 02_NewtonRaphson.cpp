//
//  main.cpp
//  integrazione
//
//  Created by Pietro Scapolo on 21/10/22.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

//calcolo derivata
double dx(double x){
    double h= -sin(x)-1;
    return h;
}

//calcolo funzione
double funzione(double x){
    double f=cos(x)-x;
    return f;
}


int main(){

    double epsilon;
    //input precisione
    cout<<"Questo algoritmo calcola la radice della funzione: cosx=x con il metodo delle tangenti di Newton-Raphson."<<endl<<"Inserisci la precisione"<<endl;
    cin>>epsilon;
    
    //punto iniziale
    double x0=0;
    
    do {
        double xi= x0 - funzione(x0)/dx(x0);
        x0=xi;
    } while (abs(funzione(x0))>epsilon);
    
    cout<<"La radice è: "<<x0<<endl;
    
    return 0;
}
