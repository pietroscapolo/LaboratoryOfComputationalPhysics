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

//funzione punto medio
double puntomedio(double sup, double inf) {
    double x=(sup+inf)/2.;
    return (x);
}

//funzione calcolo funzione nel punto medio
double funzione(double sup, double inf) {
    double f = cos((sup+inf)/2.) - (sup+inf)/2.;
    return(f);
}

int main(){
    double inf=0;
    double sup=1;
    double epsilon;
    
    //intro
    cout<<"Questo algoritmo calcola la radice della funzione cosx=x con il metodo della bisezione."<<endl<<"Che precisione vuoi?"<<endl;
    cin>>epsilon;
    
    //metodo bisezione
    do{
        if((funzione(sup,sup)*funzione(inf,sup))>0){
            sup=puntomedio(inf,sup);
        } else if ((funzione(sup,sup)*funzione(inf,sup))<0){
            inf=puntomedio(inf,sup);
        }
        cout<<puntomedio(inf,sup)<<endl;
    }while(abs(funzione(inf,sup))>epsilon);
    
    return 0;
}
