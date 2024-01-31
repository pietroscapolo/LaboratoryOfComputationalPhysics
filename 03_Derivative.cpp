//
//  main.cpp
//  derivazione
//
//  Created by Pietro Scapolo on 14/10/22.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

int main(){
    int n= 100;
    double xn = 2*M_PI;
    double h= xn / n;
    
    double *x = new double [n];
    double *y = new double [n];
    double *dx = new double [n];
    double *dx2 = new double [n];
    double *dx4 = new double [n];

    
    ofstream filedx("funzione_errore.txt");
    
    
    //calcolo seno e coseno
    for(int i=0;i<n;i++){
        x[i]= sin(i*h);
        y[i]= cos(i*h);
    }
    
    
    //input quale metodo
    int metodo;
    cout<<"Ciao! Scegli il metodo di derivazione:"<<endl<<"1. O(h)"<<endl<<"2. O(h^2)"<<endl<<"3. O(h^4)"<<endl;
    cin>>metodo;
    
    //calcoli
    switch(metodo){
        case 1:
            //metodo O(h)
            for(int i=0;i<n-1;i++){
                dx[i]= (x[i+1]-x[i])/h;
                filedx << i*h <<"    "<< dx[i]-y[i]<<endl;
            }
            break;
            
        case 2:
            //metodo O(h^2)
            for(int i=1;i<n-1;i++){
                dx2[i]= (x[i+1]-x[i-1])/(2*h);
                filedx << i*h <<"    "<< dx2[i]-y[i]<<endl;
            }
            break;
            
        case 3:
            //metodo O(h^4)
            for(int i=2;i<n-2;i++){
                dx4[i]= (-x[i+2] +(8*x[i+1])-(8*x[i-1])+x[i-2])/(12*h);
                filedx << i*h <<"    "<< dx4[i]-y[i]<<endl;
            }
            break;
            
        default:
            filedx.close();
    }
    
    delete  [] x;
    return 0;
}
