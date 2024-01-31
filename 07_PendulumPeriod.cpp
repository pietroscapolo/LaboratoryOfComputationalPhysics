//
//  main.cpp
//  Pendolo semplice
//
//  Created by Pietro Scapolo on 04/11/22.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

//Accelerazione di gravità
double g=9.806;
//lunghezza pendolo
double l;

double f(double x){
    return g/l*sin(x);
}

int main (){
    
    vector<double> x;
    vector<double> v;
 
    //Valori iniziali
    double dt;
    cout<<"Inserisci precisione"<<endl;
    cin>>dt;
    cout<<"Inserisci lunghezza del pendolo"<<endl;
    cin>>l;
    double theta0;
    cout<<"Inserisci angolo iniziale"<<endl;
    cin>>theta0;
    x.push_back(theta0);
    v.push_back(0);
    

    
    //Variabili per Runge-Kutta
    double k1, k2, k3, k4, j1, j2, j3, j4;

    
    //CALCOLO VETTORI SPAZIO E VELOCITA
    /*
    ofstream fout("Data.txt");
    int n=10000;
    for(int i=0;i<n;i++){
        k1= v[i];
        j1= x[i];
        k2= v[i]-(dt/2*f(j1));
        j2= (x[i]+dt/2*k1);
        k3= v[i]-(dt/2*f(j2));
        j3= (x[i]+(dt/2*k2));
        k4= v[i]-(dt*f(j3));
        j4= (x[i]+dt*(k3));
        
        x[i+1]= x[i] + dt/6.*(k1+2*k2+2*k3+k4);
        v[i+1]= v[i] + dt/6.*(-f(j1)-2*f(j2)-2*f(j3)-f(j4));
        
        fout << i*dt <<"   "<<x[i]<<endl;
        }
    fout.close();
     */
    
    //CALCOLO PERIODO
    int i=-1;
    do{
        i=i+1;
        k1= v[i];
        j1= x[i];
        k2= v[i]-(dt/2*f(j1));
        j2= (x[i]+dt/2*k1);
        k3= v[i]-(dt/2*f(j2));
        j3= (x[i]+(dt/2*k2));
        k4= v[i]-(dt*f(j3));
        j4= (x[i]+dt*(k3));
        
        x[i+1]= x[i] + dt/6.*(k1+2*k2+2*k3+k4);
        v[i+1]= v[i] + dt/6.*(-f(j1)-2*f(j2)-2*f(j3)-f(j4));
    }while(x[i]*x[i+1]>0);
    
    //CALCOLO PERIODO TRAMITE INTERPOLAZIONE
    double t = (-x[i]/(x[i+1]-x[i]) +i)*4*dt;
    cout<<"il periodo è: "<<t<<"    "<<2*M_PI *sqrt(l/g)<<endl;
    //cout<<t/(2*M_PI *sqrt(l/g))<<endl;
    return 0;
}
