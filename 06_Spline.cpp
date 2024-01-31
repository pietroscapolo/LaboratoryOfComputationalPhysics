//
//  main.cpp
//  Spline
//
//  Created by Pietro Scapolo on 28/10/22.
//
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

int main (){
    
    double inf=0, sup=2*M_PI;
    int N=5;
    double h=(sup-inf)/N;
    
    double * b=new double [N];
    
    for (int i=0; i<N; i++){
        b[i] = 6/h*(sin((sup-inf)/N*(i+1)) - 2*sin((sup-inf)/N*i) + sin((sup-inf)/N*(i-1)));
    }
    
    double * w=new double [N];
    double * t=new double [N];
    double * v=new double [N];

    t[0]=0;
    v[0]=0;
    v[1]=h;
    w[1]=4*h;
    
    for (int i=1; i<N; i++){
        v[i]=h;
        w[i]=4*h-t[i-1]*v[i-1];
        t[i]=h/w[i];
    }
    
    double * y=new double [N];
    y[1]=b[1]/w[1];
    
     for (int i=2; i<N; i++){
         y[i]=(b[i]-v[i-1]*y[i-1])/w[i];
     }
    
    double * z=new double [N];
     z[N-1]=y[N-1];
    for (int i=N-2; i>0; i--){
        z[i]=y[i]-t[i]*z[i+1];
    }
    
    ofstream fout("dati.txt");
    for (int i=0; i<1001; i++){
        for (int j=0; j<N+1;j++){
            if ((sup-inf)/1000*i>(sup-inf)/N*j && (sup-inf)/1000*i<(sup-inf)/N*(j+1)){
                
//fout<<(e-s)/1000*i<<"\t"<<z[j+1]/(6*h)*pow((e-s)/1000*i-(e-s)/N*j,3)-z[j]/(6*h)*pow((e-s)/1000*i-(e-s)/N*(j+1),3)+(-z[j+1]*h/6+sin((e-s)/N*(j+1))/h)*((e-s)/1000*i-(e-s)/N*j)+(z[j]*h/6-sin((e-s)/N*j)/h)*((e-s)/1000*i-(e-s)/N*(j+1))<<endl;
                
    //output polinomio
fout<<(sup-inf)/1000*i<<"\t"<<sin((sup-inf)/1000*i)-(z[j+1]/(6*h)*pow((sup-inf)/1000*i-(sup-inf)/N*j,3)-z[j]/(6*h)*pow((sup-inf)/1000*i-(sup-inf)/N*(j+1),3)+(-z[j+1]*h/6+sin((sup-inf)/N*(j+1))/h)*((sup-inf)/1000*i-(sup-inf)/N*j)+(z[j]*h/6-sin((sup-inf)/N*j)/h)*((sup-inf)/1000*i-(sup-inf)/N*(j+1)))<<endl;
            }
            
        }
        
    }
    return 0;
}

