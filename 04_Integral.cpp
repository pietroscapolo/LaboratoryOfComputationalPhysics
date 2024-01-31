//
//  main.cpp
//  integrazione
//
//  Created by Pietro Scapolo on 14/10/22.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

int main(){
    double inf=0;
    double sup=1;
    double value = double(exp(sup)-sup);
  
    ofstream filer("Rettangoli.txt");
    ofstream filet("Trapezi.txt");
    ofstream files("Simpson.txt");
    
    for(int i = 0;i<1000;i++){
        
        //divisione di n in scala logaritmica
        int n=round(9+pow(999991,(double)i/999));
        
        //spessore dx
        double h=(sup-inf)/n;
    
        //dichiarazione valori integrali
        double r=0;
        double t=0;
        double s=0;
            
        //rettangoli
        for(int j=0; j<n; j++){
            r += (h*exp(j*h)); //rettangoli naif
            //r += 0.5*(exp((j+1)*h)+exp(j*h))*h; //rettangoli
        }
        
        //trapezi
        t=0.5*(exp(inf)+exp(sup))*h;
        for(int j=1; j<n; j++){
            t+=exp(j*h)*h;
        }
        
        //metodo Simpson
            s=(1./3.)*h*(exp(inf)+exp(sup));
            for(int j=1; j<n; j++){
                    s+=((j%2+1)*2)/3.*exp(j*h)*h;
                }

        //stampa
        filer << n <<"    "<< abs(r-value) <<endl;
        filet << n <<"    "<< abs(t-value) <<endl;
        if (n%2==0){
            files << n <<"    "<< abs(s-value) <<endl;
            }
        }
    
    //chiusura files
    filer.close();
    filet.close();
    files.close();
    
    return 0;
}
