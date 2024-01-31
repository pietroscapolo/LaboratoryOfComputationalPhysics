#include <fstream>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <array>

using namespace std;

int n=1000;
int p;
double L=10;
double h= L/(n-1);

double V(int i){
    double v=0;
    switch (p) {
        case 1:{ v=0;
            break;
        }
        case 2:{
            double a=3, b=7;
            if((i*h)<b && (i*h)>a){ v=2;
            } else{v=0.;
            }
            break;
        }
        case 3:{ v= pow(i*h-L/2,2);
            break;
        }
    }
    return v;
}

double numerov(double E){
      double ya = 0; //sarebbe y_1=0;
      double yb = 0.1;
      double yc=0;
      for(int i=1; i<n-1; i++){
      yc = ((10*pow(h,2)*2*(-E)+24)*yb+(2*(-E)*pow(h,2)-12)*ya)/(12-2*(-E)*pow(h,2));
      ya = yb;
      yb = yc;
      }
      return yc;
    }

void stampaonda(double E,int j){
    vector <double> y;
    y.push_back(0);
    y.push_back(0.0001);

    ofstream file;
    file.open("Stato_" + to_string(j) +".txt");

    for(int i=1; i<n; i++){
        y.push_back(((10*pow(h,2)*-(E-V(i))+24)*y[i]+(-(E-V(i))*pow(h,2)-12)*y[i-1])/(12+(E-V(i+1))*pow(h,2)));
        file << i<<"    "<<y[i]<<endl;
    }
    file.close();
}

double bisezione(double inf, double sup){
    double m;
    do{
        m =(sup+inf)/2.;
        if((numerov(inf)*numerov(m))>0){
            sup=m;
        } else if ((numerov(inf)*numerov(m))<0){
            inf=m;
        }
    }while(abs(numerov(m))<1e-9);
    return m;
}


int main() {
    double Emin=0;
    double Emax=10;
    double DeltaE=0.001;
    double e=0;


    cout<<"Ciao! Scegli il potenziale:"<<endl<<"1. Buca"<<endl<<"2. Barriera"<<endl<<"3. Armonico"<<endl;
    cin>>p;

    int j=1;
    for(double i=Emin;i<Emax; i=i+DeltaE){
        //cout<<numerov(i)<<endl;;
        if ((numerov(i)*numerov(i+DeltaE))<0){
            e=bisezione(i,i+DeltaE);
            cout<<e<<"  "<<0.5*(pow(M_PI,2)/pow(L,2)*pow(j,2))<<endl;
            //stampaonda(e,j);
            j++;
        }
    }
    return 0;
}



