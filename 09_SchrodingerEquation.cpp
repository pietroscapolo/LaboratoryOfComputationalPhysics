//
//  main.cpp
//  SchroedingerTD
//
//  Created by Paolo Umari on 05/12/19.
//  Copyright (c) 2019 unipd. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;

void solve_tridiagonal(int n, complex<double> *d, complex<double> *u, complex<double> *l,complex<double> *b, complex<double> *a ){

    complex<double>* alfa = new complex<double>[n];
    complex<double>* beta = new complex<double>[n];

    alfa[0]=-d[0]/u[0];
    beta[0]=b[0]/u[0];

    for(int i=1;i < n; i++){
        alfa[i]=(-l[i-1]/(u[i]*alfa[i-1])-d[i]/u[i]);
        beta[i]=b[i]/u[i]+l[i-1]*beta[i-1]/(u[i]*alfa[i-1]);

    }


    a[n-1]=(b[n-1]/l[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+d[n-1]/l[n-2]);

    for(int i=n-1;i>0;i--){
        a[i-1]=a[i]/alfa[i-1]-beta[i-1]/alfa[i-1];
    }


    free(alfa);
    free(beta);


}


int main(int argc, const char * argv[])
{
    double L,x0,q,sigma;
    double a,b, V0;
    double dt;
    int Nx,Nsteps,Nprint;
    complex<double> pic (0.,1.0);

    cout << "Lunghezza L\n";
    cin >> L;
    cout << "Posizione x0\n";
    cin >> x0;
    cout << "Momento q\n";
    cin >> q;
    cout << "Sigma \n";
    cin >> sigma;
    cout << "Limite potenziale a\n";
    cin >> a;
    cout << "Limite potenziale b\n";
    cin >> b;
    cout << "Valore potenziale V0\n";
    cin >> V0;
    cout << "Numero punti asse x Nx\n";
    cin >> Nx;
    cout << "Time step dt \n";
    cin >> dt;
    cout << "Numero time steps\n";
    cin >> Nsteps;
    cout << "Scrivi ogni numero passi\n";
    cin >> Nprint;


    double h=L/(Nx-1);
    int Nmat=Nx-2;
    double norm;

    complex<double> *psi0 = new complex<double>[Nmat];
    complex<double> *psi1 = new complex<double>[Nmat];
    complex<double> *d = new complex<double>[Nmat];
    complex<double> *u = new complex<double>[Nmat];
    complex<double> *l = new complex<double>[Nmat];
    complex<double> *f = new complex<double>[Nmat];
    double * V = new double [Nmat];
//SETTA PSI0 e NORMALIZZA

//NMAT E' IL NUMERO DI GRID STEPS INTERNI DOVE LA FUNZIONE D'ONDA PUO'
//ESSERE DIVERSA DA ZERO
//Nx E' IL NUMERO TOTALE DI GRID STEPS COMPRESI GLI ESTREMI

    norm=0.;
    for(int i=1;i<Nx-1;i++){
        double x=h*i;

        psi0[i-1]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        norm+=pow(abs(psi0[i-1]),2);

    }

    norm=norm*L/Nx;
    norm=sqrt(norm);
    cout << "Norma" <<norm <<'\n';
    for(int i=0;i<Nmat;i++){
        psi0[i]=psi0[i]/norm;
    }

    //SETTA MATRICE M

    for(int i=0;i<Nmat;i++){
        u[i]=1.;
        l[i]=1.;
    }

//Mi costruisco V
for(int i=0;i<Nmat;i++){
  double x=h*(i+1);
  if(x>=a && x<=b){
    V[i] = V0;
  }
  else{V[i] = 0;}
}

// QUI CODARE LA PARTE DIAGONALE DELLA MATRICE
// STARE ATTENTI AGLI ESTREMI
//usiamo unita di misura m=1 htagliato = 1

double m=1.;
double hcut = 1.;
for(int i=0;i<Nmat;i++){
  d[i] = pic*4.*m*pow(h,2)/(hcut*dt) -2. -2.*m*pow(h,2)*V[i]/pow(hcut,2);
}



    ofstream fileg;
    string nomeg;
    nomeg = (string) "psi_tutta.dat";
    fileg.open(nomeg,ios::out);
    fileg.precision(10);

    //LOOP SU PASSI
    for(int n=0;n<Nsteps;n++){


// QUI CODARE L'UPDATE DELLA FUNZIONE D'ONDA
//CHE VIENE MESSA IN psi1

f[0] = -psi0[1]+2.*psi0[0]+(pic*4.*m*pow(h,2)/(hcut*dt))*psi0[0] + (2.*m*pow(h,2)*V[0]/pow(hcut,2)) * psi0[0];
f[Nmat-1] = 2.*psi0[Nmat-1]-psi0[Nmat-2]+ (pic*4.*m*pow(h,2)/(hcut*dt))*psi0[Nmat-1];
for(int i=1;i<Nmat-1;i++){
  f[i] = -psi0[i+1]+2.*psi0[i]-psi0[i-1]+(pic*4.*m*pow(h,2)/(hcut*dt))*psi0[i] + (2.*m*pow(h,2)*V[i]/pow(hcut,2)) * psi0[i];
}

solve_tridiagonal(Nmat, d, u, l, f, psi1);

        if(n%Nprint==0){
            for(int i=0;i<Nmat;i++){
                double x=h*(i+1);
                fileg << x << "  " << pow(abs(psi0[i]),2) << '\n';

            }

            fileg << '\n' << '\n';

        }


        norm=0.;
        for(int i=0;i<Nmat;i++){
            psi0[i]=psi1[i];
            norm+=pow(abs(psi0[i]),2);
        }
        norm=norm*L/Nx;
      //  cout << "Passo " << n << "Norma " << norm << '\n';

    }

    fileg.close();
    free(psi0);
    free(psi1);
    free(d);
    free(u);
    free(l);
    free(f);

    return 0;
}


//./tdschroedinger.exe < parametri.txt SERVE per autocompilare programma da un file con stdin
