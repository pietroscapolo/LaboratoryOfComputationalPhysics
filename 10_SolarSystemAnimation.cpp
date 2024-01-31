//
//  main.cpp
//  Sistema solare
//
//  Created by Pietro Scapolo on 11/11/22.
//
#include <fstream>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <filesystem>
#include <chrono>

using namespace std;

//struttura pianeta
struct p {
    double m;
    double r[3], v[3];
};

//costante di gravità
const double G=6.674e-11;

//lettura dati
p leggi(ifstream &file)
{
    //importazione dati
    p p;
    file >> p.m;
    file >> p.r[0] >> p.r[1] >> p.r[2];
    file >> p.v[0] >> p.v[1] >> p.v[2];
    
    //da km a m
    for (int i=0; i<3; i++){
        p.r[i] = p.r[i]*1000;
        p.v[i] = p.v[i]*1000;
    }
    return p;
}

//ALLOCAZIONE PIANETI IN VETTORE SISTEMA
vector<p> importa(ifstream &file){
    int n;
    file >> n;
    vector<p> s;
    //per ogni pianeta, riempe la struttura e la mette nel vettore s (sistema)
    for(int i=0; i<n; i++){
        p p = leggi(file);
        s.push_back(p);
    }
    return s;
}

//COPIA S
vector<p> vecpmatch(vector<p> &s)
{
    vector<p> s1;
    for (int i=0; i< s.size(); i++) s1.push_back(s[i]);
    return s1;
}

//CENTRO DI MASSA
p CM(vector<p> &s)
{
    p cm;
    //massa totale
    for (int i=0; i< s.size(); i++) cm.m += s[i].m;
    //parametri pesati
    for (int i=0; i< s.size(); i++){
        for (int k=0; k<3; k++){
            cm.r[k] += (s[i].m * s[i].r[k]);
            cm.v[k] += (s[i].m * s[i].v[k]);
        }
    }
    //rapporto con massa totale
    for (int k=0; k<3; k++){
        cm.r[k] = cm.r[k]/cm.m;
        cm.v[k] = cm.v[k]/cm.m;
    }
    return cm;
}

//FUNZIONE FORZA
array<double,3> F(vector<p> &s, p &p){
    //vettore differenza coordinate
    double r[3]{0};
    //forza
    array<double,3> forza{0};
    
    //calcolo della forza
    for (int j=0; j < s.size(); j++){
        if(s[j].m !=p.m)
        {
            //calcola la differenza tra la posizione di due pianeti
            for (int k=0; k<3; k++) r[k] = s[j].r[k] -p.r[k];
            //modulo distanza
            double modr = sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));
            //componente della forza
            double f = (G * s[j].m *p.m)/pow(modr,3);
            //allocazione vettore forza
            for (int k=0; k<3; k++) forza[k] += f*r[k];
        }
    }
    return forza;
}

//ENERGIA POTENZIALE
long double U(vector<p> &s){
    double r[3];
    long double U = 0;
    for (int j = 0; j < s.size(); j++){
        for (int i = 0; i < s.size(); i++){
            if( i != j ){
                for (int k = 0; k < 3; k++) r[k] = s[j].r[k] - s[i].r[k];
                double modr = sqrt( pow(r[0],2) + pow(r[1],2) + pow(r[2],2) );
                U += -(G * s[j].m * s[i].m)/(modr)*0.5;
            }
        }
    }
   return U;
}

//ENERGIA CINETICA
long double T(vector<p> &s){
    long double K = 0;
    for (int j = 0; j < s.size(); j++){
        double modv2 = ( pow(s[j].v[0],2) + pow(s[j].v[1],2) + pow(s[j].v[2],2) );
        K += (0.5 * s[j].m )*(modv2);
    }
   return K;
}

//VERLET da controllare
vector<p> verlet(vector<p> &s, double dt){
    p cm=CM(s);
    std::vector<std::array<double,3>> forze;
    for (size_t i = 0; i < s.size(); i++)
    {
        //for (size_t k = 0; k < 3; k++) s[i].v[k] -= cm.v[k];
        std::array<double,3> f=F(s,s[i]);
        forze.push_back(f);
        for (size_t k = 0; k < 3; k++) s[i].r[k] = s[i].r[k] + s[i].v[k]*dt + (0.5/s[i].m)* f[k] *dt*dt; //r
    }
    for (size_t i = 0; i < s.size(); i++)
    {
        std::array<double,3> f1=F(s,s[i]);
        for (size_t k = 0; k < 3; k++) s[i].v[k] = s[i].v[k] + (0.5/s[i].m)*(forze[i][k] + f1[k])*dt; //v
    }
    return s;
}

//scrivi i nomi dei files
vector<string> filenames(vector<p> &s)
{
    vector<string> fnames;
    for (int i = 0; i < s.size(); i++){
        char buffer[30];
        sprintf(buffer, "object%d.dat" ,i);
        fnames.push_back(buffer);
    }
    return fnames;
}

//crea files
void writeappfiles(vector<p> &s, vector<string> &fnames, double t){
   for (int i = 0; i < s.size(); i++){
        ofstream file;
        file.open(fnames[i], ios_base::app);
        file << t << "\t" << s[i].r[0] << "\t" << s[i].r[1] << "\t" << s[i].r[2] <<"\n";
        file.close();
   }
}

int main(int argc, const char * argv[])
{
    ifstream file;
    file.open("sistema.dat", ios::in);
    vector<p> s = importa(file);
    double giorno = 8.64e4;
    double dt = giorno/1.;
    double anni = 50.;
    long double t = giorno*365.*anni; //tempo in anni
    long int n = t/dt;
    auto fnames = filenames(s);
    cout<< n <<"\n";
    vector<p> s1;
    vector<long double> Ktot,Utot;
    int printcounter = 0;
    writeappfiles (s, fnames, 0*dt);
    auto start_time = chrono::high_resolution_clock::now();
    for (long int i = 0; i < n; i++)
    {
        
        s1 = verlet(s,dt);
        if( printcounter == 1 ){
            writeappfiles (s1, fnames, i*dt);
            printcounter = 0;
            Utot.push_back(U(s1));
            Ktot.push_back(T(s1));
        }
        s = vecpmatch(s1);
        printcounter += 1;
    }

  
    //stampo energie totali
    ofstream efile;
    efile.open("Etot.dat", ios::out);
    efile.precision(18);
    for (size_t i = 0; i < Utot.size(); i++){
        efile << Ktot[i] + Utot[i] <<"\n";
    }
    efile.close();
    auto current_time = chrono::high_resolution_clock::now();
    cout << "Program has been running for " << chrono::duration_cast<chrono::seconds>(current_time - start_time).count() << " seconds" << endl;
    return 0;
}

