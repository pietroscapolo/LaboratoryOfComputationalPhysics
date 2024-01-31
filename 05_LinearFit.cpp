//QUESTO PROGRAMMA è GIUSTO
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

int main(){
  double datix, datiy, sigmax, sigmay, sigmai, sigmatot, sigmamenouno, chi, chitot, chifin, dy, dytot;
  vector <double> datax, datay, si, yattesi, chiq, errpost, incertezzex, incertezzey, incertezzetot;
  ifstream datinputx ("/Users/pietroscapolo/Library/CloudStorage/OneDrive-UniversitàdegliStudidiPadova/C++/Interpolazione/x.txt");
  ifstream datinputy ("/Users/pietroscapolo/Library/CloudStorage/OneDrive-UniversitàdegliStudidiPadova/C++/Interpolazione/y.txt");
  ifstream sigxinput ("/Users/pietroscapolo/Library/CloudStorage/OneDrive-UniversitàdegliStudidiPadova/C++/Interpolazione/sigmax.txt");
  ifstream sigyinput ("/Users/pietroscapolo/Library/CloudStorage/OneDrive-UniversitàdegliStudidiPadova/C++/Interpolazione/sigmay.txt");

  //INSERIMENTO DATI NEI VETTORI
  while(datinputx >> datix){
    datax.push_back(datix);
  }
  //cout << datax.size() << endl;
  datinputx.close();
  while(datinputy >> datiy){
    datay.push_back(datiy);
  }
  datinputy.close();
  while(sigxinput >> sigmax){
    incertezzex.push_back(sigmax);
  }
  sigxinput.close();
  while(sigyinput >> sigmay){
    incertezzey.push_back(sigmay);
  }
  sigyinput.close();

  int x;

  cout << "\n\n\nSCEGLI IL METODO DA UTILIZZARE \nMetodo 1 = minimi quadrati (incertezze sulle x trascurabili e sulle y tutte uguali) \nMetodo 2 = incertezze sulle x trascurabili e sulle y incertezze diverse \nMetodo 3 = si considerano le incertezze sia sulle x che sulle y" << endl << endl << endl;

  cout << "Metodo 1, 2 o 3? \n";
  cin >> x;
  cout << endl;

  if(x==1){
  //CALCOLI FORMA SEMPLIFICATA
  double sommax, sommay, sommax2, sommay2, sommaxy, delta, b, a, sigmaa, sigmab;
  int i=0;
  while(i<datax.size()){
    sommax += datax.at(i);
    sommax2 += pow(datax.at(i), 2);
    sommay += datay.at(i);
    sommay2 += pow(datay.at(i), 2);
    sommaxy += datax.at(i)*datay.at(i);
     i++;
  }
  delta = datax.size()*sommax2 - sommax*sommax;
  b = (datax.size()*sommaxy-sommax*sommay)/delta;
  a = (sommax2*sommay-sommax*sommaxy)/delta;
  sigmaa = incertezzey.at(1)*sqrt(sommax2/delta);
  sigmab = incertezzey.at(1)*sqrt(datax.size()/delta);
   cout << "parametri: \n \t b: " << b << "\t\t\t± " << sigmab << "\n \t a: " << a << "\t\t\t± " << sigmaa << endl << endl;
     //CHI QUADRO
  for(int w=0; w<datax.size(); w++){
  float y;
  y=datax.at(w)*b+a;
  yattesi.push_back(y);
  //cout << yattesi.at(w) << endl;
  }
  for(int i=0; i<yattesi.size() && i<datay.size() && i<incertezzey.size(); i++){
    chi = pow((datay.at(i)-yattesi.at(i))/(incertezzey.at(i)), 2);
    chiq.push_back(chi);
    //cout << chiq.at(i) << endl;
  }
  for(int y=0; y<chiq.size(); y++){
    chitot += chiq.at(y);
  }
  //ERRORE A POSTERIORI
  for(int z=0; z<yattesi.size() && z<datay.size(); z++){
    float num= pow((datay.at(z)-yattesi.at(z)), 2);
    errpost.push_back(num);
  }
  for(int a=0; a<errpost.size(); a++){
    dytot += errpost.at(a);
  }

  dy = sqrt (dytot/(datax.size()-2));
  cout << "chi quadro \t \t \t \t \t" << chitot << endl << endl;
  cout << "errore a posteriori \t \t" << dy <<endl;
  }

  if(x==2){
   //FORMA SEMPLIFICATA CON SY =/=
   double sommax, sommay, sommax2, sommay2, sommaxy, sommaxi, sommayi, sommax2i, sommay2i, sommaxyi, delta, sigmaitot;
  long double b, a, sigmaa, sigmab;
  int i=0;
  while(i<datax.size()){
    sommax += datax.at(i)/(pow(incertezzey.at(i), 2));
    sommax2 += pow(datax.at(i), 2)/(pow(incertezzey.at(i), 2));
    sommay += datay.at(i)/(pow(incertezzey.at(i), 2));
    sommay2 += pow(datay.at(i), 2)/(pow(incertezzey.at(i), 2));
    sommaxy += datax.at(i)*datay.at(i)/(pow(incertezzey.at(i), 2));
    sigmaitot += 1/(pow(incertezzey.at(i), 2));
     i++;
  }
  delta = sigmaitot*sommax2 - sommax*sommax;
  //cout << delta << endl << endl;
  b = (sigmaitot*sommaxy-sommax*sommay)/delta;
  a = (sommax2*sommay-sommax*sommaxy)/delta;
  sigmaa = pow((sommax2/delta), 0.5);
  sigmab = pow((sigmaitot/delta), 0.5);

  cout << "SEMPLIFICATO con sy diversi \nparametri: \n \t b: " << b << "\t\t\t± " << sigmab << "\n \t a: " << a << "\t\t\t± " << sigmaa << endl << endl;
  ofstream risultati ("semplificato.txt");
  risultati << "SEMPLIFICATO con sy diversi \nparametri: \n \t b: " << b << "\t\t\t± " << sigmab << "\n \t a: " << a << "\t\t\t± " << sigmaa << endl << endl;
  //CHI QUADRO
  for(int w=0; w<datax.size(); w++){
  float y;
  y=datax.at(w)*b+a;
  yattesi.push_back(y);
  //cout << yattesi.at(w) << endl;
  }
  for(int i=0; i<yattesi.size() && i<datay.size() && i<incertezzey.size(); i++){
    chi = pow((datay.at(i)-yattesi.at(i))/(incertezzey.at(i)), 2);
    chiq.push_back(chi);
    //cout << chiq.at(i) << endl;
  }
  for(int y=0; y<chiq.size(); y++){
    chitot += chiq.at(y);
  }
  //ERRORE A POSTERIORI
  for(int z=0; z<yattesi.size() && z<datay.size(); z++){
    float num= pow((datay.at(z)-yattesi.at(z)), 2);
    errpost.push_back(num);
  }
  for(int a=0; a<errpost.size(); a++){
    dytot += errpost.at(a);
  }

  dy = sqrt (dytot/(datax.size()-2));
  cout << "chi quadro \t \t \t \t \t" << chitot << endl << endl;
  cout << "errore a posteriori \t \t" << dy <<endl;
  }
  if(x==3){
  //CALCOLI FORMA GENERALE
  double sommax, sommay, sommax2, sommay2, sommaxy, delta, b, a, sigmaa, sigmab;
  int i=0;
  while(i<datax.size()){
    sommax += datax.at(i);
    sommax2 += pow(datax.at(i), 2);
    sommay += datay.at(i);
    sommay2 += pow(datay.at(i), 2);
    sommaxy += datax.at(i)*datay.at(i);
     i++;
  }
  delta = datax.size()*sommax2 - sommax*sommax;
  b = (datax.size()*sommaxy-sommax*sommay)/delta;
  a = (sommax2*sommay-sommax*sommaxy)/delta;
  sigmaa = incertezzey.at(1)*sqrt(sommax2/delta);
  sigmab = incertezzey.at(1)*sqrt(datax.size()/delta);

  double sumx, sumy, sumx2, sumy2, sumxy, sums, sums2, Delta, m, q, sigmam, sigmaq, alfa;
  i=0;
  while(i<incertezzex.size() && i<incertezzey.size()){
    alfa = sqrt(pow(incertezzey.at(i), 2)+pow(b, 2)*pow(incertezzex.at(i), 2));
    incertezzetot.push_back(alfa);
    i++;
  }
  i=0;
  while(i<incertezzetot.size()){
    sigmatot += (pow(incertezzetot.at(i), (-2)));
    i++;
  }
  i=0;
  while(i<datax.size()){
    sumx += (datax.at(i))/pow(incertezzetot.at(i), 2);
    sumx2 += (pow(datax.at(i), 2)/pow(incertezzetot.at(i), 2));
    sumy += (datay.at(i))/pow(incertezzetot.at(i), 2);
    sumy2 += (pow(datay.at(i), 2)/pow(incertezzetot.at(i), 2));
    sumxy += (datay.at(i)*datax.at(i))/pow(incertezzetot.at(i), 2);
   i++;
  }
  Delta = (sigmatot*sumx2)-pow(sumx, 2);
  m = (1/Delta)*(sigmatot*sumxy-sumx*sumy);
  q = (1/Delta)*(sumx2*sumy-sumx*sumxy);
  sigmam= sqrt((1/Delta)*sigmatot);
  sigmaq = sqrt((1/Delta)*sumx2);

  //CHI QUADRO
  for(int w=0; w<datax.size(); w++){
  float y;
  y=datax.at(w)*m+q;
  yattesi.push_back(y);
  }
  for(int i=0; i<yattesi.size() && i<datay.size() && i<incertezzey.size(); i++){
    chi = pow((datay.at(i)-yattesi.at(i))/(incertezzey.at(i)), 2);
    chiq.push_back(chi);
  }
  for(int y=0; y<chiq.size(); y++){
    chitot += chiq.at(y);
  }

  //ERRORE A POSTERIORI
  for(int z=0; z<yattesi.size() && z<datay.size(); z++){
    float num= pow((datay.at(z)-yattesi.at(z)), 2);
    errpost.push_back(num);
  }
  for(int a=0; a<errpost.size(); a++){
    dytot += errpost.at(a);
  }

  dy = sqrt (dytot/(datax.size()-2));

  //OUTPUT
  cout << "retta interpolante \t\t y =  " << m << " x + " << q << endl << endl;
  cout << "parametri: \n \t m: " << m << "\t\t\t± " << sigmam << "\n \t q: " << q << "\t\t\t± " << sigmaq << endl << endl;
  cout << "chi quadro \t \t \t \t \t" << chitot << endl << endl;
  cout << "errore a posteriori \t \t" << dy <<endl;

  ofstream fout("risultati.txt");
  fout << "retta interpolante \t\t y =  " << m << " x + " << q << endl << endl;
  fout << "parametri: \n \t m: " << m << "\t\t\t± " << sigmam << "\n \t q: " << q << "\t\t± " << sigmaq << endl << endl;
  fout << "chi quadro \t \t \t \t \t" << chitot << endl << endl;
  fout << "errore a posteriori \t \t" << dy <<endl;
  fout.close();
  }
  return 0;
}
