#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <string.h>
#include <time.h>
#include <vector>
#include <list>
#include <algorithm>
#include <functional>
#include "omp.h"

#define _ELECTRON 1
#define _PROTON 2

using namespace std;

//Constants expressed in SI units
const double e0 = 8.85e-12;					//dielectric constant in vacuum
const double electronmass = 9.109e-31;		//electron mass
const double protonmass = 1.67e-27;			//proton mass
const double h2mass = 3.346e-27;			//H2+ ion mass
const double kb = 1.38e-23;					//boltzmann constant
const double q = 1.6e-19;					//elementary charge
const double pi = 3.14159265359;			//pi
const double uma = 1.66053886e-27;			//1/12 of C mass

//Function for obtaining sign of a number
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

//Function for converting a number into a string
template <typename T>
std::string NumberToString(T Number)
{
	std::ostringstream ss;
	ss << Number;
	return ss.str();
}

//Function for generating a number from a uniform distribution in [0,1]
inline double uniform() {
	return rand() / (1. + RAND_MAX);
}

//Function for generating a number from a uniform distribution in [0,L]
inline double uniform(double L) {
	return L * rand() / (1. + RAND_MAX);
}

//Function for generating a velocity from a normal distribution at a given temperature
double velfromnormal(double kTeV, double m) {
	double x, y;
	do {
		x = rand() / (1.0 + RAND_MAX);
	} while (x == 0);
	y = rand() / (1.0 + RAND_MAX);
	return sqrt(kTeV / m * q) * sqrt(-2.0 * log(x)) * cos(2.0 * pi * y);
}

//Function for calculating the mass of a particle
double calcMass(int particletype) {
	double mass = 0.0;
	if (particletype == _ELECTRON) {
		mass = electronmass;
	}
	else if (particletype == _PROTON) {
		mass = protonmass;
	}
	return mass;
}

//Function for calculating the charge of a particle
double calcCharge(int particletype) {
	double charge = 0.0;
	if (particletype == _ELECTRON) {
		charge = -q;
	}
	else if (particletype == _PROTON) {
		charge = +q;
	}
	return charge;
}

//Class Particle for storing main particle parameters (x,y,posx,posy)
class Particle {
public:
	Particle() {
		_wei = 1e6;
	}
	Particle(double wei) {
		_wei = wei;
	}
	Particle(double wei, int particletype) {
		setParticle(particletype);
		_wei = wei;
	}
	Particle(double wei, int particletype, double T) {
		_wei = wei;
		setParticle(particletype);
		startMBparticle(T);
	}
	Particle(double wei, int particletype, double T, double L) {
		_wei = wei;
		setParticle(particletype);
		startMBparticle(T);
		setPosition(uniform(L));
	}
	void setPosition(double x) {
		posx = x;
		posz = 0.0;
	}
	int getType() {
		return _particletype;
	}
	void setParticle(int particletype) {
		mass = calcMass(particletype);
		charge = calcCharge(particletype);
		_particletype = particletype;
	}

	void loadParticle(int particletype, double x, double z,
					double velx, double velz, double wei) {
		mass = calcMass(particletype);
		charge = calcCharge(particletype);
		posx = x;
		posz = z;
		vx = velx;
		vz = velz;
		weight = wei;
		_particletype = particletype;
	}
	double posx, posz, vx, vz, acceleration, mass, charge, weight;
	void startMBparticle(double T);
private:
	double _wei;
	int _particletype;
};

//Function in Particle class for generating a particle from a Maxwell Boltzmann distribution at a given temperature
void Particle::startMBparticle(double T) {
	vx = velfromnormal(T, mass);
	vz = velfromnormal(T, mass);
	weight = _wei;
}

//Function for generating a thermal particle in a region of size L
void createThermalParticle(vector<Particle>& part, vector<double>& field, vector<size_t>& index, const double T,
	const double dx, const double dt, size_t& indexmax, size_t& numpart, double mpw, double L, int particletype)
{
	Particle newpart(mpw, particletype, T, L);																																											//all'array di particelle totali; inizializza la velocità per leap-frog

	int j;																																														//(per calcolo accelerazione vedere sotto)
	double localfield;
	j = (int)(floor(newpart.posx / dx) + 0.1);
	localfield = field[j] * ((j + 1) * dx - newpart.posx) / dx + field[j + 1] * (newpart.posx - j * dx) / dx;
	newpart.acceleration = newpart.charge * localfield / newpart.mass;
	newpart.vx = newpart.vx + newpart.acceleration * 0.5 * dt;
	if (indexmax != 0) {
		part[index[indexmax - 1]] = newpart;
		indexmax = indexmax - 1;
	}
	else {
		part[numpart] = newpart;
		numpart = numpart + 1;
	}
}

//Function for generating a thermal particle between xmin and xmax
void createThermalParticlePerTimestep(vector<Particle>& part, vector<double>& field, vector<size_t>& index, const double T,
	const double dx, const double dt, size_t& indexmax, size_t& numpart, double mpw, double xmin, double xmax, int particletype)
{
	Particle newpart(mpw, particletype, T);
	newpart.posx = xmin + uniform(xmax-xmin);
	newpart.posz = 0.0;

	int j;																																														//(per calcolo accelerazione vedere sotto)
	double localfield;
	j = (int)(floor(newpart.posx / dx) + 0.1);
	localfield = field[j] * ((j + 1) * dx - newpart.posx) / dx + field[j + 1] * (newpart.posx - j * dx) / dx;
	newpart.acceleration = newpart.charge * localfield / newpart.mass;
	newpart.vx = newpart.vx + newpart.acceleration * 0.5 * dt;
	if (indexmax != 0) {
		part[index[indexmax - 1]] = newpart;
		indexmax = indexmax - 1;
	}
	else {
		part[numpart] = newpart;
		numpart = numpart + 1;
	}
}

//Function for calculating the particle numerical density
void numericdensityParticle(vector<double>& partnumdensity, const int particletype, vector<Particle>& part,
	const size_t numpart, const int numcells, const double dx, const double L) {
	int i;
#pragma omp parallel for shared (partnumdensity) private (i) schedule (static)																															//la particella viene distribuita tra i due nodi in mezzo ai quali si trova
	for (i = 0; i < numcells; i++) {
		partnumdensity[i] = 0;
	}
	int j;
	for (i = 0; i < numpart; i++) {
		//if (part[i].posx < L) {
		if (part[i].weight > 0) {
			j = (int)(floor(part[i].posx / dx) + 0.1);
			if (part[i].getType() == particletype) {
				partnumdensity[j] = partnumdensity[j] + (1.0 / dx) * part[i].weight * (dx * (j + 1) - part[i].posx) / dx;
				partnumdensity[j + 1] = partnumdensity[j + 1] + (1.0 / dx) * part[i].weight * (part[i].posx - dx * j) / dx;
			}
		}
	}
	partnumdensity[0] *= 2.0;
	partnumdensity[numcells - 1] *= 2.0;
}

//Function for calculating the average particle energy
void calculateaverageeneryParticle(vector<double>& aveE, vector<double>& aveTotalE, const int particletype, vector<Particle>& part,
	const size_t numpart, const int numcells, const double dx, const double L) {
	int i;
#pragma omp parallel for shared (aveE,aveTotalE) private (i) schedule (static)
	for (i = 0; i < numcells; i++) {
		aveE[i] = 0;
		aveTotalE[i] = 0;
	}
	int j;
	vector<int> nparts(aveE.size(), 0);

	for (i = 0; i < numpart; i++) {
		//if (part[i].posx < L) {
		if (part[i].weight > 0) {
			j = (int)(floor(part[i].posx / dx) + 0.1);
			if (part[i].getType() == particletype && part[i].weight != 0) { //TODO serve il secondo if?
				aveE[j] = aveE[j]
					+ 0.5 * part[i].mass * part[i].vx * part[i].vx / q;
				aveE[j + 1] = aveE[j + 1]
					+ 0.5 * part[i].mass * part[i].vx * part[i].vx / q;
				aveTotalE[j] = aveTotalE[j]
					+ 0.5 * part[i].mass * (part[i].vx * part[i].vx + part[i].vz * part[i].vz) / q;
				aveTotalE[j + 1] = aveTotalE[j + 1]
					+ 0.5 * part[i].mass * (part[i].vx * part[i].vx + part[i].vz * part[i].vz) / q;
				nparts[j]++;
				nparts[j + 1]++;
			}
		}
	}
	//la particella viene distribuita tra i due nodi in mezzo ai quali si trova
//#pragma omp parallel for shared (electronaveE, electronaveE, nelecs, nions) schedule (static)
	for (int j = 0; j < numcells; j++) {
		if (nparts[j] > 0.0) aveE[j] = aveE[j] / nparts[j];
		if (nparts[j] > 0.0) aveTotalE[j] = aveTotalE[j] / nparts[j];
	}

}

//Function for setting the particle charge density to 0
void chargeDensityReset(vector<double>& rho,
	const int numcells, const double dx, const double L)
{
	int i;
#pragma omp parallel for shared (rho) private (i) schedule (static)																													
	for (i = 0; i < numcells; i++) {
		rho[i] = 0.0;
	}
}

//Function for calculating the charge density
void depositChargeDensity(vector<double>& rho, vector<double>& particlenumdensity, int particletype,
	const int numcells, const double dx, const double L)
{
	double mycharge = calcCharge(particletype);
	for (int i = 0; i < numcells; i++) {
		rho[i] = rho[i] + particlenumdensity[i] * mycharge;
	}
}

//Direct Poisson solver (Cavenago)
void poisson2mc(vector<double>& rho, double* vnew, double* ez, double* eztot, double* po, const double dx, const int numcells, const int potSX, const int potDX) {
	int ii, nz0;
	double dzz, dza;
	double vdiff0, vstep0, ez0, ez8;
	vstep0 = 0;

	// algorithm start
	dzz = dx / e0;
	dza = dx;
	nz0 = numcells - 1; //index of the last valid element in the vector
	po[0] = 0;
	ez[0] = rho[0] * dzz / 2;

	// integrate for u (that is po and ez)
	for (ii = 1; ii <= nz0; ii++) {
		ez[ii] = ez[ii - 1] + dzz * rho[ii];
	}

	for (ii = 1; ii <= nz0; ii++) {
		po[ii] = po[ii - 1] - dza * ez[ii - 1];
	}

	//adjust for boundary for combinations of Neumann and Dirichlet
	if ((potSX == 999) && (potDX == 999))
		vstep0 = 0; //wrong; two Neumann boundaries are impossible
	else if (potDX == 999) { //Neumann at right end (so Direchlet at left end)
		ez8 = (ez[nz0 - 1] + ez[nz0]) / 2;
		vstep0 = -dza * ez8;
		for (ii = 0; ii <= nz0; ii++) {
			vnew[ii] = potSX + ii * vstep0 + po[ii]; //add v (Laplace part) to u (the po term)
			eztot[ii] = ez8 + ez[ii];
		}
	}
	else if (potSX == 999) { //Neumann at left end (so Direchlet at right end)
		vdiff0 = potDX - po[nz0];
		for (ii = 0; ii <= nz0; ii++) {
			vnew[ii] = vdiff0 + po[ii];
			eztot[ii] = ez[ii];
		}
	}
	else {// no Neuman, (so Direchlet at left and right end)
		vdiff0 = (double)(potDX - potSX) - po[nz0];  //part of potential increase assigned to v
		vstep0 = vdiff0 / nz0;
		ez0 = -vstep0 / dza;						 //part of electric field assigned to v
		for (ii = 0; ii <= nz0; ii++) {
			vnew[ii] = potSX + ii * vstep0 + po[ii];
			eztot[ii] = ez0 + ez[ii];
		}
	}
}

//Jacobi Poisson solver
void pot(double* vold, double* vnew, vector<double>& rho, const double dx, const int numcells) {
	vnew[0] = vold[1] + 0.5 * std::pow(dx, 2) * rho[0] / e0;
	vnew[numcells - 1] = 0;

	int i;
#pragma omp parallel for shared(vold,vnew,rho) private (i) schedule(static)
	for (i = 1; i < numcells - 1; i++) {
		vnew[i] = 0.5 * (vold[i + 1] + vold[i - 1]) + 0.5 * std::pow(dx, 2) * rho[i] / e0;
	}
}

//Function for calculating the electric field
void efieldcalculator(double* vnew, vector<double>& field, const double dx, const int numcells) {
	field[0] = 0.;	
	field[numcells - 1] = -1.0 * (vnew[numcells - 1] - vnew[numcells - 2]) / dx;

	int i;
#pragma omp parallel for shared (vnew, field) private (i) schedule (static)
	for (i = 1; i < numcells - 1; i++) {
		field[i] = -0.5 * (vnew[i + 1] - vnew[i - 1]) / dx;
	}
}

//Function for calculating particle acceleration
void giveacceleration(vector<double>& field, vector<Particle>& part, const size_t numpart, const double dx, const double L) {
	int i,j;
	double localfield;
#pragma omp parallel for shared (field, part) private (i,j,localfield) schedule (static)
	for (i = 0; i < numpart; i++) {
		if (part[i].posx > 0 && part[i].posx < L) {
			j = (int)(floor(part[i].posx / dx) + 0.1);
			localfield = field[j] * ((j + 1) * dx - part[i].posx) / dx + field[j + 1] * (part[i].posx - j * dx) / dx;
			part[i].acceleration = part[i].charge * localfield / part[i].mass;
		}
	}
}

//Function for moving particles
void move(vector<Particle>& part, vector<size_t>& index, const size_t numpart,
	const double dt, const double L, size_t& indexmax,
	bool sey_me, size_t& sey_counter, size_t& out_counter, size_t& z_out_counter, bool mirrorleft = false, double Zmax = 1e20)
{
	sey_counter = 0;
	out_counter = 0;
	
	int i;
#pragma omp parallel for shared (part) private (i) schedule (static)
	for (i = 0; i < numpart; i++) {
		if (part[i].posx < L) {
			part[i].vx = part[i].vx + part[i].acceleration * dt;
			part[i].posx = part[i].posx + part[i].vx * dt;
			//Mirror particles on the left false by default
			if (mirrorleft)
				if (part[i].posx <= 0) {
					part[i].posx = -part[i].posx;
					part[i].vx = -part[i].vx;
				}
			part[i].posz = part[i].posz + part[i].vz * dt;
		}
	}
	/** SEY */
	if (sey_me) {
		double sey_electron = 1.0;
		double sey_ion = 0.0;
		for (int i = 0; i < numpart; i++) {
			if (part[i].posx >= L && part[i].weight > 0) {
				if (part[i].getType() == _ELECTRON && uniform() < sey_electron) {
					part[i].vx *= -uniform();
					part[i].posx = L - uniform(fabs(part[i].vx * dt));
					sey_counter++;
				}
				if ((part[i].getType() == _PROTON) && uniform() < sey_ion) {
					part[i].setParticle(_ELECTRON);
					part[i].vx = -fabs(velfromnormal(1.0, part[i].mass));
					part[i].posx = L - uniform(fabs(part[i].vx * dt));
					sey_counter++;
				}
			}
		}
	}

	/** OUT OF DOMAIN */
#pragma omp parallel for shared(part,index) private (i)
	for (i = 0; i < numpart; i++) {
		if (part[i].posx >= L && part[i].weight > 0) {
			part[i].weight = 0;
			index[indexmax] = i;
#pragma omp critical
			{
				indexmax++;
				out_counter++;
			}
		}
	}

	if (!mirrorleft) {
#pragma omp parallel for shared(part,index) private (i)
		for (i = 0; i < numpart; i++) {
			if (part[i].posx <= 0.0 && part[i].weight > 0) {
				part[i].weight = 0;
				index[indexmax] = i;
#pragma omp critical
				{
					indexmax++;
					out_counter++;
				}
			}
		}
	}

	/** 1D 2V check */
	for (int i = 0; i < numpart; i++){
		if (part[i].posz > Zmax && part[i].weight > 0) {
			part[i].weight = 0;
			index[indexmax] = i;
			indexmax = indexmax + 1;
			z_out_counter++;
		}
	}
}

//Thermostat
size_t thermostat(std::vector<Particle>& part, std::vector<size_t>& index, const size_t numpart,
	const double dt, const double tau, bool keepdirection,
	const double T, const int particletype,
	const double Lmin, const double Lmax, size_t& indexmax)
{
	/// Fraction of particles to thermostat
	const double fractionOfParticlesToThermostat = dt / tau;
	size_t counter = 0;
	int i;
#pragma omp parallel for shared (part) private (i) schedule (static)
	for (i = 0; i < numpart; i++) {
		if (part[i].getType() == particletype &&
			part[i].posx < Lmax && part[i].posx > Lmin && part[i].weight > 0) {
			double R = rand() / (1.0 + RAND_MAX);
			if (R < fractionOfParticlesToThermostat) {
				double v2 = velfromnormal(T, part[i].mass);
				if (keepdirection && part[i].vx != 0.0) v2 *= part[i].vx / fabs(part[i].vx);
				part[i].vx = v2;
#pragma omp critical
				counter++;
			}
		}
	}
	return counter;
}

int main(int argc, char * argv[]) {
#pragma omp parallel
	{
		int a = 0;
		if (omp_get_thread_num() == 0)
			printf("Number of threads %d\n", omp_get_num_threads());
	}

	std::cout << "****  1D PARTICLE-IN-CELL  ****\n";
	std::cout << " build " << __DATE__ << " " << __TIME__ << "\n";

	srand(time(0));
	clock_t tStart = clock();
	
	//Warm start
	int loadIt = 0;
	float loadTime = 0.0;
	if (argc > 1) {
		loadIt = atoi(argv[1]);
		if (argc > 2) {
			loadTime = atof(argv[2]);
		}
		std::cout << " loading step " << loadIt << "\n";
		std::cout << " initial time " << loadTime << "\n";
	}

	//Set simulation parameters
	const double dt = 5e-12;              //simulation timestep (s)
	const double dx = 2.5e-6;             //cell size (m)
	const double tmax = 4.2e-8;           //maximum time (s)
	const double L = 0.0005;              //domain length (m)
	double Lz_max = 0.0;                  //domain length in z-direction (m)

	double desired_n = 4.0e17;            //target plasma density (m^-3)
	const double mpw = 5e8;              //macroparticle weight
	const size_t nummaxpart = 4e6;        //max number of particles

	const int timesteps = tmax / dt + 1;  //number of timesteps
	const int numcells = L / dx + 1;      //number of cells

	long double res, norm, norm0;         //double variables for Jacobi solver
	int i;                                //iteration and timestep indexes

	//Flags for additional routines
	bool thermostat_me = false;           //if true, electrons are heated by the thermostat
	bool newpartpertimestep = true;       //if true, new MP are generated at each timestep

	int saveeverytotiterations = 10;      //defines how frequently output data are stored
	int savepartlist = 500;               //save particle list every tot iterations

	//Plasma parameters
	double Te = 10.0;				      //Electron temperature (eV)
	const double Ti = 1.0;                //Positive ion temperature (eV)

	//Define how many macroparticles will be generated at...
	size_t thermalsourceterm0 = desired_n * L / mpw;   //...first iteration
	size_t ithermalsourceterm1 = 0.25 * desired_n *    //... subsequent iterations
		                         pow(8 * q * Ti / pi / protonmass, 0.5) *
		                         dt / mpw;
	size_t ethermalsourceterm1 = 0.25 * desired_n *    //... subsequent iterations
								 pow(8 * q * Te / pi / electronmass, 0.5) *
							     dt / mpw;
	double BohmFlux = 0.5 * desired_n * pow(q * Te / protonmass, 0.5);
	size_t volumetricRate = 2 * BohmFlux * dt / mpw;   //Number of pairs re-injected per timestep
	double gen_xmin = 0;
	double gen_xmax = L;

	//Thermostat parameters
	double Ltherm = 0.0002;               //Range over which thermostat can act
	double Lmin = 0 + dx;                 //Starting position for thermostat action
	double Lmax = Lmin + (L - 2 * dx);    //Ending position for thermostat action
	double tau_i = 0.015e-6;              //characteristic time for ion heating (s)
	double tau_e = 0.00015e-6;            //characteristic time for electron heating (s)

	//Variables initialisation
	size_t numpart = 0;
	size_t indexmax = 0;
	size_t actualnumpart;
	double* vold = new double[numcells];
	double* vnew = new double[numcells];
	double* temp = new double[numcells];
	double* ez = new double[numcells];
	double* eztot = new double[numcells];
	double* po = new double[numcells];

	vector<double> field(numcells);
	vector<double> rho(numcells);
	vector<double> electronnumdensity(numcells);
	vector<double> ionnumdensity(numcells);
	vector<double> nionnumdensity(numcells);
	vector<double> electronaverageenergy(numcells);
	vector<double> ionaverageenergy(numcells);
	vector<double> nionaverageenergy(numcells);

	vector<double> electronaveragetotalenergy(numcells);
	vector<double> ionaveragetotalenergy(numcells);
	vector<double> nionaveragetotalenergy(numcells);
	vector<double> molionaveragetotalenergy(numcells);

	vector<double> molionnumdensity(numcells);
	vector<double> molionaverageenergy(numcells);

	vector<Particle> part(nummaxpart, Particle(mpw));
	vector<size_t> index(nummaxpart);

	//Decide simulation type

#define DEBYESHEATH

#ifdef DEBYESHEATH
	thermostat_me = true;
	newpartpertimestep = true;
#endif

	//Specify if Jacobi solver is used
//#define JACOBISOLVER

	//Print some parameters on screen
	std::cout << "**** SIMULATION PARAMETERS ****\n";
	std::cout << " --> timestep dt = " << dt << " s\n";
	std::cout << " --> cell size dx = " << dx << " m\n";
	std::cout << " --> domain size L = " << L << " m\n";
	std::cout << " --> target density n0 = " << desired_n << " m^-3\n";
	if (Lz_max != 0.0f)
		std::cout << " --> remove particles if z > " << Lz_max << " m \n";
	std::cout << " --> initial MP number = " << thermalsourceterm0 << "\n";
	if (newpartpertimestep) {
		std::cout << " --> " << volumetricRate << " (e,H+) pairs are injected at each iteration\n";
	}
	if (thermostat_me) {
		std::cout << " --> thermostat acts between " << Lmin << " and " << Lmax << " m\n";
		std::cout << " --> e are resampled with tau " << tau_e << "(" << dt / tau_e * (desired_n * L / mpw) << ")\n";
		std::cout << " --> ions are resampled with tau " << tau_i << "(" << dt / tau_i * (desired_n * L / mpw) << ")\n";
	}

	double t_pe = 1.0 / (pow(desired_n * q * q / electronmass / e0, 0.5) / 2.0 / 3.14159);          //inverse of plasma frequency
	double l_D = pow((e0 * kb / q / q) / (desired_n / Te / 11604 + desired_n / Ti / 11604), 0.5);   //generalised Debye length
	std::cout << " --> dt is " << dt / t_pe << " times Tau_pe (should be < 0.2)\n";
	std::cout << " --> dx is " << dx / l_D << " times the Debye length (should be < 3)\n";

	//Open convergence history file
	ofstream fout;
	vector<double> x(numcells, 0.0);
	string filename;
	string convhisto = "convergence_history";
	fout.open(convhisto);
	fout << "%t \t ne(x) \t ni(x) \t nn(x) \t rho(x) \t phi(x) \t <Ee> \t <Ei> \t <En>\n";
	fout.close();
	std::cout << " --> Convergence history data saved in ''" << convhisto << "''\n";

    //Set variables used by Poisson solver to zero
#pragma omp parallel for shared(vold,vnew,rho) private (i) schedule (static)
	for (i = 0; i < numcells; i++) {
		vold[i] = 0;
		vnew[i] = 0;
		rho[i] = 0;
		temp[i] = 0;
		ez[i] = 0;
		eztot[i] = 0;
		po[i] = 0;
	}

	//Set particle weight to zero
#pragma omp parallel for shared(part) private (i) schedule (static)
	for (i = 0; i < part.size(); i++) {
		part[i].weight = 0;
	}

	//Main iteration cycle
	//1. create MParticles
	//2. calculate charge density
	//3. solve Poisson equation
	//4. calculate electric field
	//5. move particles
	//+ additional routines (i.e. thermostat, secondary emission...)
	std::cout << "**** START MAIN ITERATION CYCLE ****\n";

	for (i = 0; i < timesteps; i++) {

		//Warm start
		if (loadIt != 0 && i == 0) {
			std::cout << "****   Loading particles   ****\n";
			filename = "partlist" + NumberToString(loadIt) + ".txt";
			ifstream loadfile(filename);
			string temp;
			int m = 0;
			const char sep[2] = { '\t' };
			vector<string> line(5);

			while (getline(loadfile, temp)) {
				int k = 0;
				for (size_t p = 0, q = 0; p != temp.npos; p = q) {
					line[k] = temp.substr(p + (p != 0), (q = temp.find(sep, p + 1)) - p - (p != 0));
					k++;
				}
				//type posx posz vx vz
				part[m].loadParticle( atoi(line[0].c_str()), atof(line[1].c_str()), atof(line[2].c_str()),
									  atof(line[3].c_str()), atof(line[4].c_str()), mpw);
				m++;
			}
			loadfile.close();
			numpart = m;

			cout << " --> " << numpart << " particles are loaded\n";
		}//end warm start

		//Generate particles at first iteration
		if (loadIt == 0 && i == 0) {
			for (int j = 0; j < thermalsourceterm0; j++)
				createThermalParticle(part, field, index, Te, dx, dt, indexmax, numpart, mpw, L, _ELECTRON);
			for (int j = 0; j < thermalsourceterm0; j++)
				createThermalParticle(part, field, index, Ti, dx, dt, indexmax, numpart, mpw, L, _PROTON);
		}
		else {
			if (newpartpertimestep) {
				for (int j = 0; j < volumetricRate; j++)
					createThermalParticlePerTimestep(part, field, index, Te, dx, dt, indexmax, numpart, mpw, gen_xmin, gen_xmax, _ELECTRON);
				for (int j = 0; j < volumetricRate; j++)
					createThermalParticlePerTimestep(part, field, index, Ti, dx, dt, indexmax, numpart, mpw, gen_xmin, gen_xmax, _PROTON);
			}
		}

		numericdensityParticle(electronnumdensity, _ELECTRON, part, numpart, numcells, dx, L);
		numericdensityParticle(ionnumdensity, _PROTON, part, numpart, numcells, dx, L);

		calculateaverageeneryParticle(electronaverageenergy, electronaveragetotalenergy, _ELECTRON, part, numpart, numcells, dx, L);
		calculateaverageeneryParticle(ionaverageenergy, ionaveragetotalenergy, _PROTON, part, numpart, numcells, dx, L);

		chargeDensityReset(rho, numcells, dx, L);
		depositChargeDensity(rho, electronnumdensity, _ELECTRON, numcells, dx, L);
		depositChargeDensity(rho, ionnumdensity, _PROTON, numcells, dx, L);

#ifdef JACOBISOLVER
		p = 0;
		int poissonmaxiterations = 300000;
		do{
			p++;
			norm = 0.;
			temp = vold;
			vold = vnew;
			vnew = temp;
			pot(vold, vnew, rho, dx, numcells);
			for (int k = 0; k < numcells; k++){
				res = vnew[k] - vold[k];
				norm = norm + res*res;
			}
			norm=pow(norm,0.5);
			if (p == 1){
				norm0 = norm;
			}
		}while (norm / norm0 > 1.0e-29 && p < poissonmaxiterations);
#else
		poisson2mc(rho, vnew, ez, eztot, po, dx, numcells, 0, 0); //sx dx
#endif

		efieldcalculator(vnew, field, dx, numcells);
		giveacceleration(field, part, numpart, dx, L);
		if (i == 0 && loadIt == 0)
			for (int j = 0; j < numpart; j++)
				part[j].vx = part[j].vx + part[j].acceleration * 0.5 * dt;

		size_t sey_counter = 0;
		size_t out_counter = 0;
		size_t z_out_counter = 0;
		
		if (Lz_max != 0.0f)
			move(part, index, numpart, dt, L, indexmax, false, sey_counter, out_counter, z_out_counter, false, Lz_max);
		else
			move(part, index, numpart, dt, L, indexmax, false, sey_counter, out_counter, z_out_counter, false);

		if (thermostat_me) {
			thermostat(part, index, numpart, dt, tau_e, true, Te, _ELECTRON, Lmin, Lmax, indexmax);
			thermostat(part, index, numpart, dt, tau_i, true, Ti, _PROTON, Lmin, Lmax, indexmax);
		}

		actualnumpart = 0;
		for (int m = 0; m < numpart; m++)
			if (part[m].weight > 0)
				actualnumpart++;

		//Save output data
		if (i % saveeverytotiterations == 0) {
			filename = "iteration" + NumberToString(i + loadIt) + ".txt";
			fout.open(filename.c_str());
			fout << "%x \t ne(x) \t ni(x) \t nn(x) \t rho(x) \t phi(x) \t <Ee> \t <Ei> \t <En> \t <Ee_tot> \t <Ei_tot> \t <En_tot>\n";
			for (int j = 0; j < numcells; j++) {
				fout << j * dx << "\t" << electronnumdensity[j] << "\t" << ionnumdensity[j] << "\t" << nionnumdensity[j]
					<< "\t" << rho[j] << "\t" << vnew[j]
					<< "\t" << electronaverageenergy[j] << "\t" << ionaverageenergy[j] << "\t" << nionaverageenergy[j]
					<< "\t" << electronaveragetotalenergy[j] << "\t" << ionaveragetotalenergy[j] << "\t" << nionaveragetotalenergy[j] << "\n";
			}
			fout.close();
		}
		//Update convergence history file
		{
			int xsave = (int)(0.5 * L / dx);
			fout.open("convergence_history", std::ios_base::app);
			fout << dt * (i + loadIt) << "\t" << electronnumdensity[xsave] << "\t" << ionnumdensity[xsave] << "\t" << nionnumdensity[numcells - 1]
				<< "\t" << rho[xsave] << "\t" << vnew[xsave]
				<< "\t" << electronaverageenergy[xsave] << "\t" << ionaverageenergy[xsave] << "\t" << nionaverageenergy[numcells-1]
				<< "\t" << electronaveragetotalenergy[xsave] << "\t" << ionaveragetotalenergy[xsave] << "\t" << nionaveragetotalenergy[numcells - 1] << "\n";
			fout.close();
		}
		//Save particle list
		if (i % savepartlist == 0 && i > 0) {
			filename = "partlist" + NumberToString(i+loadIt) + ".txt";
			ofstream savelist(filename, std::ios_base::app);
			for (int j = 0; j < numpart; j++) {
				if (part[j].weight > 0) {
					savelist << part[j].getType()
						<< "\t" << part[j].posx
						<< "\t" << part[j].posz
						<< "\t" << part[j].vx
						<< "\t" << part[j].vz << "\n";
				}
			}
			savelist.close();
		}
		//Print some parameters on screen
		std::cout << "   > it = " << i << "/" << timesteps << ", t = "
			<< (dt * (i + loadIt))/1.0E-6 << "us, NMP = " << actualnumpart
			<< " / " << numpart << ", MP out = " << out_counter
			<< ", MP Zout = " << z_out_counter << " \r";
	}//end main cycle

	std::cout << "**** END MAIN ITERATION CYCLE ****\n";
	//At the end of the time cycle, save total execution time
	fout.open("executiontime.txt");
	fout << "Total execution time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " s\n";
	fout.close();
	return 0;
}
