#include "df_ergodic.h"
#include "coord.h"
#include <cmath>
#include <omp.h>
//#include"gnuplt.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
//#include "Pi.h"
#include <gsl/gsl_complex.h>
#include "math_sphharm.h"
#include <format>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <chrono>

inline std::istream& skip_char(std::istream& i, const char c) {
	char x;
	i >> x;
	if (x == c) return i;
	return i.putback(x);
}

inline void SwallowRestofLine(std::istream& from)
{
	char c;
	do from.get(c); while (from.good() && c != '\n');
}

void read_yanc_binary(std::ifstream& is, double& mass, std::vector<coord::PosCar>& Xs,
	std::vector<coord::VelCar>& Vs) {
	int N = Xs.size();
	Xs.clear();
	Xs.resize(N);
	//printf("Reading %d particles in binary\n", N);
	char option[100];
	do {
		skip_char(is, '#') >> option;
		SwallowRestofLine(is);
	} while (strcmp(option, "start_data"));
	int w;
	is.read((char*)&w, sizeof(int));// printf("%d\n", w);
	std::vector<float> ms(N), xs(3 * N), vs(3 * N);
	is.read((char*)&ms[0], sizeof(float) * N);
	is.read((char*)&xs[0], sizeof(float) * 3 * N);
	is.read((char*)&vs[0], sizeof(float) * 3 * N);
	is.close();
	mass = ms[N - 1];
	for (int i = 0; i < N; i++) {
		Xs[i].x = xs[3 * i]; Xs[i].y = xs[3 * i + 1]; Xs[i].z = xs[3 * i + 2];
		Vs[i].vx = vs[3 * i]; Vs[i].vy = vs[3 * i + 1]; Vs[i].vz = vs[3 * i + 2];
	}
}

class rsquared_rho : public math::IFunction {
	potential::Isochrone& IsoPot;
public:
	rsquared_rho(potential::Isochrone& _IsoPot) : IsoPot(_IsoPot) {}
	virtual void evalDeriv(const double r, double* val, double* der, double* der2) const {
		if (val) *val = r * r * IsoPot.density(r);
	}
	virtual unsigned int numDerivs() const { return 1; }
};

class radius_getter : public math::IFunction {
	math::IFunction& rsqrho;
	double xi, Mtot;
public:
	radius_getter(math::IFunction& _rsqrho, double _xi) : rsqrho(_rsqrho), xi(_xi) {
		Mtot = math::integrate(rsqrho, 0, 100, 1e-6);
	}
	virtual void evalDeriv(const double r, double* val, double* deriv, double* deriv2) const {
		if (val) *val = math::integrate(rsqrho, 0, r, 1e-6) - xi * Mtot;
		if (deriv) *deriv = rsqrho.value(r);
	}
	virtual unsigned int numDerivs() const { return 1; }
};

void createGrid(std::vector<double>& grid) {


	for (int i = 0; i < grid.size(); i++) {
		grid[i] = pow(10, 2 * (i - 50.0) / double(grid.size()));
	}
}

void createEqualGrid(std::vector<double>& grid) {
	potential::Isochrone IsoPot(1, 1);
	rsquared_rho rsqRho(IsoPot);
	for (int i = 0; i < grid.size(); i++) {
		double xi = (double)i / grid.size();
		radius_getter rget(rsqRho, xi);
		grid[i] = math::findRoot(rget, 0, 200, 1e-6);
	}
}
std::vector<double> get_potcoefs(std::vector<coord::PosCar> Xs, double massi, double r, math::SphHarmIndices ind) {
	std::vector<double> pot(ind.size(), 0.0);

	for (int i = 0; i < Xs.size(); i++) {
		double mass = massi;
		coord::PosCar xv = Xs[i];
		double ri = sqrt(pow_2(xv.x) + pow_2(xv.y) + pow_2(xv.z));
		double costi = xv.z / ri;
		double sinti = sqrt(1 - costi * costi);
		double phi = atan2(xv.y, xv.x);
		double tau = costi / (1 + sinti);
		double K2 = 4 * M_PI * mass;
		for (int m = ind.mmin(); m <= ind.mmax; m++) { 
			int lmin = ind.lmin(m);
			if (lmin > ind.lmax) {
				continue;
			}
			double sinmti = sin(abs(m) * phi);
			double cosmti = cos(abs(m) * phi);
			int size = ind.lmax - abs(m) + 1;
			double* arr = new double[size];
			math::sphHarmArray(ind.lmax, abs(m), tau, arr);
			for (int l = lmin; l <= ind.lmax; l++) {
				double K1 = K2 / (2 * l + 1);
				int n = ind.index(l, m);
				double Ylm = m == 0 ? arr[l - abs(m)] :
					(m < 0 ? arr[l - abs(m)] * sinmti * M_SQRT2 :
					arr[l - abs(m)] * cosmti * M_SQRT2);
				if (ri > r) {
					
					pot[n] += -K1 * Ylm * pow(r / ri, l) / ri ;
				}
				else {
					pot[n] += -K1 * Ylm * pow(ri / r, l + 1) / ri;
				}
			}
			delete[] arr;
		}
	}
	return pot;
}

std::vector<std::vector<double>> get_dipolecoefs(std::vector<coord::PosCar> Xs, double massi, std::vector<double> r, math::SphHarmIndices ind, bool in) {
	std::vector<std::vector<double>> pot(r.size(),std::vector<double>(ind.size(), 0.0));
	for (int k = 0; k < r.size(); k++) {
		//std::cout << " " << k;
		for (int i = 0; i < Xs.size(); i++) {
			double mass = massi;
			coord::PosCar xv = Xs[i];
			double ri = sqrt(pow_2(xv.x) + pow_2(xv.y) + pow_2(xv.z));
			if (ri > r[k] && in) {
				continue;
			}
			if (ri < r[k] && not in) {
				continue;
			}
			//std::cout << "hi";
			double costi = xv.z / ri;
			double sinti = sqrt(1 - costi * costi);
			double phi = atan2(xv.y, xv.x);
			double tau = costi / (1 + sinti);
			//double* arr = static_cast<double*>(_malloca(size * sizeof(double)));
			//double* arr = tmptrig +2 * ind.mmax;
			for (int m = ind.mmin(); m <= ind.mmax; m++) {
				int lmin = ind.lmin(m);
				if (lmin > ind.lmax) {
					continue;
				}
				int size = ind.lmax - m + 1;
				double* arr = new double[size];
				double sinmti = sin(abs(m) * phi);
				double cosmti = cos(abs(m) * phi);
				math::sphHarmArray(ind.lmax, abs(m), tau, arr);
				for (int l = lmin; l <= ind.lmax; l++) {
					int n = ind.index(l, m);
					double Ylm = m == 0 ? arr[l - abs(m)] :
						(m > 0 ? arr[l - abs(m)] * cosmti * M_SQRT2 :
							arr[l - abs(m)] * sinmti * M_SQRT2);
					//std::cout << Ylm;
					pot[k][n] += Ylm;
				}
				delete[] arr;
			}
		}
	}

	return pot;
}
std::vector<std::vector<double>> getBaseCoefs(std::vector<coord::PosCar> Xpos, std::vector<double> r, math::SphHarmIndices ind, double mass) {
	std::vector<std::vector<double>> integral(r.size(), std::vector<double>(ind.size(), 0.0));
	for (int j = 0; j < r.size(); j++) {
		std::vector<double> coefs = get_potcoefs(Xpos, mass,r[j],ind);
		for (int m = 0; m < coefs.size(); m++) {
			integral[j][m] += coefs[m] ;
		}
	}
	return integral;
}


int main()
{
	math::SphHarmIndices ind(2, 2, coord::ST_NONE);
	std::string suffix;

	std::vector<double> grid(500);
	std::vector<std::vector<double>> TotalDipolePower;
	std::vector<std::vector<double>> TotalQuadPower;
	std::vector<std::vector<double>> TotalDipole;
	std::vector<std::vector<int>> TotalCount;
	int N = 1000; //Particles
	int Nr = 1000; // Number of realizations
	int T = 250;
	int dT = 25;
	int Nsteps = T / dT + 1;
	int ExpNum = N * Nr / grid.size();

	for (int k = 0; k < Nsteps; k++) {
		std::vector<double> DipoleCoef(grid.size(), 0.0);
		std::vector<double> dipole_power(grid.size(), 0.0);
		std::vector<double> quad_power(grid.size(), 0.0);
		if (k < 10) {
			suffix = ".000" + std::to_string(k);
		}
		else {
			suffix = ".00" + std::to_string(k);
		}

		std::cout << "########## T = " << k * dT << " ##########\n";

#pragma omp parallel
		{
			std::vector<double> local_dipole_power(grid.size(), 0.0);
			std::vector<double> local_quad_power(grid.size(), 0.0);
			std::vector<double> local_dipole_coef(grid.size(), 0.0);
#pragma omp for
			for (int i = 0; i < Nr; i++) {
				std::string filename = "D:\\config_files\\IsochICs_" + std::to_string(i) + suffix;
				std::ifstream PtrFile(filename, std::ios::binary);
				double mass;
				std::vector<coord::PosCar> Xs(N);
				std::vector<coord::VelCar> Vs(N);
				// Read binary data from file

				read_yanc_binary(PtrFile, mass, Xs, Vs);
				std::vector<std::vector<double>> TempDipoleInCoef = get_dipolecoefs(Xs, mass, grid, ind, true);
				std::vector<std::vector<double>> TempDipoleOutCoef = get_dipolecoefs(Xs, mass, grid, ind, false);
				for (int j = 0; j < grid.size(); j++) {
					double rhoyin = TempDipoleInCoef[j][1];
					double rhoxin = TempDipoleInCoef[j][3];
					double rhozin = TempDipoleInCoef[j][2];
					double rhoyout = TempDipoleOutCoef[j][1];
					double rhoxout = TempDipoleOutCoef[j][3];
					double rhozout = TempDipoleOutCoef[j][2];
					double costheta = (rhoyin * rhoyout + rhoxin * rhoxout + rhozin * rhozout)/ (sqrt(pow_2(rhoyin) + pow_2(rhoxin) + pow_2(rhozin)) * sqrt(pow_2(rhoyout) + pow_2(rhoxout) + pow_2(rhozout)));
					local_dipole_coef[j] += costheta / (double)Nr;
				}
				

				std::vector<std::vector<double>> TempCoef = getBaseCoefs(Xs, grid, ind, mass); //getCoefs(Xs, grid, 1.0, ind, mass);
				for (int s = 0; s < grid.size(); s++) {
					double temp_c1 = pow_2(TempCoef[s][1]) + pow_2(TempCoef[s][2]) + pow_2(TempCoef[s][3]);
					local_dipole_power[s] += temp_c1 / (double)Nr;
					double temp_c2 = pow_2(TempCoef[s][4]) + pow_2(TempCoef[s][5]) + pow_2(TempCoef[s][6]) + pow_2(TempCoef[s][7]) + pow_2(TempCoef[s][8]);
					local_quad_power[s] += temp_c2 / (double)Nr;
				}

			}
#pragma omp critical
			{
				for (int j = 0; j < grid.size();j++) {
					dipole_power[j] += local_dipole_power[j];
					quad_power[j] += local_quad_power[j];
					DipoleCoef[j] += local_dipole_coef[j];
				}
			}
		}
		
		TotalDipolePower.push_back(dipole_power);
		TotalQuadPower.push_back(quad_power);
		TotalDipole.push_back(DipoleCoef);

	}

	for (int i = 0; i < TotalDipolePower.size(); i++) {
		std::string outFileName = "PowerRadialVariation" + std::to_string(i) + ".csv";
		std::ofstream outfile(outFileName);
		for (int j = 0; j < grid.size(); j++) {
			outfile << grid[j] << "," << TotalDipolePower[i][j] << "," << TotalQuadPower[i][j] << '\n';
		}
		outfile.close();
	}
	
	
	for (int i = 0; i < TotalDipole.size(); i++) {
		std::ofstream DipoleFile("DipoleVar" + std::to_string(i) + ".csv");
		for (int j = 0; j < grid.size(); j++) {
			DipoleFile << grid[j] << "," << TotalDipole[i][j] << '\n';
		}
		DipoleFile.close();
	}
	return 0;
}