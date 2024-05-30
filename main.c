#include <stdio.h>
#include <math.h>

double Pi = 3.141592;

struct Weight
{
	double empty;
	double payload;
	double gross;
	double fuselage;
	double rotor;
	double battery;
	double wing;
	double seat;
};

struct Rotor
{
	double N;
	double diameter;
	double radius;
	double sigma;
	double A_blade;
	double A_disk;
	double rpm;
	double omega;
};

struct MainWing
{
	double span;
	double c_root;
	double c_tip;
	double A;
};

struct MissionProfile
{
	double t;
	double h;
	double rho;
	double vi;
	double vc;
	double v_tip;
	double T;
	double P_i;
	double P_0;
	double P_req;
	double E;
};


double Get_rho(double rho, double h)
{
	return rho * exp(-0.0296 * h / 304.8);
}

int main()
{
	double rho = 1.225;
	double g = 0;

	struct Weight W;
	struct Rotor rotor;
	struct MainWing mainwing;

	W.gross = 2100;
	W.empty = 1600;
	W.payload = 500;

	rotor.N = 4;
	rotor.diameter = 3;
	rotor.radius = 0.5 * rotor.diameter;
	rotor.A_blade = 0;
	rotor.sigma = (rotor.A_blade) / (Pi * pow(rotor.radius, 2));

	mainwing.span = 7;
	mainwing.c_tip = 1.2;
	mainwing.c_root = 1.8;
	mainwing.A = (mainwing.c_tip + mainwing.c_root) * mainwing.span;

	double WL = W.gross / mainwing.A;
	double DL = W.gross / (rotor.A_disk * 4);

	struct MissionProfile hover;
	hover.T = W.gross * g;
	hover.t = 0;
	hover.rho = Get_rho();

}

