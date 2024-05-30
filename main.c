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
	double root;
	double tip;
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
	double T_rotor;
	double P_i;
	double P_0;
	double P_req;
	double E;
};


double Get_rho(double rho, double h)
{
	return rho * exp(-0.0296 * h / 304.8);
}

double Get_omega(double rpm)
{
	return rpm * 2 * Pi / 60;
}

double Get_A_wing(double root, double tip, double span)
{
	return (root + tip) * span;
}

double Get_vi(double T_rotor, double rho, double A_disk)
{
	return sqrt(T_rotor / (2 * rho * A_disk));
}

double Get_P_i(double T_rotor, double vi)
{
	return T_rotor * vi;
}

double Get_P_0(double sigma, double Cd0, double rho, double A_disk, double v_tip)
{
	return ((sigma * Cd0) / 8) * rho * A_disk * pow(v_tip, 3);
}

int main()
{
	double rho = 1.225;
	double g = 0;

	double Cd0 = 0.03008;

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
	rotor.A_disk = Pi * pow(rotor.radius, 2);
	rotor.sigma = rotor.A_blade / rotor.A_disk;
	rotor.rpm = 2000;
	rotor.omega = Get_omega(rotor.rpm);
	mainwing.span = 7;
	mainwing.tip = 1.2;
	mainwing.root = 1.8;
	mainwing.A = Get_A_wing(mainwing.root, mainwing.tip, mainwing.span);

	double WL = W.gross / mainwing.A;
	double DL = W.gross / (rotor.A_disk * 4);

	struct MissionProfile hover;
	hover.T = W.gross * g;
	hover.T_rotor = hover.T / rotor.N;
	hover.t = 30 * 60;
	hover.h = 0;
	hover.rho = Get_rho(rho, hover.h);

	hover.vi = Get_vi(hover.T_rotor, hover.rho, rotor.A_disk);

	double v_tip = rotor.omega * (rotor.radius);

	hover.P_i = Get_P_i(hover.T_rotor, hover.vi);
	hover.P_0 = Get_P_0(rotor.sigma, Cd0, hover.rho, rotor.A_disk, v_tip);

	hover.P_req = hover.P_i + hover.P_0;
	hover.E = hover.P_req * (hover.t / 3600);

	struct MissionProfile climb;
	climb.T = W.gross * g;
	
}

