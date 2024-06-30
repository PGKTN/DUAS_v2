#include <stdio.h>
#include <math.h>

double Pi = 3.141592;



//struct Rotor
//{
//	double N;
//	double diameter;
//	double radius;
//	double sigma;
//	double A_blade;
//	double A_disk;
//	double rpm;
//	double omega;
//};
//
//struct MainWing
//{
//	double span;
//	double root;
//	double tip;
//	double A;
//};
//
//struct MissionProfile
//{
//	double t;
//	double h;
//	double rho;
//	double v;
//	double vi;
//	double vc;
//	double v_tip;
//	double L;
//	double T;
//	double Tc;
//	double T_rotor;
//	double P_i;
//	double P_0;
//	double P_req;
//	double E;
//	double TW;
//};


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
	return (root + tip) * span * 0.5;
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

double Cfc(double Re_L, double M)
{

}

int main()
{
	double rho = 1.225;
	double g = 9.81;

	double payload_req = 500;

	double GW = 2100;


	////////////////
	/// Fuselage ///
	////////////////
	double range_max = 0.0;
	double W_fuse = 0.026 * pow(GW, 0.943) * pow(range_max, 0.654);
	printf("W_fuse = %lf [kg]\n", W_fuse);

	////////////
	/// Wing ///
	////////////
	double span_wing = 14;
	double tip_wing = 1.2;
	double root_wing = 1.8;
	double A_wing = (root_wing + tip_wing) * span_wing / 2;
	double AR = pow(span_wing, 2) / A_wing;



	//printf("%lf\n", pow(S_w, 0.758));
	//printf("%lf\n", pow(W_fw, 0.0035));
	//printf("%lf\n", pow(AR / pow(cos(Rambda), 2), 0.6));
	//printf("%lf\n", pow(q, 0.006));
	//printf("%lf\n", pow((100 * tc) / (cos(Rambda)), -0.3));
	//printf("%lf\n", pow((N * w_dg), 0.49));


	///////////
	/// Rod ///
	///////////
	double W_rod = 100;


	////////////////////
	/// Landing_gear ///
	////////////////////
	double W_landing_gear = 100;


	/////////////
	/// Motor ///
	/////////////
	double W_motor = 40 * 4;



	////////////////////////
	/// Mission Analysis ///
	////////////////////////
	double N_rotor = 4;

	double T = GW * g;
	double T_rotor = T / N_rotor;

	double d_rotor = 3;
	double r_rotor = d_rotor / 2;
	double A_blade = 0;
	double A_disk = Pi * pow(r_rotor, 2);
	double sigma = 0.1696;

	double rpm_climb = 1500;
	double omega_climb = Get_omega(rpm_climb);
	double v_tip_climb = omega_climb * r_rotor;

	double DL = GW * g / (A_disk * N_rotor);  // Loading per disk. disk는 총 4개.

	double Cd0_climb = 0.03407;
	double rho_climb = Get_rho(rho, 0);
	double t_climb = 3 * 60;
	double vi = Get_vi(T_rotor, rho_climb, A_disk);
	double v_climb = 3.333333;
	double P_i_climb = Get_P_i(T_rotor, vi + v_climb);
	double P_0_climb = Get_P_0(sigma, Cd0_climb, rho_climb, A_disk, v_tip_climb);
	double P_climb = P_i_climb + P_0_climb;
	double E_climb = P_climb * (t_climb / 3600);
	printf("E_climb = %lf [Wh]\n", E_climb);

	double range_cruise = 50 * 1000;
	double v_cruise = 55.55555556;
	double t_cruise = range_cruise / v_cruise;
	double rho_cruise = Get_rho(rho, 600);

	double fe_fuse = 1.4 * pow(((GW * 0.453592) / 1000), 2 / 3);		// lb -> kg
	double D_fuselage = fe_fuse * 0.5 * rho * pow(v_cruise, 2);

	double fe_hub = 0.85 * pow(((GW * 0.453592) / 1000), 2 / 3);		// lb -> kg
	double D_hub = fe_hub * 0.5 * rho * pow(v_cruise, 2);

	double alpha = 4 * Pi / 180;
	double Cl = 0.9307;
	double L_wing = 0.5 * rho_cruise * pow(v_cruise, 2) * A_wing * Cl;
	//printf("L_wing = %lf [N]\n", L_wing);
	double Cd = 0.0053;
	double D_wing = 0.5 * rho_cruise * pow(v_cruise, 2) * A_wing * Cd;

	double Tc = D_hub * N_rotor + D_fuselage + D_wing * cos(alpha) + L_wing * sin(alpha);
	double T1 = GW - L_wing * sin(alpha) + D_wing * sin(alpha);

	//printf("Tc = %lf [N]\n", Tc);
	//printf("T1 = %lf [N]\n", T1);

	double P_cruise = (Tc + T1) * v_cruise;
	printf("P_cruise = %lf [W]\n", P_cruise);
	double E_cruise = P_cruise * (t_cruise / 3600);
	//printf("E_cruise = %lf [Wh]\n", E_cruise);


	//double t_des = 9 * 60;
	//double v_des = -1.111111;
	//double Cd0_des = 0.03407;
	//double rho_des = Get_rho(rho, 600);
	//double v_tip_des = v_tip_climb;

	//double P_i_des = Get_P_i(T_rotor, vi+v_des);
	//double P_0_des = Get_P_0(sigma, Cd0_des, rho_des, A_disk, v_tip_des);
	//double P_des = P_i_des + P_0_des;
	//double E_des = P_des * (t_des / 3600);
	//printf("E_des = %lf [Wh]\n", E_des);

	//double E_total = E_climb + E_cruise + E_des;
	//printf("E_total = %lf [Wh]\n", E_total);

	//double capacity_bat = 250;
	//printf("capacity_bat = %lf [Wh/kg]\n", capacity_bat);

	/////////////////
	///// Battery ///
	/////////////////
	//double W_battery = E_total / capacity_bat;
	//printf("W_battery = %lf [kg]\n", W_battery);


	///////////////
	///// Wing ///
	//////////////
	//double S_w = 41.56; // m²
	//double W_fw = W_battery; // fuel weight. 연료가 없으므로 배터리 중량으로 대체.
	//double Rambda = 3.5 * Pi / 180; // 도 -> radian

	//double q = 0.5 * rho_cruise * pow(v_cruise, 2);
	//double rambda = tip_wing / root_wing; // taper ratio.
	//double tc = 0.12044; // wing thickness chord ratio.
	//double N = 2; // The Ultimate Load. 안전계수?
	//double w_dg = GW; // 초기 설계 총중량.
	//double W_wing = 0.036 * pow(S_w, 0.758) * pow(W_fw, 0.0035) * pow(AR / pow(cos(Rambda), 2), 0.6) * pow(q, 0.006) * pow(rambda, 0.04) * pow((100 * tc) / (cos(Rambda)), -0.3) * pow((N * w_dg), 0.49);
	//printf("W_wing = %lf [kg]\n", W_wing);


	//double empty_weight = W_fuse + W_wing + W_rod + W_motor + W_battery + W_landing_gear;
	//printf("empty_weight = %lf [kg]\n", empty_weight);
	//double payload_calc = GW - empty_weight;
	//printf("payload_calc = %lf [kg]\n", payload_calc);
	//return 0;

	// 
		//struct MissionProfile climb;
	//climb.T = GW * g;
	//printf("climb.T = %lf [N]\n", climb.T);
	//climb.T_rotor = climb.T / rotor.N;
	//printf("climb.T_rotor = %lf [N]\n", climb.T_rotor);
	//climb.t = 3 * 60;
	//printf("climb.t = %lf [s]\n", climb.t);
	//climb.h = 600;
	//printf("climb.h = %lf [m]\n", climb.h);
	//climb.rho = Get_rho(rho, climb.h);
	//printf("climb.rho = %lf [kg/s^3]\n", climb.rho);


	//climb.vi = Get_vi(climb.T_rotor, climb.rho, rotor.A_disk);
	//printf("climb.vi = %lf [m/s]\n", climb.vi);
	//climb.vc = 3.333333;
	//printf("climb.vc = %lf [m/s]\n", climb.vc);

	//climb.P_i = Get_P_i(climb.T_rotor, climb.vi + climb.vc);
	//printf("climb.P_i = %lf [W]\n", climb.P_i);
	//climb.P_0 = Get_P_0(rotor.sigma, Cd0, climb.rho, rotor.A_disk, v_tip);
	//printf("climb.P_0 = %lf [W]\n", climb.P_0);
	//climb.P_req = climb.P_i + climb.P_0;
	//printf("climb.P_req = % lf[W]\n", climb.P_req);
	//climb.E = climb.P_req * (climb.t / 3600);
	//printf("climb.E = %lf [Wh]\n", climb.E);


	//rotor.N = 4;
	//rotor.diameter = 3;
	//rotor.radius = 0.5 * rotor.diameter;
	//rotor.A_blade = 0;
	//rotor.A_disk = Pi * pow(rotor.radius, 2);
	////rotor.sigma = rotor.A_blade / rotor.A_disk;
	//rotor.sigma = 0.1696;
	//rotor.rpm = 2000;
	//rotor.omega = Get_omega(rotor.rpm);
	//double v_tip = rotor.omega * (rotor.radius);




	//struct MissionProfile hover;
	//hover.T = GW * g;
	//hover.T_rotor = hover.T / rotor.N;
	//hover.t = 30 * 60;
	//hover.h = 0;
	//hover.rho = Get_rho(rho, hover.h);

	//hover.vi = Get_vi(hover.T_rotor, hover.rho, rotor.A_disk);


	//hover.P_i = Get_P_i(hover.T_rotor, hover.vi);
	//hover.P_0 = Get_P_0(rotor.sigma, Cd0, hover.rho, rotor.A_disk, v_tip);

	//hover.P_req = hover.P_i + hover.P_0;
	//hover.E = hover.P_req * (hover.t / 3600);

	//struct MissionProfile cruise;

	//cruise.v = 42.864538;
	//printf("cruise.v = %lf [m/s]\n", cruise.v);
	//cruise.h = 600;
	//printf("cruise.h = %lf [m]\n", cruise.h);
	//cruise.rho = Get_rho(rho, cruise.h);
	//printf("cruise.rho = %lf [kg/m^3]\n", cruise.rho);

	//double nu = 1.4207E-5;
	//printf("nu = %lf\n", nu);

	//double reynolds = cruise.v * 1.5 / nu;
	//printf("reynolds = %lf\n", reynolds);








	////double L_wing = 0.5 * cruise.rho * pow(cruise.v, 2) * mainwing.A * CL;
	////printf("L_wing = %lf [N]\n", L_wing);



	//while (abs((GW * g) - L_wing) > pow(10, -12))
	//{
	//	if ((GW * g) > L_wing)
	//		mainwing.span = mainwing.span * 1.001;
	//	else
	//		mainwing.span = mainwing.span * 0.999;
	//	//printf("mainwing.span = %lf\n", mainwing.span);
	//	mainwing.A = Get_A_wing(mainwing.root, mainwing.tip, mainwing.span);
	//	//printf("mainwing.A = %lf\n", mainwing.A);
	//	L_wing = 0.5 * cruise.rho * pow(cruise.v, 2) * mainwing.A * CL;
	//	//printf("L_wing = %lf\n", L_wing);
	//}

	//printf("GW*g - L_wing = %lf\n", (GW * g) - L_wing);
	//printf("GW*g = %lf\n", GW * g);
	//printf("L_wing = %lf [N]\n", L_wing);



	//struct MissionProfile descent;
	//descent.T = GW * g;
	//printf("descent.T = %lf [N]\n", descent.T);
	//descent.T_rotor = descent.T / rotor.N;
	//printf("descent.T_rotor = %lf [N]\n", descent.T_rotor);
	//descent.t = 9 * 60;
	//printf("descent.t = %lf [s]\n", descent.t);
	//descent.h = 600;
	//printf("descent.h = %lf [m]\n", descent.h);
	//descent.rho = Get_rho(rho, descent.h);
	//printf("descent.rho = %lf [kg/s^3]\n", descent.rho);


	//descent.vi = Get_vi(descent.T_rotor, descent.rho, rotor.A_disk);
	//printf("descent.vi = %lf [m/s]\n", descent.vi);
	//descent.vc = -1.111111;
	//printf("descent.vc = %lf [m/s]\n", descent.vc);

	//descent.P_i = Get_P_i(descent.T_rotor, descent.vi + descent.vc);
	//printf("descent.P_i = %lf [W]\n", descent.P_i);
	//descent.P_0 = Get_P_0(rotor.sigma, Cd0, descent.rho, rotor.A_disk, v_tip);
	//printf("descent.P_0 = %lf [W]\n", descent.P_0);
	//descent.P_req = descent.P_i + descent.P_0;
	//printf("descent.P_req = % lf[W]\n", descent.P_req);
	//descent.E = descent.P_req * (descent.t / 3600);
	//printf("descent.E = %lf [Wh]\n", descent.E);

	//double E_total = climb.E + cruise.E + descent.E;
	//printf("E_total = %lf [Wh]\n", E_total);




}