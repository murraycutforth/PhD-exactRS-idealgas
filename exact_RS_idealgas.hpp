/*
 *	This file contains a class to solve for the exact solution of the Riemann Problem for the one dimensional Euler
 *	equations with ideal gas equation of state
 *
 *	Author: Murray Cutforth
 *	Date:	06/10/2016
 */



#ifndef EXACT_RS_IDEALGAS_H
#define EXACT_RS_IDEALGAS_H


#include <blitz/array.h>




class exact_rs_idealgas {

	public:

	blitz::Array<double,1> W_L;
	blitz::Array<double,1> W_STAR_L;
	blitz::Array<double,1> W_STAR_R;
	blitz::Array<double,1> W_R;

	double S_L;
	double S_R;
	double S_STAR;
	double S_HL;
	double S_TL;
	double S_HR;
	double S_TR;
	double P_STAR;

	bool left_rarefaction;
	bool left_shock;
	bool right_rarefaction;
	bool right_shock;

	double gamma_L;
	double gamma_R;

	double v_L;
	double v_R;

	

	// Top-level functions which will usually be called outside this class

	exact_rs_idealgas(double gamma_L, double gamma_R);

	void solve_RP (blitz::Array<double,1> W_L, blitz::Array<double,1> W_R);

	blitz::Array<double,1> sample_solution (double S);

	void solve_2D_RP (blitz::Array<double,1> W_L, blitz::Array<double,1> W_R);

	blitz::Array<double,1> sample_2D_solution (double S);


	
	// Functions used to solve for p_star iteratively

	double find_p_star_newtonraphson ();

	double total_pressure_function (double p_star);

	double total_pressure_function_deriv (double p_star);

	double f (double p_star, double rho, double p, double gamma);

	double f_deriv (double p_star, double rho, double p, double gamma);



	// Functions to find the state inside a rarefaction fan

	blitz::Array<double,1> left_rarefaction_fan_state (double S);

	blitz::Array<double,1> right_rarefaction_fan_state (double S);



	// Misc functions

	double Q_K (double p_star, double rho, double p, double gamma);

	

	// Equation of state functions

	double a (double rho, double p, double gamma);

};



#endif
