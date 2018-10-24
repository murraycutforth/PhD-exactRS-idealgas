/*
 *	This file contains a class to solve for the exact solution of the Riemann Problem for the one dimensional Euler
 *	equations with ideal gas equation of state
 *
 *	Author: Murray Cutforth
 *	Date:	06/10/2016
 */

#ifndef EXACT_RS_IDEALGAS_H
#define EXACT_RS_IDEALGAS_H

#include <Eigen/Dense>

class exact_rs_idealgas {

	public:
	
	const double gamma_L;
	const double gamma_R;
	
	double S_STAR;
	double P_STAR;
	double rho_star_L;
	double rho_star_R;
	
	double S_L;
	double S_R;
	double S_HL;
	double S_TL;
	double S_HR;
	double S_TR;

	exact_rs_idealgas (double gamma_L, double gamma_R);
	
	
	// Functions used to generate exact solutions to Riemann problems

	void solve_RP (const Eigen::Vector3d& W_L, const Eigen::Vector3d& W_R);

	Eigen::Vector3d sample_solution (const Eigen::Vector3d& W_L, const Eigen::Vector3d& W_R, double S);
	
	
	// Function called from 2D Godunov-type methods
	
	void celledge_primitives_2D (
		
		const Eigen::Vector4d& W_L, 
		const Eigen::Vector4d& W_R,
		Eigen::Vector4d& soln
	);


	// Functions used to solve for p_star iteratively

	double find_p_star_newtonraphson (
	
		const double rho_L,
		const double u_L,
		const double p_L,
		const double rho_R,
		const double u_R,
		const double p_R
	);

	double total_pressure_function (

		const double p_star,
		const double rho_L,
		const double u_L,
		const double p_L,
		const double rho_R,
		const double u_R,
		const double p_R
	);

	double total_pressure_function_deriv (

		const double p_star,
		const double rho_L,
		const double p_L,
		const double rho_R,
		const double p_R
	);

	double f (double p_star, double rho, double p, double gamma);

	double f_deriv (double p_star, double rho, double p, double gamma);


	// Functions to find the state inside a rarefaction fan

	void set_left_rarefaction_fan_state (const Eigen::Vector3d& W_L, double S, Eigen::Vector3d& W);

	void set_right_rarefaction_fan_state (const Eigen::Vector3d& W_R, double S, Eigen::Vector3d& W);


	// Misc functions

	double Q_K (double p_star, double rho, double p, double gamma);

	
	// Equation of state functions

	double a (double rho, double p, double gamma);
};

#endif
