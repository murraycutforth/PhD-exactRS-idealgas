#include "exact_RS_idealgas.hpp"


#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>



exact_rs_idealgas :: exact_rs_idealgas (double gamma_L, double gamma_R)
:
	W_L	(3),
	W_STAR_L	(3),
	W_STAR_R	(3),
	W_R	(3),
	S_L	(0.0),
	S_R	(0.0),
	S_STAR	(0.0),
	S_HL	(0.0),
	S_TL	(0.0),
	S_HR	(0.0),
	S_TR	(0.0),
	P_STAR	(0.0),
	left_rarefaction	(false),
	left_shock	(false),
	right_rarefaction	(false),
	right_shock	(false),
	gamma_L	(gamma_L),
	gamma_R	(gamma_R)
{}









void exact_rs_idealgas :: solve_RP (blitz::Array<double,1> W_L_in, blitz::Array<double,1> W_R_in)
{
	assert(W_L_in.extent(blitz::firstDim) == 3);
	assert(W_R_in.extent(blitz::firstDim) == 3);
	assert(W_L_in(0) >= 0.0);
	assert(W_L_in(2) >= 0.0);
	assert(W_R_in(0) >= 0.0);
	assert(W_R_in(2) >= 0.0);

	W_L = W_L_in;
	W_R = W_R_in;

	
	// Calculate p_star

	P_STAR = find_p_star_newtonraphson();
	W_STAR_L(2) = P_STAR;
	W_STAR_R(2) = P_STAR;

	
	// Calculate u_star

	S_STAR = 0.5*(W_L(1)+W_R(1)) + 0.5*(f(P_STAR,W_R(0),W_R(2),gamma_R) - f(P_STAR,W_L(0),W_L(2),gamma_L));
	W_STAR_L(1) = S_STAR;
	W_STAR_R(1) = S_STAR;


	// Solution now depends on character of 1st and 3rd waves

	if (P_STAR > W_L(2))
	{
		left_shock = true;
		
		W_STAR_L(0) = W_L(0)*((P_STAR/W_L(2))+((gamma_L-1.0)/(gamma_L+1.0)))/(1.0+((P_STAR/W_L(2))*(gamma_L-1.0)/(gamma_L+1.0)));
		S_L = W_L(1) - (Q_K(P_STAR,W_L(0),W_L(2),gamma_L)/W_L(0));

	}
	else
	{
		left_rarefaction = true;

		W_STAR_L(0) = W_L(0)*std::pow(P_STAR/W_L(2), 1.0/gamma_L);

		double a_L = a(W_L(0), W_L(2), gamma_L);
		double a_star_L = a_L*std::pow(P_STAR/W_L(2), (gamma_L-1.0)/(2.0*gamma_L));

		S_HL = W_L(1) - a_L;
		S_TL = S_STAR - a_star_L;
	}

	if (P_STAR > W_R(2))
	{
		right_shock = true;
		
		W_STAR_R(0) = W_R(0)*((P_STAR/W_R(2))+((gamma_R-1.0)/(gamma_R+1.0)))/(1.0+((P_STAR/W_R(2))*(gamma_R-1.0)/(gamma_R+1.0)));

		S_R = W_R(1) + (Q_K(P_STAR,W_R(0),W_R(2),gamma_R)/W_R(0));
	}
	else
	{
		right_rarefaction = true;

		W_STAR_R(0) = W_R(0)*std::pow(P_STAR/W_R(2), 1.0/gamma_R);

		double a_R = a(W_R(0),W_R(2),gamma_R);
		double a_star_R = a_R*std::pow(P_STAR/W_R(2), (gamma_R-1.0)/(2.0*gamma_R));

		S_HR = W_R(1) + a_R;
		S_TR = S_STAR + a_star_R;
	}


}








blitz::Array<double,1> exact_rs_idealgas :: sample_solution (double S)
{
	blitz::Array<double,1> W (3);

	
	// Find appropriate part of solution

	if (S < S_STAR)
	{
		// To the left of the contact

		if (P_STAR > W_L(2))
		{
			assert(left_shock);
			
			if (S < S_L)
			{
				W = W_L;
			}
			else
			{
				W = W_STAR_L;
			}
		}
		else
		{
			assert(left_rarefaction);
			
			if (S < S_HL)
			{
				W = W_L;
			}
			else
			{
				if (S > S_TL)
				{
					W = W_STAR_L;
				}
				else
				{
					W = left_rarefaction_fan_state(S);
				}
			}
		}
	}
	else
	{
		// To the right of the contact

		if (P_STAR > W_R(2))
		{
			assert(right_shock);
			
			if (S > S_R)
			{
				W = W_R;
			}
			else
			{
				W = W_STAR_R;
			}
		}
		else
		{
			assert(right_rarefaction);
			
			if (S > S_HR)
			{
				W = W_R;
			}
			else
			{
				if (S < S_TR)
				{
					W = W_STAR_R;
				}
				else
				{
					W = right_rarefaction_fan_state(S);
				}
			}
		}
	}

	return W;
}













double exact_rs_idealgas :: find_p_star_newtonraphson ()
{
	double p_star;
	double p_star_next;
	double p_L = W_L(2);
	double p_R = W_R(2);
	double rho_L = W_L(0);
	double rho_R = W_R(0);
	double u_L = W_L(1);
	double u_R = W_R(1);
	double TOL = 1e-6;


	// First we set the initial guess for p_star using the PVRS approximation

	double a_L = a(rho_L, p_L, gamma_L);
	double a_R = a(rho_R, p_R, gamma_R);
	p_star_next = 0.5*(p_L+p_R) - 0.125*(u_R - u_L)*(rho_L + rho_R)*(a_L + a_R);
	p_star_next = std::max(TOL, p_star_next);

	
	// Now we use the Newton-Raphson algorithm

	int num_iter = 0;

	do {
		p_star = p_star_next;

		p_star_next = p_star - total_pressure_function(p_star)/total_pressure_function_deriv(p_star);

		num_iter++;

	} while ((fabs(p_star_next - p_star)/(0.5*(p_star+p_star_next)) > TOL));

	return p_star_next;
}






double exact_rs_idealgas :: total_pressure_function (double p_star)
{
	return	f(p_star, W_L(0), W_L(2), gamma_L)
		+ f(p_star, W_R(0), W_R(2), gamma_R)
		+ W_R(1) - W_L(1);
}




double exact_rs_idealgas :: total_pressure_function_deriv (double p_star)
{
	return 	f_deriv (p_star, W_L(0), W_L(2), gamma_L)
		+ f_deriv (p_star, W_R(0), W_R(2), gamma_R);
}






double exact_rs_idealgas :: f (double p_star, double rho, double p, double gamma)
{
	double A = 2.0/((gamma+1.0)*rho);
	double B = p*(gamma-1.0)/(gamma+1.0);

	if (p_star > p)
	{
		return (p_star - p)*sqrt(A/(p_star + B));
	}
	else
	{
		return (2.0*a(rho,p,gamma)/(gamma-1.0))*(std::pow(p_star/p, (gamma-1.0)/(2.0*gamma)) - 1.0);
	}
}




double exact_rs_idealgas :: f_deriv (double p_star, double rho, double p, double gamma)
{
	double A = 2.0/((gamma+1.0)*rho);
	double B = p*(gamma-1.0)/(gamma+1.0);

	if (p_star > p)
	{
		return sqrt(A/(B+p_star))*(1.0 - ((p_star-p)/(2.0*(B+p_star))));
	}
	else
	{
		return (1.0/(rho*a(rho,p,gamma)))*std::pow(p_star/p, -(gamma+1.0)/(2.0*gamma));
	}
}














blitz::Array<double,1> exact_rs_idealgas :: left_rarefaction_fan_state (double S)
{
	assert(left_rarefaction == true);

	double a_L = a(W_L(0),W_L(2),gamma_L);
	double rho = W_L(0)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), 2.0/(gamma_L - 1.0));
	double u = (2.0/(gamma_L+1.0))*(a_L + S + ((gamma_L-1.0)/2.0)*W_L(1));
	double p = W_L(2)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), (2.0*gamma_L)/(gamma_L-1.0));

	blitz::Array<double,1> result (3);
	result(0) = rho;
	result(1) = u;
	result(2) = p;

	return result;
}








blitz::Array<double,1> exact_rs_idealgas :: right_rarefaction_fan_state (double S)
{
	assert(right_rarefaction == true);

	double a_R = a(W_R(0),W_R(2),gamma_R);
	double rho = W_R(0)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), 2.0/(gamma_R - 1.0));
	double u = (2.0/(gamma_R+1.0))*(- a_R + S + ((gamma_R-1.0)/2.0)*W_R(1));
	double p = W_R(2)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), (2.0*gamma_R)/(gamma_R-1.0));

	blitz::Array<double,1> result (3);
	result(0) = rho;
	result(1) = u;
	result(2) = p;

	return result;
}




















double exact_rs_idealgas :: Q_K (double p_star, double rho, double p, double gamma)
{
	double A = 2.0/((gamma+1.0)*rho);
	double B = p*(gamma-1.0)/(gamma+1.0);
	return sqrt((p_star+B)/A);
}
	

















double exact_rs_idealgas :: a (double rho, double p, double gamma)
{
	return sqrt(gamma*(p/rho));
}





