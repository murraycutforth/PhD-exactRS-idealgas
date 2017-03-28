#include "exact_RS_idealgas.hpp"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>



exact_rs_idealgas :: exact_rs_idealgas (double gamma_L, double gamma_R)
:
	gamma_L	(gamma_L),
	gamma_R	(gamma_R),
	S_STAR	(0.0),
	P_STAR	(0.0),
	rho_star_L	(0.0),
	rho_star_R	(0.0),
	S_L	(0.0),
	S_R	(0.0),
	S_HL	(0.0),
	S_TL	(0.0),
	S_HR	(0.0),
	S_TR	(0.0)
{}









void exact_rs_idealgas :: solve_RP (const blitz::Array<double,1>& W_L, const blitz::Array<double,1>& W_R)
{
	assert(W_L(0) >= 0.0);
	assert(W_L(2) >= 0.0);
	assert(W_R(0) >= 0.0);
	assert(W_R(2) >= 0.0);

	
	// Calculate p_star

	P_STAR = find_p_star_newtonraphson(W_L(0), W_L(1), W_L(2), W_R(0), W_R(1), W_R(2));
	
	
	// Calculate u_star

	S_STAR = 0.5*(W_L(1)+W_R(1)) + 0.5*(f(P_STAR,W_R(0),W_R(2),gamma_R) - f(P_STAR,W_L(0),W_L(2),gamma_L));


	// Solution now depends on character of 1st and 3rd waves

	if (P_STAR > W_L(2))
	{
		// Left shock
		
		rho_star_L = W_L(0)*((P_STAR/W_L(2))+((gamma_L-1.0)/(gamma_L+1.0)))/(1.0+((P_STAR/W_L(2))*(gamma_L-1.0)/(gamma_L+1.0)));
		S_L = W_L(1) - (Q_K(P_STAR,W_L(0),W_L(2),gamma_L)/W_L(0));

	}
	else
	{
		// Left rarefaction

		rho_star_L = W_L(0)*std::pow(P_STAR/W_L(2), 1.0/gamma_L);

		double a_L = a(W_L(0), W_L(2), gamma_L);
		double a_star_L = a_L*std::pow(P_STAR/W_L(2), (gamma_L-1.0)/(2.0*gamma_L));

		S_HL = W_L(1) - a_L;
		S_TL = S_STAR - a_star_L;
	}

	if (P_STAR > W_R(2))
	{
		// Right shock
				
		rho_star_R = W_R(0)*((P_STAR/W_R(2))+((gamma_R-1.0)/(gamma_R+1.0)))/(1.0+((P_STAR/W_R(2))*(gamma_R-1.0)/(gamma_R+1.0)));

		S_R = W_R(1) + (Q_K(P_STAR,W_R(0),W_R(2),gamma_R)/W_R(0));
	}
	else
	{
		// Right rarefaction

		rho_star_R = W_R(0)*std::pow(P_STAR/W_R(2), 1.0/gamma_R);

		double a_R = a(W_R(0),W_R(2),gamma_R);
		double a_star_R = a_R*std::pow(P_STAR/W_R(2), (gamma_R-1.0)/(2.0*gamma_R));

		S_HR = W_R(1) + a_R;
		S_TR = S_STAR + a_star_R;
	}
}





blitz::Array<double,1> exact_rs_idealgas :: sample_solution (const blitz::Array<double,1>& W_L, const blitz::Array<double,1>& W_R, double S)
{
	blitz::Array<double,1> W (3);

	
	// Find appropriate part of solution and return primitives

	if (S < S_STAR)
	{
		// To the left of the contact

		if (P_STAR > W_L(2))
		{
			// Left shock
			
			if (S < S_L)
			{
				W = W_L;
			}
			else
			{
				W(0) = rho_star_L;
				W(1) = S_STAR;
				W(2) = P_STAR;
			}
		}
		else
		{
			// Left rarefaction
			
			if (S < S_HL)
			{
				W = W_L;
			}
			else
			{
				if (S > S_TL)
				{
					W(0) = rho_star_L;
					W(1) = S_STAR;
					W(2) = P_STAR;
				}
				else
				{
					set_left_rarefaction_fan_state(W_L, S, W);
				}
			}
		}
	}
	else
	{
		// To the right of the contact

		if (P_STAR > W_R(2))
		{
			// Right shock
			
			if (S > S_R)
			{
				W = W_R;
			}
			else
			{
				W(0) = rho_star_R;
				W(1) = S_STAR;
				W(2) = P_STAR;
			}
		}
		else
		{
			// Right rarefaction
			
			if (S > S_HR)
			{
				W = W_R;
			}
			else
			{
				if (S < S_TR)
				{
					W(0) = rho_star_R;
					W(1) = S_STAR;
					W(2) = P_STAR;
				}
				else
				{
					set_right_rarefaction_fan_state(W_R, S, W);
				}
			}
		}
	}

	return W;
}



void exact_rs_idealgas :: celledge_primitives_2D (
		
	const blitz::Array<double,1>& W_L, 
	const blitz::Array<double,1>& W_R,
	blitz::Array<double,1>& soln
)
{
	const double v_L = W_L(2);
	const double v_R = W_R(2);
	
	
	// Calculate p_star

	P_STAR = find_p_star_newtonraphson(W_L(0), W_L(1), W_L(3), W_R(0), W_R(1), W_R(3));
	
	
	// Calculate u_star

	S_STAR = 0.5*(W_L(1)+W_R(1)) + 0.5*(f(P_STAR,W_R(0),W_R(3),gamma_R) - f(P_STAR,W_L(0),W_L(3),gamma_L));
	
	
	if (S_STAR > 0.0)
	{
		// Look at states to the left of the contact
		
		soln(2) = v_L;
		
		if (P_STAR > W_L(3))
		{
			// Left shock
			
			S_L = W_L(1) - (Q_K(P_STAR,W_L(0),W_L(3),gamma_L)/W_L(0));
			
			if (S_L > 0.0)
			{
				soln = W_L;
			}
			else
			{
				soln(0) = W_L(0)*((P_STAR/W_L(3))+((gamma_L-1.0)/(gamma_L+1.0)))/(1.0+((P_STAR/W_L(3))*(gamma_L-1.0)/(gamma_L+1.0)));
				soln(1) = S_STAR;
				soln(3) = P_STAR;
			}
		}
		else
		{
			// Left rarefaction
			
			double a_L = a(W_L(0), W_L(3), gamma_L);
			double a_star_L = a_L*std::pow(P_STAR/W_L(3), (gamma_L-1.0)/(2.0*gamma_L));
			S_HL = W_L(1) - a_L;
			S_TL = S_STAR - a_star_L;
			
			if (S_HL > 0.0)
			{
				soln = W_L;
			}
			else
			{
				if (S_TL < 0.0)
				{
					soln(0) = W_L(0)*std::pow(P_STAR/W_L(3), 1.0/gamma_L);
					soln(1) = S_STAR;
					soln(3) = P_STAR;
				}
				else
				{
					double a_L = a(W_L(0),W_L(3),gamma_L);
					soln(0) = W_L(0)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1)), 2.0/(gamma_L - 1.0));
					soln(1) = (2.0/(gamma_L+1.0))*(a_L + ((gamma_L-1.0)/2.0)*W_L(1));
					soln(3) = W_L(3)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1)), (2.0*gamma_L)/(gamma_L-1.0));
				}
			}
		}
	}
	else
	{
		// Look at states to the right of the contact
		
		soln(2) = v_R;
		
		if (P_STAR > W_R(3))
		{
			// Right shock
			
			S_R = W_R(1) + (Q_K(P_STAR,W_R(0),W_R(3),gamma_R)/W_R(0));
			
			if (S_R < 0.0)
			{
				soln = W_R;
			}
			else
			{
				soln(0) = W_R(0)*((P_STAR/W_R(3))+((gamma_R-1.0)/(gamma_R+1.0)))/(1.0+((P_STAR/W_R(3))*(gamma_R-1.0)/(gamma_R+1.0)));
				soln(1) = S_STAR;
				soln(3) = P_STAR;
			}
		}
		else
		{
			// Right rarefaction
			
			double a_R = a(W_R(0),W_R(3),gamma_R);
			double a_star_R = a_R*std::pow(P_STAR/W_R(3), (gamma_R-1.0)/(2.0*gamma_R));
			S_HR = W_R(1) + a_R;
			S_TR = S_STAR + a_star_R;
			
			if (S_HR < 0.0)
			{
				soln = W_R;
			}
			else
			{
				if (S_TL > 0.0)
				{
					soln(0) = W_R(0)*std::pow(P_STAR/W_R(3), 1.0/gamma_R);
					soln(1) = S_STAR;
					soln(3) = P_STAR;
				}
				else
				{
					double a_R = a(W_R(0),W_R(3),gamma_R);

					soln(0) = W_R(0)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1)), 2.0/(gamma_R - 1.0)); 
					soln(1) = (2.0/(gamma_R+1.0))*(- a_R + ((gamma_R-1.0)/2.0)*W_R(1));
					soln(3) = W_R(3)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1)), (2.0*gamma_R)/(gamma_R-1.0));
				}
			}
		}
	}
}











double exact_rs_idealgas :: find_p_star_newtonraphson (
	
	const double rho_L,
	const double u_L,
	const double p_L,
	const double rho_R,
	const double u_R,
	const double p_R
)
{
	double p_star;
	double p_star_next;
	const double TOL = 1e-6;


	// First we set the initial guess for p_star using the PVRS approximation

	const double a_L = a(rho_L, p_L, gamma_L);
	const double a_R = a(rho_R, p_R, gamma_R);
	p_star_next = 0.5*(p_L+p_R) - 0.125*(u_R - u_L)*(rho_L + rho_R)*(a_L + a_R);
	p_star_next = std::max(TOL, p_star_next);

	
	// Now use the Newton-Raphson algorithm

	do {
		p_star = p_star_next;

		p_star_next = p_star - total_pressure_function(p_star,rho_L,u_L,p_L,rho_R,u_R,p_R)/total_pressure_function_deriv(p_star,rho_L,p_L,rho_R,p_R);

	} while ((fabs(p_star_next - p_star)/(0.5*(p_star+p_star_next)) > TOL));

	return p_star_next;
}






double exact_rs_idealgas :: total_pressure_function (

	const double p_star,
	const double rho_L,
	const double u_L,
	const double p_L,
	const double rho_R,
	const double u_R,
	const double p_R
)
{
	return	f(p_star, rho_L, p_L, gamma_L)
		+ f(p_star, rho_R, p_R, gamma_R)
		+ u_R - u_L;
}




double exact_rs_idealgas :: total_pressure_function_deriv (

	const double p_star,
	const double rho_L,
	const double p_L,
	const double rho_R,
	const double p_R
)
{
	return 	f_deriv (p_star, rho_L, p_L, gamma_L)
		+ f_deriv (p_star, rho_R, p_R, gamma_R);
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














void exact_rs_idealgas :: set_left_rarefaction_fan_state (const blitz::Array<double,1>& W_L, double S, blitz::Array<double,1>& W)
{
	double a_L = a(W_L(0),W_L(2),gamma_L);
	W(0) = W_L(0)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), 2.0/(gamma_L - 1.0));
	W(1) = (2.0/(gamma_L+1.0))*(a_L + S + ((gamma_L-1.0)/2.0)*W_L(1));
	W(2) = W_L(2)*std::pow((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(W_L(1) - S), (2.0*gamma_L)/(gamma_L-1.0));
}





void exact_rs_idealgas :: set_right_rarefaction_fan_state (const blitz::Array<double,1>& W_R, double S, blitz::Array<double,1>& W)
{
	double a_R = a(W_R(0),W_R(2),gamma_R);
	W(0) = W_R(0)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), 2.0/(gamma_R - 1.0));
	W(1) = (2.0/(gamma_R+1.0))*(- a_R + S + ((gamma_R-1.0)/2.0)*W_R(1));
	W(2) = W_R(2)*std::pow((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(W_R(1) - S), (2.0*gamma_R)/(gamma_R-1.0));
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





