#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "/usr/local/cpp/src/stdafx.h"
#include "/usr/local/cpp/src/optimization.h"
#include "/usr/local/cpp/src/ap.h"
#include "/usr/local/cpp/src/ap.cpp"
#include "/usr/local/cpp/src/alglibinternal.cpp"
#include "/usr/local/cpp/src/alglibmisc.cpp"
#include "/usr/local/cpp/src/linalg.cpp"
#include "/usr/local/cpp/src/solvers.cpp"
#include "/usr/local/cpp/src/specialfunctions.cpp"
#include "/usr/local/cpp/src/integration.cpp"
#include "/usr/local/cpp/src/optimization.cpp"

using namespace alglib;

typedef std::vector <double> Vec;
typedef std::vector <Vec> Mat;


//svr for function f(x) = (6x-2)^2 * sin (12x-4)

double findFunctionValue (double x)
{
    double z = (6*x - 2)*(6*x - 2);
    z *= sin (12*x - 4);
    return z;
}

int main ()
{
	const int noOfDataPoints = 21;
	double X [noOfDataPoints], y[noOfDataPoints];
	double k = 0;

	for (int i = 0; i < noOfDataPoints; i++)
	{
    	 	X[i] = k;
    		y[i] = findFunctionValue (k);
  		k+= 0.05;
	}
	/*printf ("X values\n");
	for (int i = 0; i < noOfDataPoints; i++)
		printf ("%lf\n", X[i]);
	printf ("Y values\n");
	for (int i = 0; i < noOfDataPoints; i++)
		printf ("%lf\n", y[i]);*/
	

//Build correlation matrix

	double sigma = 0.2;
	double psi [noOfDataPoints][noOfDataPoints];

	for (int i = 0; i < noOfDataPoints; i++)
    		for (int j = 0; j < noOfDataPoints; j++)
        		psi[i][j] = exp ((-1* pow(abs(X[i] - X[j]), 2))/(pow(sigma,2)));

	int e = 1;
	int C = 1000;
	double xi = pow(10, -6);

	Mat Psi;
	Vec r1 (2*noOfDataPoints);
	int i = 0, j = 0;

	for (int i = 0; i < noOfDataPoints; i++)
	{
		for (int j = 0; j < noOfDataPoints; j++)
		{
			r1[j] = psi[i][j];
			r1[j+noOfDataPoints] = -1*psi[i][j];
		}
		Psi.push_back (r1);
	}
	for (int i = noOfDataPoints; i < 2*noOfDataPoints; i++)
	{
		for (int j = 0; j < noOfDataPoints; j++)
		{
			r1[j] = -1*psi[i-noOfDataPoints][j];
			r1[j+noOfDataPoints] = psi[i-noOfDataPoints][j];
		}
		Psi.push_back(r1);
	}
	/*for (int i = 0; i < 2*noOfDataPoints; i++)
	{	
		printf ("row %d\n", i);		
		for (int j = 0; j < 2*noOfDataPoints; j++)
		{
					
			printf ("%lf ", Psi[i][j]);
		}
		printf ("\n");
	}*/

	//alglib for quadratic programming

	double psiarr [1764];
	int pos = 0;
	for (int i = 0; i < 2*noOfDataPoints; i++)
	{
		for (int j = 0; j < 2*noOfDataPoints; j++)
		{
			psiarr[pos] = Psi[i][j];
			pos++;
		}
	}	


	real_2d_array A;
	A.setcontent (42, 42, psiarr);
	//printf ("printing a\n");
	//for (int i = 0; i < 42; i++)
	//	for (int j = 0; j < 42; j++)
	//		printf ("%lf\n", A[i][j]);

	double p [2*noOfDataPoints];
	for (int i = 0; i < noOfDataPoints; i++)
	{
		p[i] = e - y[i];
		p[i+noOfDataPoints] = e + y[i];
	}
	real_1d_array B;
	B.setcontent (42, p);
	//printf ("printing b\n");
	//for (int i = 0; i < 42; i++)
	//	printf ("%lf\n", B[i]);

	double sarr [2*noOfDataPoints] = {1};
	real_1d_array s;
	s.setcontent (42, sarr);
	

	double constraints [43];
	for (int i = 0; i < 21; i++)
		constraints [i] = -1;
	for (int i = 21; i < 42; i++)
		constraints [i] = 1;
	constraints [42] = 0;	
	real_2d_array Con; //constraints
	Con.setcontent (1, 43, constraints);
	//printf ("printing constraints\n");
	//for (int i = 0; i < 43; i++)
	//	printf ("%lf\n", Con[0][i]);

	alglib::integer_1d_array ct = "[0]";

	double ubval = C/noOfDataPoints;

	double lowbnd [2*noOfDataPoints] = {0};
	real_1d_array bndl;
	bndl.setcontent (42, lowbnd);
	//printf ("printing low bounds\n");
	//for (int i = 0; i < 42; i++)
	//	printf ("%lf\n", bndl[i]);

	
	double uppbnd [2*noOfDataPoints];
	for (int i = 0; i < 2*noOfDataPoints; i++)
		uppbnd[i] = ubval;
	real_1d_array bndu;
	bndu.setcontent (42, uppbnd);
	
	//printf ("printing upp bounds\n");
	//for (int i = 0; i < 42; i++)
	//	printf ("%lf\n", bndu[i]);

	real_1d_array x; //to store results

	minqpstate state;
	minqpreport rep;

	printf ("Variables made\n");

	minqpcreate (noOfDataPoints*2, state);
	minqpsetquadraticterm (state, A);
	minqpsetlinearterm (state, B);
	minqpsetbc (state, bndl, bndu);
	minqpsetlc (state, Con, ct);
	minqpsetscaleautodiag (state);
	minqpsetalgobleic (state, 0.0, 0.0, 0.0, 0);
	minqpoptimize (state);
	minqpresults (state, x, rep);
	printf ("Printing\n");
	printf("%s\n", x.tostring(45).c_str());
	printf ("report\n");
	printf ("%ld\n", rep.terminationtype);

	

	double alpha_pm [noOfDataPoints];
	for (int i = 0; i < noOfDataPoints; i++)
	{
		alpha_pm[i] = x[i] - x[i+noOfDataPoints];
		std::cout << alpha_pm[i];
		printf ("\n");
	}

	//find indices of SVs
	Vec sv_i;
	for (int i = 0; i < noOfDataPoints; i++)
	{
		if (abs(alpha_pm[i]) > xi)
			sv_i.push_back(i);
	}
	printf ("sv_i");
	for (int i = 0; i < sv_i.size(); i++)
		printf ("%lf ", sv_i[i]);

	//find SV midway between 0 and C
	int sv_mid_i = 0;
	double sub = (double)C/(2*noOfDataPoints);
	printf ("c/2n is %lf\n", sub);
	double sv_mid = abs(abs(alpha_pm[0]) - sub);
	for (int i = 1; i < noOfDataPoints; i++)
	{
		
		double val = abs(abs(alpha_pm[i])-sub);
		printf ("printing val %lf\n", val);
		if (val < sv_mid)
		{
			sv_mid = val;
			sv_mid_i = i;
		}
	}
	printf ("Sv mid %lf\n", sv_mid);
	printf ("sv mid i %d\n", sv_mid_i);

	//find mu
	int sign;
	if (alpha_pm[sv_mid_i] > 0)
		sign = 1;
	else if (alpha_pm[sv_mid_i] < 0)
		sign = -1;
	else
		sign = 0;
	double mu = y[sv_mid_i] - (e*sign);
	for (int i = 0; i < sv_i.size(); i++)
	{
		int ind = sv_i[i];		
		mu -= (alpha_pm[ind]*Psi[ind][sv_mid_i]);
	}	
	printf ("mu is %lf\n", mu);

	//points at which to plot prediction
	double xval [101];
	double psipred [noOfDataPoints][1];
	k = 0;
	for (int i = 0; i < 101; i++)
	{
		xval[i] = k;
		k += 0.01;
	}
	double pred [noOfDataPoints];
	for (int i = 0; i < 101; i++)
	{
		for (int j = 0; j < noOfDataPoints; j++)
			psipred[j][0] = exp ((-1* pow((xval[i] - X[j]), 2))/(pow(sigma,2)));
		pred[i] = mu;
		for (int j = 0; j < noOfDataPoints; j++)		
			pred[i] += alpha_pm[j] * psipred[j][0];
		printf ("pred %d is %lf as against %lf\n", i, pred[i], y[i]);
	}

	return 0;
}
