/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

/** @file
 * @brief Computational routines for calculating NSB entropy.
 * This file contains the computational routines for calculating the entropy of
 * a random variable from a hist1d vector of word counts using the Nemenman-
 * Shafee-Bialek (NSB) method.
 */

#include "toolkit_c.h"
#include <stdexcept>
#include <sstream>
#include <vector>

#ifdef DEBUG
#include <iostream>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

/**
 * @brief Contains the members, classes, and functions used to calculate NSB entropy.
 * This namespace contains all the members, classes, and functions used to calculate
 * NSB entropy. The routines entropy_nsb and variance_nsb, which follow the format
 * for the toolkit, uses objects in this namespace to perform the actual calculations
 * of NSB entropy and its variance. Where possible, variable names follow the
 * notation used in "Inference of Entropies of Discrete Random Variables with Unknown
 * Cardinalities" by Ilya Nemenman (arXiv:physics/0207009v1), and references to
 * particular equations are noted. Other variable names are kept from the original
 * nsb-entropy code from which this code was distilled
 * (http://nsb-entropy.sourceforge.net/).
 */
namespace nsb_entropy
{
	/**
	 * NSB computational warnings values.
	 */
	enum NSB_WARN
	{
		NSB_OK,
		NSB_NOCOINC,
		NSB_NOEVENTS,
		NSB_ALLCOINC,
		NSB_SER_B0,
		NSB_NR_B0_SIGN,
		NSB_NR_B0,
		NSB_NR_BCL_SIGN,
		NSB_NR_BCL,
		NSB_INT_NOCONV,
	};

	/**
	 * NSB computational warnings names.
	 */
	const char *p_warn_names[] = 
	{
		"All OK. Recovered.",
		"No coincidences.",
		"Number of events is less than or equal to one.",
		"All data coincide.",
		"Series expansion for B0 did not converge. Will probably recover.",
		"Newton-Raphson ERROR in B0 calculation: wrong sign. May recover.",
		"Newton-Raphson search for B0 did not converge. Will probably recover.",
		"Newton-Raphson ERROR in Bcl calculation: wrong sign. Possibly serious.",
		"Newton-Raphson search for Bcl did not converge.",
		"Numerical evaluation of an integral did not converge.",
	};

	const double psi0_1 = gsl_sf_psi(1); //additive inverse of Euler's constant (C_gamma, see equation 26)

	//forward declarations
	double g_int_1(double, void*);
	double g_int_S(double, void*);
	double g_int_S2(double, void*);

	//polygamma functions
	inline double psi(double x) {return gsl_sf_psi(x);}
	inline double psi1(double x) {return gsl_sf_psi_1(x);}
	inline double psi2(double x) {return gsl_sf_psi_n(2,x);}
	inline double psi3(double x) {return gsl_sf_psi_n(3,x);}
	inline double psi4(double x) {return gsl_sf_psi_n(4,x);}

	/**
	 * @brief Factorial fraction.
	 * Given b[0]=1/a, we have b[i]=b[i-1]/(i+a).
	 */
	inline void factorial_fraction(const int size, const double &a, double *b)
	{
		if(size>0)
		{
			b[0] = 1/a;
			for(int i=1; i<size; i++)
				b[i] = b[i-1]/(i+a);
		}
	}

	/**
	 * @brief Sum a power series.
	 * We always sum a series backwards, since the last terms are usually
	 * smaller, and this reduces the loss in precision. Sum a[i]*x^(i+n),
	 * where n is an integer.
	 */
	inline const double series(const int size, const double *a, const double &x, const int n, double result=0.0)
	{
		for(int i=size-1; i>=0; i--)
			result += a[i]*gsl_pow_int(x,i+n);
		return result;
	}

	/**
	 * @brief Sum a power series of a function.
	 * We always sum a series backwards, since the last terms are usually
	 * smaller, and this reduces the loss in precision. Sum a[i]*x^(i+n)*(*f)(i+m,y),
	 * where n is an integer, and y is a double.
	 */
	inline const double series_f(const int size, const double *a, const double &x, const int n, double (*f)(int, double), const double &y, const int m, double result=0.0)
	{
		for(int i=size-1; i>=0; i--)
			result += a[i]*gsl_pow_int(x,i+n)*((*f)(i+m,y));
		return result;
	}

	/**
	 * @brief Calculate the dot product f(a)*b.
	 * For vectors a and b of doubles, and function f that accepts a double,
	 * calculate f(a)*b in Matlab/Octave vector notation. The return value is
	 * added to the last argument, and the size of the arrays is the first
	 * argument.
	 */
	inline const double dot_f_a_b(const int size, const double *a, const double *b, double (*f)(double), double result=0.0)
	{
		for(int i=0; i<size; i++)
			result += (*f)(a[i])*b[i];
		return result;
	}

	/**
	 * @brief Calculate the dot product f(a+a0)*b.
	 * For double a0, vectors a and b of doubles, and function f that accepts a
	 * double, calculate f(a+a0)*b in Matlab/Octave vector notation. The return
	 * value is added to the last argument, and the size of the arrays is the
	 * first argument.
	 */
	inline const double dot_f_aplusa0_b(const int size, const double *a, const double *b, const double a0, double (*f)(double), double result=0.0)
	{
		for(int i=0; i<size; i++)
			result += (*f)(a[i]+a0)*b[i];
		return result;
	}

	/**
	 * @brief Calculate the dot product (f(a+a0).*(b+b0))'*c.
	 * For doubles a0 and b0, vectors a, b, and c of doubles, and function f that
	 * accepts a double, calculate (f(a+a0).*(b+b0))'*c in Matlab/Octave vector
	 * notation. The return value is added to the last argument, and the size of
	 * the arrays is the first argument.
	 */
	inline const double dot_f_aplusa0_bplusb0_c(const int size, const double *a, const double *b, const double *c, const double a0, const double b0, double (*f)(double), double result=0.0)
	{
		for(int i=0; i<size; i++)
			result += (*f)(a[i]+a0)*(b[i]+b0)*c[i];
		return result;
	}

	/**
	 * Subleading asymptotic behavior of Psi (or Digamma) function valid
	 * for positive argument. We aim at extremely large x, so we use the
	 * Stirling form for psi, and not the Lancocs one. The Stirling form
	 * has the log(x) term in it, and subtraction of the logarithm can be
	 * done with no loss of precision.
	 */
	inline const double psi_asymp(const double &x)
	{
		const double b[]= { //even Bernoulli coefficients (b_1=-0.5 not shown here)
			0.166666666666667, -0.0333333333333333, 0.0238095238095238, -0.0333333333333333,
			0.0757575757575758, -0.253113553113553, 1.16666666666667, -7.09215686274510,
			54.9711779448622, -529.124242424242, 6192.12318840580, -86580.2531135531,
			1425517.16666667, -27298231.0678161, 601580873.900642, -15116315767.0922,
			429614643061.167, -13711655205088.3, 4.88332318973593e+14, -19296579341940068.0};

		const int bn = 20; //length of the above array
		const double asx = 10.0; //value of x beyond which asymptotic is believed to work
		const double xx = GSL_MAX(x, asx+ (x-floor(x))); //the value to go into the asymptotic formula
		const long recur = lround(GSL_MAX(ceil(asx-x), 0.0));	//number of recursions needed to get to that value

		double f = -series(bn, b, gsl_pow_2(1/xx), 1);
		f -= 0.5/xx; //adding the only odd order term

		//accounting for recursion that brought x to asymptotic regime
		if(recur>0)
			for(long i=recur-1; i>=0; i--) //subtract off smallest value first to improve precision
				f -= 1/(x+i);

		return f;
	}

	/**
	 * Inverse of the error function in double precision.
	 */
	inline const double dierfc(const double y)
	{
		double s, t, u, w, x, z;
	
		z = y;
		if(y>1)
			z = 2 - y;
		w = 0.916461398268964 - log(z);
		u = sqrt(w);
		s = (log(u) + 0.488826640273108) / w;
		t = 1.0 / (u + 0.231729200323405);
		x = u * (1.0 - s * (s * 0.124610454613712 + 0.5)) - 
			((((-0.0728846765585675 * t + 0.269999308670029) * t + 
			0.150689047360223) * t + 0.116065025341614) * t + 
			0.499999303439796) * t;
		t = 3.97886080735226 / (x + 3.97886080735226);
		u = t - 0.5;
		s = (((((((((0.00112648096188977922 * u + 
			1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
			7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
			0.00339721910367775861) * u - 0.011274916933250487) * u - 
			0.0118598117047771104) * u + 0.0142961988697898018) * u + 
			0.0346494207789099922) * u + 0.00220995927012179067;
		s = ((((((((((((s * u - 0.0743424357241784861) * u - 
			0.105872177941595488) * u + 0.0147297938331485121) * u + 
			0.316847638520135944) * u + 0.713657635868730364) * u + 
			1.05375024970847138) * u + 1.21448730779995237) * u + 
			1.16374581931560831) * u + 0.956464974744799006) * u + 
			0.686265948274097816) * u + 0.434397492331430115) * u + 
			0.244044510593190935) * t - 
			z * exp(x * x - 0.120782237635245222);
		x += s * (x * s + 1.0);
		if(y>1.0)
			x = -x;
		return x;
	}

	/**
	 * This is a class to encapsulate a GSL spline.
	 */
	class Spline
	{
	private:
		gsl_spline *p_spline;
		gsl_interp_accel *p_search_states;
		const int size;
	public:
		const double eval(const double &xnew) const {return gsl_spline_eval(p_spline,xnew,p_search_states);}

		//use linear interpolation to conserve memory
		Spline(const double *x, const double *y, const int points, const gsl_interp_type *T=gsl_interp_linear): size(points)
		{
			int status;

			p_spline = gsl_spline_alloc(T,size);
			if(p_spline==NULL)
				throw std::bad_alloc();
			status = gsl_spline_init(p_spline,x,y,size);
			if(status)
				throw std::bad_alloc();
			p_search_states = gsl_interp_accel_alloc();
			if(p_search_states==NULL)
				throw std::bad_alloc();
		}

		~Spline()
		{
			gsl_spline_free(p_spline);
			gsl_interp_accel_free(p_search_states);
		}
	};

	/**
	 * This is a class to encapsulate the GSL integration workspace.
	 */
	class Workspace
	{
	private:
		gsl_integration_workspace *workspace;
		const int size;
	public:
		const int get_size() const {return size;}
		gsl_integration_workspace* get_workspace() const {return workspace;}

		Workspace(const int intervals=10000): workspace(NULL), size(intervals)
		{
			workspace = gsl_integration_workspace_alloc(intervals);
			if(workspace==NULL)
				throw std::bad_alloc();
		}

		~Workspace() {gsl_integration_workspace_free(workspace);}
	};

	/**
	 * This is a class to encapsulate the interpolation splines. These splines
	 * relate the Dirichlet parameter beta (kappa=K*beta) to the a priori entropy
	 * xi. Unfortunately, the rationale for how this is done is not made explicit
	 * in Nemenman (2002), thus we simply accept here that it "works", as
	 * evidenced by direct comparison to the original NSB entropy code, and based
	 * on the following (edited) comments from Nemenman. "This function is a lot
	 * different from that in the Octave code. We split the range of kappa into
	 * three asymptotic regimes (kappa<=1, 1<kappa<=K, kappa>K), and use the
	 * asymptotic expansion of kappa (maybe equation 12) in each of the regimes
	 * as an equispaced variable (easily analytically related to kappa). We then
	 * build three different (overlapped for better precision) splines to relate
	 * xi to this asymptotic variable (and thus to kappa) in each regime."
	 */
	class Interpolation
	{
	private:
		Spline *low; //spline for kappa<=1 (remember kappa=K*beta)
		Spline *med; //spline for 1<kappa<=K
		Spline *high; //spline for kappa>K
		const double K; //cardinality of the space
		const double maximum_entropy; //maximum entropy is log(K)
		const double low_kappa_threshold; //lower threshold for kappa (K is the upper threshold)
		const long nodes; //the number of nodes in the splines (roughly)
		const double slope; //slope for the kappa<=1 regime
		double low_xi_threshold; //lower threshold for xi
		double high_xi_threshold; //upper threshold for xi

		const double low_var(const double &kappa) const {return slope*kappa;} //low range transformed variable
		const double low_kappa(const double &var) const {return var/slope;} //low range reverse transform
		const double med_var(const double &kappa) const {return log(kappa) - psi0_1;} //medium range transformed variable
		const double med_kappa(const double &var) const {return exp(var + psi0_1);} //medium range reverse transform
		const double high_var(const double &kappa) const {return maximum_entropy - K/(2*kappa);} //high range transfomed variable
		const double high_kappa(const double &var) const {return 0.5*K/(maximum_entropy-var);} //high range reverse transform
	public:
		/**
		 * Function xi(kappa,beta), see equation 4, is the first moment (a priori
		 * entropy) of the Dirichlet family of priors. This function returns xi
		 * given its argument kappa and data member K.
		 */
		const double xi(const double &kappa) const
		{
			if(gsl_isinf(kappa))
				return maximum_entropy;
			return psi(kappa+1.0) - psi(kappa/K+1.0);
		}

		/**
		 * This function returns the derivative with respect to kappa of member
		 * function xi given its argument kappa and data member K.
		 */
		const double d_xi(const double &kappa) const {return psi1(kappa+1.0) - 1/K*psi1(kappa/K+1.0);}

		/**
		 * This function returns the inverse of member function xi, that is it
		 * returns kappa as a function of its argument xi.
		 */
		const double inv_xi(const double &xi) const
		{
			if((xi<0.0) || (xi>maximum_entropy))
				throw std::range_error("Input xi lies outside the allowed range [0,log(K)].");
			if(xi<low_xi_threshold)
				return low_kappa(low->eval(xi));
			else if(xi<high_xi_threshold)
				return med_kappa(med->eval(xi));
			else if(xi<maximum_entropy)
				return high_kappa(high->eval(xi));
			else
				return GSL_POSINF;
		}

		/**
		 * This is the constructor for the Interpolation object. It handles
		 * creation of the splines between xi and the asymptotic expansion
		 * of kappa.
		 */
		Interpolation(const double &cardinality, const double precision):
			K(cardinality), maximum_entropy(log(K)), low_kappa_threshold(1.0), nodes((long)(1/precision)), slope((K-1)/K*psi1(1.0))
		{
			double *x; //transformed variable
			double *y; //xi
			double range; //range of the transformed variable
			double d_tr; //step size for transformed varible

			low_xi_threshold = xi(low_kappa_threshold); //lower threshold for xi (approaches 1 as K goes to infinity)
			high_xi_threshold = xi(K); //upper threshold for xi (approaches log(K)-psi(2) as K goes to infinity)

			x = new double[nodes];
			y = new double[nodes];

			//fill low range
			range = 1.1*low_var(low_kappa_threshold); //extend slightly beyond threshold
			d_tr = range/nodes;
			for(long i=0; i<nodes; i++)
			{
				x[i] = i*d_tr;
				y[i] = xi(low_kappa(x[i]));
			}

			//create the low range spline
			low = new Spline(y,x,nodes);

			//fill medium range
			range = (1.1*med_var(K) - 0.9*med_var(low_kappa_threshold));
			d_tr = range/nodes;
			x[0] = 0.9*med_var(low_kappa_threshold);
			for(long i=1; i<nodes; i++)
			{
				x[i] = x[i-1]+d_tr;
				y[i] = xi(med_kappa(x[i]));
			}

			//create the medium range spline
			med = new Spline(y,x,nodes);

			//fill high range
			range = maximum_entropy - 0.9*high_var(K);
			d_tr = range/nodes;
			x[0] = 0.9*high_var(K);
			for(long i=1; i<nodes; i++)
			{
				x[i] = x[i-1]+d_tr;
				y[i] = xi(high_kappa(x[i]));
			}

			//create the high range spline
			high = new Spline(y,x,nodes);

			//clean up
			delete[] x;
			delete[] y;
		} //end of Interpolation constructor
		
		~Interpolation()
		{
			delete low;
			delete med;
			delete high;
		}
	}; //end of class Interpolation

	/**
	 * @brief Handles all NSB bookkeeping and calculations.
	 * This is the main class to do (and keep) all NSB calculations.
	 */
	class Computation
	{
	private:
		double N; /**< Total number of words observed (independent samples) */
		const double K; /**< Maximum number of possible words (bins) */
		double K1; /**< Unique number of words observed (non-empty bins) */
		double K2; /**< Number of words observed at least twice (non-empty non-singleton bins) */

		const double *nx; /**< Array of unique, non-zero word counts (unique, non-zero bin occupancy values) */
		const double *kx; /**< Array of counts of nx word counts (bins with nx occupancy) */
		double *nxng1; /**< Array of unique, non-zero, non-singleton word counts (unique, non-zero, non-singleton bin occupancies) */
		double *kxng1; /**< Array of counts of nxng1 word counts (bins with nxng1 occupancy) */

		const int size_nx; /**< Length of nx (and kx) */
		int size_nxng1; /**< Length of nxng1 (and kxng1) */

		const Interpolation *p_splines; /**< Pointer to the interpolation splines for given K and precision_error */

		const double maximum_entropy; /**< Maximum entropy (log(K)) */
		const double precision_error;	/**< Allowed precision error */
		const int find_var; /**< Flag to indicate whether to calculate the variance of the entropy (vSnsb) */

		double Snsb, vSnsb; /**< NSB entropy and its variance (nats) */
		double Scl; /**< Entropy at the saddle point */
		double Sas, dSas; //asymptotic value of entropy and its standard deviation or variance or derivative
		double Bcl, xicl, dxicl; //saddle point values for beta (or kappa) and xi and its standard deviation or variance or derivative
		double mlog; //value of the -log(evidence) at the saddle

		static const int maxcounter = 40; /**< Maximum number of Newton-Raphson iterations */
		NSB_WARN warncode; /**< Error code (0 is no errors) */

		/**
		 * Report and record warnings.
		 */
		void warning(NSB_WARN warn)
		{
			warnings.push_back(p_warn_names[warn]);
			warncode = warn;
		}

		/**
		 * Calculate nx and kx for samples with nx>1, and K1 and K2.
		 */
		void make_ng1() 
		{
			int *p_ng1; //temporary variable to store result of nx>1
			p_ng1 = new int[size_nx];
			for(int i=0; i<size_nx; i++)
				p_ng1[i] = (nx[i]>1.05); //find nx>1 (use 1.05 to make sure that reals don't cause problem)

			for(int i=0; i<size_nx; i++)
				size_nxng1 += p_ng1[i]; //number of elements in nx greater than 1

			if(size_nxng1) //there exist words spoken more than once
			{
				nxng1 = new double[size_nxng1];
				kxng1 = new double[size_nxng1];

				int j = 0;
				for(int i=0; i<size_nx; i++)
					if(p_ng1[i])
					{
		  				nxng1[j] = nx[i];
						kxng1[j] = kx[i];
						j++;
					}
			}
			delete[] p_ng1;

			for(int i=0; i<size_nx; i++)
				K1 += kx[i];
			for(int i=0; i<size_nxng1; i++)
				K2 += kxng1[i];
		}

		/**
		 * This function finds the position of the minimum of the a posteriori
		 * evidence and the variance around it. The integration variable is xi
		 * (a priori entropy), and values for Bcl (classical value of B), xicl
		 * (classical value of xi), and dxicl (standard deviation near the
		 * classical value) are calculated.
		 */
		const NSB_WARN max_evidence()
		{
			status.push_back("Finding the saddle point.");

			if(round(K1)==round(N)) //no coincidences
			{
				Bcl = xicl = GSL_POSINF;
				warning(NSB_NOCOINC);
				dxicl = GSL_NAN;
			}
			else if(round(K1)==1) //all data coincides
			{
				Bcl = xicl = 0.0;
				warning(NSB_ALLCOINC);
				dxicl = GSL_NAN;
			}
			else //some non-trivial value of B0 and Bcl have to be calculated
			{
				double B0 = 0.0;
				//summing the series
				{ 
					const int order = 10; //calculate B to this order in epsilon=(N-K1)/N
					const double N2 = gsl_pow_2(N);
					const double N3 = gsl_pow_3(N);
					const double N4 = gsl_pow_4(N);
					const double N5 = gsl_pow_5(N);
					const double N6 = gsl_pow_6(N);
					const double N7 = gsl_pow_7(N);
					const double N8 = gsl_pow_8(N);
					const double N9 = gsl_pow_9(N);
					const double N10 = gsl_pow_int(N,10);
					const double N11 = gsl_pow_int(N,11);

					const double ovrN = 1/N;
					const double Nm = N-1.0;
					const double Nm2 = gsl_pow_2(Nm);
					const double Nm3 = gsl_pow_3(Nm);
					const double Nm4 = gsl_pow_4(Nm);
					const double Nm5 = gsl_pow_5(Nm);
					const double Nm6 = gsl_pow_6(Nm);
					const double Nm7 = gsl_pow_7(Nm);
					const double Nm8 = gsl_pow_8(Nm);
					const double Nm9 = gsl_pow_9(Nm);
					const double Nm10 = gsl_pow_int(Nm,10);

					//coefficients of the expansion of B in powers of epsilon (calculated by Mathematica)
					const double b[] = {
						(-1.0 + N)/(2.0*N), //b(-1) (see equation 19)
						(-2.0 + ovrN)/3.0, //b(0) (see equation 20)
						(2.0 + N - N2)/(9.0*N - 9.0*N2), //b(1) (see equation 21)
						(2.0*(2.0 - 3.0*N - 3.0*N2 + 2.0*N3))/(135.0*Nm2*N), //b(2)
						(4.0*(22.0 + 13.0*N - 12*N2 - 2.0*N3 + N4))/(405.0*Nm3*N), //b(3)
						(4.0*(-40.0 + 58.0*N + 65.0*N2 - 40.0*N3 - 5.0*N4 + 2.0*N5))/(1701.0*Nm4*N), //b(4)
						(4.0*(-9496.0 - 6912.0*N + 5772.0*N2 + 2251.0*N3 - 1053.0*N4 - 87.0*N5 + 29.0*N6))/(42525.0*Nm5*N), //b(5)
						(16.0*(764.0 - 1030.0*N - 1434.0*N2 + 757.0*N3 + 295.0*N4 - 111.0*N5 - 7.0*N6 + 2.0*N7))/(18225.0*Nm6*N), //b(6)
						(16.0*(167000.0 + 142516.0*N - 108124.0*N2 - 66284.0*N3 + 26921.0*N4 + 7384.0*N5 - 2326.0*N6 - 116.0*N7 + 29.0*N8))/(382725*Nm7*N), //b(7)
						(16.0*(-17886224.0 + 22513608.0*N + 37376676.0*N2 - 17041380.0*N3 - 11384883.0*N4 + 3698262.0*N5 + 846930.0*N6 - 229464.0*N7 - 9387.0*N8 + 2086*N9))/(37889775.0*Nm8*N), //b(8)
						(16.0*(-4166651072.0 - 3997913072.0*N + 2783482560.0*N2 + 2290151964.0*N3 - 803439834.0*N4 - 395614251.0*N5 + 108055443.0*N6 + 20215218.0*N7 - 4805712.0*N8 - 165395.0*N9 + 33079.0*N10))/(795685275*Nm9*N), //b(9)
						(32.0*(52543486208.0 - 62328059360.0*N - 118489458160.0*N2 + 47185442088.0*N3 + 44875379190.0*N4 - 12359832987.0*N5 - 5400540075.0*N6 + 1272974916.0*N7 + 200644800.0*N8 - 42495955.0*N9 - 1255067.0*N10 + 228194.0*N11))/(14105329875.0*Nm10*N) //b(10)
					};

					//calculate the value of B0 as a series exansion in powers of epsilon=(N-K1)/N (see equation 18)
					B0 = N*series(order+2, b, (N-K1)/N, -1);
				} //clear all B0 expansion related variables

				if(B0<0) //bad but can still recover precision with Newton-Raphson root finding
				{
					B0 = precision_error; //assign precision error to B0
					warning(NSB_SER_B0);
				}

#ifdef DEBUG
				std::cout << "Value of B0 as a series expansion: " << B0 << "\n";
#endif

				//use Newton-Raphson to polish the value of B0 and find B0 to the desired precision
				{
					double dB0, F, dF = 0.0;
					int counter = 0;
					do{
						counter++;
						F = K1/B0 + psi(B0) - psi(B0+N); //see equation 15
						dF = -K1/gsl_pow_2(B0) + psi1(B0) - psi1(B0+N); //derivative of F
						dB0 = -F/dF;
						B0 += dB0;
						if(B0<=0.0)
						{
							B0 = precision_error;
							dB0 = 0.0;
							warning(NSB_NR_B0_SIGN);
						}
					} while ((fabs(dB0)>fabs(B0*precision_error)) && (counter<=maxcounter));

					if(counter==maxcounter)
						warning(NSB_NR_B0);
				} //clear all Newton-Raphson related variables

#ifdef DEBUG
				std::cout << "Starting value of Bcl (also B0 after Newton-Raphson polish): " << B0 << "\n";
#endif

				Bcl = B0; //take the calculated value as the first approximation for Bcl
				{
					const int order_K = 4; //number of terms in the series (orders up to 10 are in the nsb-entropy-oct distribution in the file other_orders_K.m)

					//temporary variables
					const double pg1B0 = psi1(B0);
					const double pg1NB0 = psi1(N+B0);
					const double denum = K1/gsl_pow_2(B0) - pg1B0 + pg1NB0; //short for denumerator which probably means denominator
					const double pg2B0 = psi2(B0);
					const double pg2NB0 = psi2(N+B0);
					const double pg21 = psi2(1);
					const double pg3B0 = psi3(B0);
					const double pg3NB0 = psi3(N+B0);
					const double pg4B0 = psi4(B0);
					const double pg4NB0 = psi4(N+B0);

					const double f0 = dot_f_a_b(size_nxng1, nxng1, kxng1, &psi);
					const double d1f0 = dot_f_a_b(size_nxng1, nxng1, kxng1, &psi1);
					const double d2f0 = dot_f_a_b(size_nxng1, nxng1, kxng1, &psi2);
					const double d3f0 = dot_f_a_b(size_nxng1, nxng1, kxng1, &psi3);

					const double B02 = gsl_pow_2(B0);
					const double B03 = gsl_pow_3(B0);
					const double B04 = gsl_pow_4(B0);
					const double B05 = gsl_pow_5(B0);

					const double PI_2 = gsl_pow_2(M_PI);
					const double PI_4 = gsl_pow_2(PI_2);

					//all 4 expansion orders
					double b[order_K];
					b[0] = B02*(M_EULER*K2 + f0)/(B02*denum);
					const double b02 = gsl_pow_2(b[0]);
					const double b03 = gsl_pow_3(b[0]);
					const double b04 = gsl_pow_4(b[0]);
					b[1] = (K2*PI_2*B0 - (6.0*K1*b02)/B03 - 3*b02*pg2B0 + 3.0*b02*pg2NB0 - 6.0*B0*d1f0)/(-6.0*denum);
					const double b12 = gsl_pow_2(b[1]);
					b[2] = (K2*PI_2*b[0] + (6.0*K1*b03)/B04 -(12*K1*b[0]*b[1])/B03 + 3.0*K2*B02*pg21 - 6.0*b[0]*b[1]*pg2B0 + 6.0*b[0]*b[1]*pg2NB0 - b03*pg3B0 + b03*pg3NB0 - 6.0*b[0]*d1f0 - 3.0*B02*d2f0)/(-6.0*denum);
					b[3] = -(-(K2*PI_4*B03)/90.0 + (K1*b04)/B05 - (K2*PI_2*b[1])/6.0 - (3.0*K1*b02*b[1])/B04 + (K1*b12)/B03 + (2.0*K1*b[0]*b[2])/B03 - K2*B0*b[0]*pg21 + ((b12 + 2.0*b[0]*b[2])*pg2B0)/2.0 - ((b12 + 2*b[0]*b[2])*pg2NB0)/2.0 + (b02*b[1]*pg3B0)/2.0 - (b02*b[1]*pg3NB0)/2.0 + (b04*pg4B0)/ 24.0 - (b04*pg4NB0)/24.0 +  b[1]*d1f0 + B0*b[0]*d2f0 + (B03*d3f0)/6.0)/(-denum);

					Bcl += series(order_K, b, 1/K, 1); //get the expansion
				} //clear all Bcl expansion related variables

#ifdef DEBUG
				std::cout << "Value of Bcl after adding expansion terms: " << Bcl << "\n";
#endif

				const double Ksq = gsl_pow_2(K);
				//use Newton-Raphson to polish the value of Bcl and find Bcl to the desired precision
				{
					double dBcl, F, dF = 0.0;
					int counter = 0;
					do{
						counter++;
						F = 1/K*dot_f_aplusa0_b(size_nxng1, nxng1, kxng1, Bcl/K, &psi) - K2/K*psi(1.0+Bcl/K) + K1/Bcl + psi(Bcl) - psi(Bcl+N);
						dF = 1/Ksq*dot_f_aplusa0_b(size_nxng1, nxng1, kxng1, Bcl/K, &psi1) - K2/Ksq*psi1(1+Bcl/K) - K1/gsl_pow_2(Bcl) + psi1(Bcl) - psi1(Bcl +N);
						dBcl = -F/dF;
						Bcl += dBcl;
						if(Bcl<=0.0)
						{
							dBcl = 0.0;
							Bcl = precision_error;
							warning(NSB_NR_BCL_SIGN);
						}
					} while ((fabs(dBcl)>fabs(Bcl*precision_error)) && (counter<=maxcounter));

					if(counter==maxcounter)
					{
						warning(NSB_NR_BCL);
						errors.push_back(p_warn_names[NSB_NR_BCL]);
					}
					else if((warncode==NSB_NR_B0) || (warncode==NSB_NR_B0_SIGN) || (warncode==NSB_SER_B0))
						warning(NSB_OK); //method seems to have recovered
				} //clear all Newton-Raphson related variables

#ifdef DEBUG
				std::cout << "Value of Bcl after Newton-Raphson polishing: " << Bcl << "\n";
#endif

				//calculate xicl and dxicl
				const double dBcl = 1/Ksq*dot_f_aplusa0_b(size_nxng1, nxng1, kxng1, Bcl/K, &psi1) - K2/Ksq*psi1(1.0+Bcl/K) - K1/gsl_pow_2(Bcl) + psi1(Bcl) - psi1(Bcl +N);
				xicl = xi(Bcl);
				dxicl = 1/sqrt(-dBcl/gsl_pow_2(d_xi(Bcl)));
			} //end of else (non-trivial value for B0)

			return warncode;
		} //end of member function max_evidence

		/**
		 * Mean value of the a posterior entropy for given value of B (it is unclear
		 * here whether B refers to beta or kappa).
		 */
		const double meanS(const double &B) const
		{
			if(gsl_isinf(B)) //infinite beta yields maximum entropy
				return maximum_entropy;
			const double ovrNB = 1/(N+B);
			const double BoK = B/K;
			const double prod = dot_f_aplusa0_bplusb0_c(size_nx, nx, nx, kx, BoK+1.0, BoK, &psi);

			//sum is written in a form to avoid the loss of precision as K approaches Inf
			return psi(N+B+1.0) - ovrNB*prod - B*ovrNB*(1.0-K1/K)*psi(BoK+1.0);
		}

		/**
		 * Mean value of the a posterior entropy squared for given value of B (it is
		 * unclear here whether B refers to beta or kappa).
		 */
		const double meanS2(const double &B) const
		{
			if(gsl_isinf(B)) //infinite beta yields maximum entropy squared
				return gsl_pow_2(maximum_entropy);

			double f = 0.0; //the variable to return later

			//temporary variables
			const double BoK = B/K;
			const double p0NB2 = psi(N+B+2.0);
			const double p1NB2 = psi1(N+B+2.0);
			const double pb1 = psi(BoK + 1.0) - p0NB2;

			//temporary arrays
			double *pnxb1;
			double *nxb;
			pnxb1 = new double[size_nx];
			nxb = new double[size_nx];

			for(int i=0; i<size_nx; i++)
			{
				nxb[i] = nx[i] + BoK;
				pnxb1[i] = psi(nxb[i] + 1.0) - p0NB2;
			}

			//-----------------------------------------------------
			//sum over all i and j (include i==j terms then correct for it)

			//ni*nj~=0 contribution
			for(int i=0; i<size_nx; i++)
				for(int j=0; j<size_nx; j++)
					f += nxb[i]*pnxb1[i]*kx[i]*nxb[j]*pnxb1[j]*kx[j] - nxb[i]*kx[i]*nxb[j]*kx[j]*p1NB2;

			//ni*b contribution
			for(int i=0; i<size_nx; i++)
				f += 2.0*B*(1-K1/K)*nxb[i]*(pnxb1[i]*pb1 - p1NB2)*kx[i];

			//b*b contribution 
			f += (1-K1/K)*(1-(K1+1.0)/K)*B*B*(pb1*pb1-p1NB2);

			//correct for overcounting
			for(int i=0; i<size_nx; i++)
				f -= gsl_pow_2(nxb[i]*pnxb1[i])*kx[i] - nxb[i]*nxb[i]*kx[i]*p1NB2;

			//-----------------------------------------------------
			//i term

			//ni contribution
			for(int i=0; i<size_nx; i++)
				f += (nxb[i]*(nxb[i]+1.0) * (gsl_pow_2(psi(nxb[i]+2.0) - p0NB2) + psi1(nxb[i]+2.0) - p1NB2))*kx[i];

			//b contribution
			f += B*(1-K1/K)*(1.0+BoK) * (gsl_pow_2(psi(2.0+BoK)-p0NB2) + psi1(BoK+2.0) - p1NB2);

			//-----------------------------------------------------
			//normalize
			f /= ((N+B)*(N+B+1.0));

			delete[] pnxb1;
			delete[] nxb;

			return f;
		}

		/**
		 * Computes the "action" (the negative logarithm of the evidence) for
		 * the integral over xi (see NSB method for calculating entropies of
		 * discrete pdfs). Does not include the contribution from the prior over
		 * xi. Note that even though the integration variable is xi, the argument
		 * of this function is B (it is unclear here whether B refers to beta or
		 * kappa).
		 */
		const double mlog_evidence(const double &B) const
		{
			double f = 0.0;

			if((B<=0.0) || gsl_isinf(B)) //infinite evidence for negative or infinite B
				return GSL_POSINF;

			const double BoK = B/K;
			if(size_nxng1) //if there are coincidences
				f += -dot_f_aplusa0_b(size_nxng1, nxng1, kxng1, BoK, &gsl_sf_lngamma) + K2*gsl_sf_lngamma(1+BoK);

			//calcuate f += -K1*log(B) + gammaln(B+N) - gammaln(B) in Matlab/Octave
			//notation but to avoid loss in precision treat different regimes of N
			//and B differently - for aymptotically large B and B/N (B>100 and
			//N<0.01*B) expand gammaln-gammaln = psi*N + psi_1/2*N^2 + ...
			const int large = (B>GSL_MAX(100.0, 100.0*N));

			//polygamma(n,B)~B^(-n) thus we expand in N/B - unclear which of
			//expansion parameters is the worst or how many series terms will be
			//needed but for N=K1 (worst case) leading term in f is
			//psi_asymp(B)/N~N/B and we need 10^(-15) precision relative to that
			//term with series expansion of the form f = leading + psi_1/2!*N^2 +
			//psi_2/3!*N^3 + ... = leading + (N/B + (N/B)^2 + ...)
			if(large)
			{
				const int nterms  = (int)ceil(fabs((-15.0 -log10(N))/log10(N/B))) + 1;
				double *ifac;
				ifac = new double[nterms];
				factorial_fraction(nterms, 2.0, ifac); //populating ifac with 1/factorial
				f += series_f(nterms, ifac, N, 2, &gsl_sf_psi_n, B, 1);
				delete[] ifac;
				f += (N-K1)*log(B) + psi_asymp(B)*N;
			}
			else //no asymptotic expansion needed
				f += -K1*log(B) + gsl_sf_lngamma(B+N) - gsl_sf_lngamma(B);

			return f;
		}

		/**
		 * This is the main computational routine for the NSB method; it is
		 * functionally equivalent to find_nsb_entropy.m in the Octave code. New
		 * values for Snsb (NSB entropy estimate), vSnsb (variance of the estimate),
		 * Scl (entropy at the saddle point), dScl (standard deviation at the saddle
		 * point), xicl and dxicl (saddle point entropy and stardard deviation) are
		 * calculated. A warning code of type NSB_WARN is returned to indicate the
		 * success or failure of the algorithm to find a solution.
		 */
		const NSB_WARN calculate()
		{
			const int maxIntegrands = 3;
			int nIntegrands; //number of integrands
			if(find_var)
				nIntegrands = 3;
			else
				nIntegrands = 2;

			//allocate three different integrands
			gsl_function integrands[maxIntegrands];
			integrands[0].function = &g_int_1;
			integrands[0].params = this;
			integrands[1].function = &g_int_S;
			integrands[1].params = this;
			if(find_var)
			{
				integrands[2].function = &g_int_S2;
				integrands[2].params = this;
			}

			const char *integrand_names[3] = {"normalization","S", "S^2"};
			double integrals[maxIntegrands]; //values of integrals (result in gsl_integration_qag)
			double abserr[maxIntegrands]; //estimate of absolute errors (abserr in gsl_integration_qag)

			//prepare for saddle point integration
			make_ng1(); //fill in n>1 arrays (now have nxng1, kxng1, K1, and K2)
			max_evidence(); //find the saddle (now have Bcl, xicl, and dxicl)
			Scl = meanS(Bcl);

			status.push_back("Prepared to do saddle point integration.");

			//set limits of integration
			double edges[] = {precision_error, maximum_entropy-precision_error}; //start with edges precision_error away from bounds
			double xilim[2]; //actual integration boundaries
			double delta = 0.0; //delta is the interval around the peak on which gaussian approximation falls to "precision_error/2" on each side

			if(warncode!=NSB_OK)
			{
				//some problem with finding saddle value - switch to full range
				warnings.push_back("Switching to integration over the full range.");

				//integrating over the whole range
				xilim[0] = edges[0];
				xilim[1] = edges[1];
				mlog = mlog_evidence(inv_xi(xilim[0])); //don't know the value at the saddle

				//now recurse through the entire domain to find smaller mlog
				for(double x=xilim[0]; x<=xilim[1]; x+=(xilim[1]-xilim[0])/100)
					mlog = GSL_MIN(mlog,mlog_evidence(inv_xi(x)));
			}
			else
			{
				std::ostringstream stream;
				stream << "Expect S=" << NAT2BIT(Scl) << " +/-" << NAT2BIT(dxicl) << ". Integrating around the peak.";
				status.push_back(stream.str());
				mlog = mlog_evidence(Bcl); //value at the saddle
				delta = dierfc(precision_error/2.0)*M_SQRT2;
			}

			Workspace ws; //create integration workspace
			for(int i=0; i<nIntegrands; i++)
			{
				if(delta>0.0) //if integrating around the peak need the limits of integration
				{
					//the value of the integrand at xicl
					const double cent = (*(integrands[i].function))(xicl,integrands[i].params);

					//if the integral was purely Gaussian the value at +/-delta would have the
					//required precision 1/sqrt(2pi)*exp(-(delta^2)/2) so as a safety we
					//require the value at the limits of integration to be ten times smaller
					const double good = 0.1*cent*exp(-(delta*delta)/2)*1/sqrt(2*M_PI);

					//increase the integration window until the integrand is small at the
					//edges by starting with xicl +/-delta*dxicl and expanding by +/-0.5*dxicl
					//while making sure that the limits don't go outside the edges
					xilim[0] = GSL_MAX(xicl-delta*dxicl, edges[0]);
					xilim[1] = GSL_MIN(xicl+delta*dxicl, edges[1]);

					double limval[2] = {(*(integrands[i].function))(xilim[0],integrands[i].params), (*(integrands[i].function))(xilim[1],integrands[i].params)};

					for(int j=0; j<2; j++) //do a loop over lower and upper limits
					{
						double window = 0.0; //limits=(delta+window)*dxicl
						while((limval[j]>good) && (xilim[j]!=edges[j]))
						{
							if(window>10.0) //fast growth for large windows
								window *= 1.2;
							else //regular growth
								window += 0.5;

							//increase the corresponding limit but make sure we are still in range
							if(!j) //for left limit
								xilim[0] = GSL_MAX(edges[0], xicl - (delta+window)*dxicl);
							else //for right limit
								xilim[1] = GSL_MIN(edges[1], xicl + (delta+window)*dxicl);
							limval[j] = (*(integrands[i].function))(xilim[j],integrands[i].params);
						}
					}
				}

				if(xilim[0]>=xilim[1]) //a crude hack to correct one possible problem
				{
					xilim[0] = edges[0];
					xilim[1] = edges[1];
				}

				std::ostringstream stream;
				stream << "Doing " << integrand_names[i] << " integral within limits: " << NAT2BIT(xilim[0]) << "<xi<" << NAT2BIT(xilim[1]) << ".";
				status.push_back(stream.str());

				int ec = gsl_integration_qag(&integrands[i], xilim[0], xilim[1], 0.0, precision_error, ws.get_size(), GSL_INTEG_GAUSS21,  ws.get_workspace(), &integrals[i], &abserr[i]);
				if(abserr[i]/integrals[i]>precision_error)
					errors.push_back(p_warn_names[NSB_INT_NOCONV]);

				if(ec!=GSL_SUCCESS) //handle error
					throw std::runtime_error("GSL integration failed.");
			} //end of the for loop over int_1, int_S, and int_S2

			Snsb = integrals[1]/integrals[0];
			if(find_var)
				vSnsb = integrals[2]/integrals[0] - gsl_pow_2(Snsb);

			std::ostringstream stream;
			if(find_var)
				stream << "Found S=" << NAT2BIT(Snsb) << " +/-" << NAT2BIT(sqrt(vSnsb)) << ".";
			else
				stream << "Found S=" << NAT2BIT(Snsb) << ".";
			status.push_back(stream.str());

			const double D = N-K1; //number of coincidences
			if(D)
			{
				Sas = -psi(1.0) - log(2.0) + 2.0*log(N) - psi(D);
				dSas = sqrt(psi1(D));
			}
			else
				Sas = dSas = GSL_POSINF;

			return warncode;
		} //end of member function calculate

		//a priori entropy xi and its derivative d_xi for given value of kappa=K*beta and its inverse inv_xi
		const double xi(const double &kappa) const {return p_splines->xi(kappa);}
		const double d_xi(const double &kappa) const {return p_splines->d_xi(kappa);}
		const double inv_xi(const double &xi) const {return p_splines->inv_xi(xi);}
	public:
		std::vector<std::string> status; /**< Status messages */
		std::vector<std::string> warnings; /**< Warning messages */
		std::vector<std::string> errors; /**< Error messages */
		const double get_Snsb() const {return Snsb;} //nsb entropy
		const double get_vSnsb() const {return vSnsb;} //nsb posterior variance
		const double get_Scl() const {return Scl;} //value of S at the saddle
		const double get_Sas() const {return Sas;} //small coincidence asymptotic
		const double get_dSas() const {return dSas;} //small coincidence asymptotic standard deviation
		const double get_Bcl() const {return Bcl;} //saddle value for beta
		const double get_xicl() const {return xicl;} //saddle value for a priori entropy
		const double get_dxicl() const {return dxicl;} //standard deviation around the saddle
		const NSB_WARN get_warncode() const {return warncode;}; //warning code

		/**
		 * Integrand for the normalization integral. The only parameter here is
		 * the value of the argument. The rest is implicit through class structures.
		 */
		const double int_1(const double &xi) const
		{
			return exp(-mlog_evidence(inv_xi(xi)) + mlog)/maximum_entropy;
		}

		/**
		 * Integrand for the S integral. The only parameter here is the value
		 * of the argument. The rest is implicit through class structures.
		 */
		const double int_S(const double &xi) const
		{
			const double kappa = inv_xi(xi);
			return exp(-mlog_evidence(kappa) + mlog)/maximum_entropy*meanS(kappa);
		}

		/**
		 * Integrand for the S^2 integral. The only parameter here is the value
		 * of the argument. The rest is implicit through class structures.
		 */
		const double int_S2(const double &xi) const
		{
			const double kappa = inv_xi(xi);
			return exp(-mlog_evidence(kappa) + mlog)/maximum_entropy*meanS2(kappa);
		}

		/**
		 * @brief Initialize a Computation object.
		 * This is the constructor for the Computation object. Its inputs reflect the
		 * format used by find_nsb_entropy.m in the package nsb-entropy-oct_1.11.tgz
		 * at http://nsb-entropy.sourceforge.net/, as this is more directly
		 * compatible with the internal toolkit data structures.
		 */
		Computation(const double *p_nx, const double *p_kx, const int size, const int samples, const double cardinality, const double precision, const int var_flag):
			N(samples), K(cardinality), K1(0.0), K2(0.0), size_nx(size), size_nxng1(0), nx(p_nx), kx(p_kx), nxng1(NULL), kxng1(NULL), maximum_entropy(log(K)), precision_error(precision), find_var(var_flag), p_splines(NULL), Snsb(0.0), vSnsb(0.0), Scl(0.0), Sas(0.0), dSas(0.0), Bcl(0.0), xicl(0.0), dxicl(0.0), mlog(0.0), warncode(NSB_OK)
		{
			//initialize messages
			status.reserve(10);
			warnings.reserve(3);
			errors.reserve(1);

			//return immediately if data is trivial
			if(K<2)
			{
				warnings.push_back("Cardinality is 1; results are trivial.");
				return;
			}

			//prepare error handler
			gsl_error_handler_t *gsl_error_handler;
			gsl_error_handler = gsl_set_error_handler_off();

			//get splines
			try
			{
				std::ostringstream stream;
				stream << "Creating spline data for K=" << K << ".";
				status.push_back(stream.str());
				p_splines = new Interpolation(K,precision_error);
			}
			catch(std::bad_alloc)
			{
				errors.push_back("Failed to properly instantiate interpolation splines.");
				return;
			}

			//calculate results
			calculate();

			//reset error handler
			gsl_set_error_handler(gsl_error_handler);
		} 

		/**
		 * The destructor for the Computation object deallocates all arrays.
		 */
		~Computation()
		{
			delete[] nx;
			delete[] kx;
			delete[] nxng1;
			delete[] kxng1;
			delete p_splines;
		}
	}; //end of class Computation

	/**
	 * Integration routine envelopes used to call the member integration routines
	 * from Computation class objects (g_ in the name stands for "global").
	 */
	double g_int_1(double xi, void *obj)
	{
		Computation *p = (Computation*)obj;
		return p->int_1(xi);
	}
	double g_int_S(double xi, void *obj)
	{
		Computation *p = (Computation*)obj;
		return p->int_S(xi);
	}
	double g_int_S2(double xi, void *obj)
	{
		Computation *p = (Computation*)obj;
		return p->int_S2(xi);
	}
} //end of namespace nsb-entropy

/**
 * @brief Estimate NSB entropy from a vector of word counts.
 * Given a hist1d vector of word counts, and an options structure, this
 * function returns an estimate of the entropy based on the NSB method. Note:
 * this function conforms to the toolkit standard for inputs and outputs.
 */
int entropy_nsb(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy)
{
	double cardinality; /**< Number of possible words (or alphabet size in NSB lingo) */
	int n_unique_wordcnt; /**< Number of unique values in in->wordcnt */
	double *p_nx; /**< Counts array */
	int *p_nx_to_wordcnt; /**< Index from p_nx to in->wordcnt */
	int *p_wordcnt_to_nx; /**< Index from in->wordcnt to p_nx */
	int *p_kx; /**< Counts of counts array */
	double *p_kx_double; /**< Counts of counts array (as doubles) */

	//get possible words from options
	if(opts->possible_words>0)
		cardinality = opts->possible_words;
	else
		switch((int)(opts->possible_words))
		{
			case 0: //Inf
				cardinality = INFINITY;
				break;
			case -1: //recommended
				cardinality = (double)in->C; //based on README in nsb-entropy-oct code
				break;
			case -2: //unique
				cardinality = (double)in->C;
				break;
			case -3: //total
				cardinality = (double)in->P;
				break;
			case -4: //possible
				cardinality = max_possible_words(in,0);
				break;
			case -5: //min_tot_pos
				cardinality = max_possible_words(in,1);
				break;
			case -6: //min_lim_tot_pos
				cardinality = MIN(1e5,max_possible_words(in,1));
				break;
		}

	//allocate memory for NSB parameters
	p_nx = (double *)malloc((in->C)*sizeof(double));
	p_nx_to_wordcnt = (int *)malloc((in->C)*sizeof(int));
	p_wordcnt_to_nx = (int *)malloc((in->C)*sizeof(int));
	p_kx = (int *)calloc(in->C,sizeof(int));

	//convert in->wordcnt to counts
	n_unique_wordcnt = UniqueDouble(in->C,in->wordcnt,p_nx,p_nx_to_wordcnt,p_wordcnt_to_nx,p_kx);

	//remove empty bin counts
	if((int)p_nx[0]==0)
	{
		for(int i=1; i<n_unique_wordcnt; i++)
		{
			p_nx[i-1] = p_nx[i];
			p_kx[i-1] = p_kx[i];
		}
		n_unique_wordcnt--;
	}

#ifdef DEBUG
	std::cout << "nx = [ ";
	for(int i=0; i<n_unique_wordcnt; i++)
		std::cout << (int)p_nx[i] << " ";
	std::cout << "]\n";
	std::cout << "kx = [ ";
	for(int i=0; i<n_unique_wordcnt; i++)
		std::cout << p_kx[i] << " ";
	std::cout << "]\n";
#endif

	//reallocate and cast memory
	p_nx = (double *)realloc(p_nx,n_unique_wordcnt*sizeof(double));
	p_kx_double = (double *)malloc(n_unique_wordcnt*sizeof(double));
	for(int i=0; i<n_unique_wordcnt; i++)
		p_kx_double[i] = (double)p_kx[i];

	//determine whether variance is requested (KLUDGE)
	int var_index, find_var = 0;
	if(!strcmp(entropy->name,"nsb"))
		for(int v=0; v<entropy->V; v++)
			if(!strcmp(entropy->ve[v].name,"nsb_var"))
			{
				var_index = v;
				find_var = 1;
				break;
			}

	//create nsb object
	nsb_entropy::Computation nsb(p_nx,p_kx_double,n_unique_wordcnt,in->P,cardinality,opts->nsb_precision,find_var);

	//copy entropy (and variance)
	entropy->value = NAT2BIT(nsb.get_Snsb());
	if(find_var)
		entropy->ve[var_index].value = NAT2BIT(NAT2BIT(nsb.get_vSnsb()));

	//copy status, warning, and error messages
	if(nsb.status.size())
	{
		entropy->messages->i = nsb.status.size();
		entropy->messages->status = (char **)malloc(entropy->messages->i*sizeof(char *));
		for(int i=0; i<entropy->messages->i; i++)
		{
			entropy->messages->status[i] = (char *)malloc((nsb.status[i].size()+1)*sizeof(char));
			strcpy(entropy->messages->status[i],nsb.status[i].c_str());
		}
	}
	if(nsb.warnings.size())
	{
		entropy->messages->j = nsb.warnings.size();
		entropy->messages->warnings = (char **)malloc(entropy->messages->j*sizeof(char *));
		for(int j=0; j<entropy->messages->j; j++)
		{
			entropy->messages->warnings[j] = (char *)malloc((nsb.warnings[j].size()+1)*sizeof(char));
			strcpy(entropy->messages->warnings[j],nsb.warnings[j].c_str());
		}
	}
	if(nsb.errors.size())
	{
		entropy->messages->k = nsb.errors.size();
		entropy->messages->errors = (char **)malloc(entropy->messages->k*sizeof(char *));
		for(int k=0; k<entropy->messages->k; k++)
		{
			entropy->messages->errors[k] = (char *)malloc((nsb.errors[k].size()+1)*sizeof(char));
			strcpy(entropy->messages->errors[k],nsb.errors[k].c_str());
		}
	}

	//copy extra computed values
	const char extras_names[][MAXCHARS] = {"Sas", "dSas", "Scl", "Bcl", "xicl", "dxicl"};
	const double extras_values[] = {NAT2BIT(nsb.get_Sas()), NAT2BIT(nsb.get_dSas()), NAT2BIT(nsb.get_Scl()), NAT2BIT(nsb.get_Bcl()), NAT2BIT(nsb.get_xicl()), NAT2BIT(nsb.get_dxicl())};
	entropy->E = sizeof(extras_names)/MAXCHARS;
	entropy->extras = (struct nv_pair *)malloc(entropy->E*sizeof(struct nv_pair));
	for(int e=0; e<entropy->E; e++)
	{
		strcpy(entropy->extras[e].name,extras_names[e]);
		entropy->extras[e].value = extras_values[e];
	}

	//free allocated memory (note that p_nx and p_kx_double are freed by the nsb object)
	free(p_nx_to_wordcnt);
	free(p_wordcnt_to_nx);
	free(p_kx);

	if(nsb.errors.size())
		return EXIT_FAILURE;
	else
		return EXIT_SUCCESS;
}

/**
 * @brief Estimate the variance in NSB entropy from a vector of word counts.
 * Given a hist1d vector of word counts, and an options structure, this
 * function returns an estimate of the variance in entropy based on the NSB
 * method. Note: this function conforms to the toolkit standard for inputs
 * and outputs, though the actual variance estimate is supplied by the
 * entropy_nsb method.
 */
int variance_nsb(struct hist1d *in,struct options_entropy *opts,struct nv_pair *variance)
{
	return EXIT_SUCCESS;
}

