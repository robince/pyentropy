/*
 *  Copyright 2010, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_sf.h>
#include "toolkit_c.h"

#define BUB_MAXWORDS_LIMIT 100000

struct bub_options{
	double lambda_N;
	int use_weighting;
	double mesh_min;
	double mesh_max;
	double mesh_inc;
	int opt_use_weighting;
	double opt_mesh_min;
	double opt_mesh_max;
	double opt_mesh_inc;
	int bag_threshold;
};

int bag1(int N,double m,double *a,double lambda_0,struct bub_options *bubopts);
int bub1(int N,double m,double *a,double lambda_0,int K,int compat,struct bub_options *bubopts);
void make_binom(int N,int meshsize,double *p,double **B_jN);
void make_weight(double *p,int meshsize,double m,int use_weighting,double *f,double *c);
double compute_error(double *a,int N,double m,int compat,struct bub_options *bubopts);

int entropy_bub(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy)
{
	double *a;
	int *h;
	double H;
	struct bub_options *bubopts;
	double eps,epsm;
	int N,i,j,max_bin;
	double m;
	int status;

	N = (*in).P;

	if(opts->possible_words_flag==0)
	{
		if(opts->bub_possible_words_strategy==0)
			m = (double)in->C; 
		else if(opts->bub_possible_words_strategy==1)
			m = (double)in->N; 
		else if(opts->bub_possible_words_strategy==2)
			m = MIN(max_possible_words(in,1),BUB_MAXWORDS_LIMIT);
		else
			m = (double)in->C;
	}
	else
		if(opts->possible_words>0)
			m = opts->possible_words;
		else
			switch((int)(opts->possible_words))
			{
				case 0: /* Inf */
					entropy->value = NAN; /* this value based on piece-meal derivation when m=INFINITY */
					return EXIT_SUCCESS;
				case -1: /* recommended */
					m = (double)in->C;
					break;
				case -2: /* unique */
					m = (double)in->C;
					break;
				case -3: /* total */
					m = (double)in->P;
					break;
				case -4: /* possible */
					m = max_possible_words(in,0);
					break;
				case -5: /* min_tot_pos */
					m = max_possible_words(in,1);
					break;
				case -6: /* min_lim_tot_pos */
					m = MIN(BUB_MAXWORDS_LIMIT,max_possible_words(in,1));
					break;
			}

	/* First find a vector */
	a = (double *)malloc((N+1)*sizeof(double));

	bubopts = (struct bub_options *)malloc(sizeof(struct bub_options));
	(*bubopts).bag_threshold = 20;
	if(N<(*bubopts).bag_threshold)
	{
		(*bubopts).lambda_N = 0;
		(*bubopts).mesh_min = 0;
		(*bubopts).mesh_max = 1;
		(*bubopts).mesh_inc = 1/(5*(double)N);
		(*bubopts).use_weighting = 0;
		status = bag1(N,m,a,(*opts).bub_lambda_0,bubopts);
	}
	else
	{
		eps = (1/(double)N)*1e-10;
		epsm = (1/m)*1e-10;
		(*bubopts).lambda_N = 1;
		(*bubopts).mesh_min = eps;
		(*bubopts).mesh_max = MIN(1,3/(double)N)-eps;
		(*bubopts).mesh_inc = MIN(1,3/(double)N)/100;
		(*bubopts).lambda_N = 1;
		(*bubopts).use_weighting = 0;
		(*bubopts).opt_mesh_min = epsm;
		(*bubopts).opt_mesh_max = MIN(1,3/m)-epsm;
		(*bubopts).opt_mesh_inc = MIN(1,3/m)/100;
		(*bubopts).opt_use_weighting = 0;
		status = bub1(N,m,a,(*opts).bub_lambda_0,(*opts).bub_K,(*opts).bub_compat,bubopts);
	}
	 
	/* Get counts of counts */
	h = (int *)calloc(N+1,sizeof(int));
	for(i=0; i<(*in).C; i++)
		h[(int)(*in).wordcnt[i]]++;

	/* Then compute entropy */
	H = 0;
	for(j=0; j<=N; j++)
		H += a[j]*h[j];
	entropy->value = NAT2BIT(H);

	free(a);
	free(h);
	free(bubopts);

	return EXIT_SUCCESS; /*	TODO return status; reports errors with staverify.m - need to look into these errors and possibly pass them to message structure */
}

int bag1(int N,double m,double *a,double lambda_0,struct bub_options *bubopts)
{
	int meshsize;
	int i,j;
	double *p,*Y,**X,**B_jN,*A_data,**D,*I0,*IN;
	double *f,c;
	gsl_vector *b,*x,*tau,*norm,*S,*work;
	gsl_matrix *A,*V;
	int *signum;
	int d_status, s_status, status = EXIT_SUCCESS;

	meshsize = floor(((*bubopts).mesh_max-(*bubopts).mesh_min)/(*bubopts).mesh_inc)+1;

	A = gsl_matrix_calloc(meshsize+(N+1)+2,N+1);
	A_data = (*A).data;
	X = VectorToMatrixDouble(&A_data[0],meshsize,N+1);
	D = VectorToMatrixDouble(&A_data[meshsize*(N+1)],N+1,N+1);
	I0 = &A_data[(meshsize+(N+1))*(N+1)];
	IN = &A_data[(meshsize+(N+1)+1)*(N+1)];

	b = gsl_vector_calloc(meshsize+(N+1)+2);
	Y = (*b).data;

	/* Make mesh */
	p = (double *)malloc(meshsize*sizeof(double));
	for(i=0; i<meshsize; i++)
		p[i] = (*bubopts).mesh_min + i*(*bubopts).mesh_inc;
	
	/* Make weighting */
	f = (double *)malloc(meshsize*sizeof(double));
	make_weight(p,meshsize,m,(*bubopts).use_weighting,f,&c);

	/* Make H */
	for(i=0; i<meshsize; i++)
		Y[i] = f[i]*XLOGX(p[i]);
	
	/* Make binomial coefficents */
	B_jN = MatrixDouble(meshsize,N+1);
	make_binom(N,meshsize,p,B_jN);

	for(i=0; i<meshsize; i++)
		for(j=0; j<N+1; j++)
			X[i][j] = f[i]*B_jN[i][j];

	/* Make D */
	for(j=0; j<N; j++)
	{
		D[j][j]=-1*sqrt((double)N)/c;
		D[j][j+1]=1*sqrt((double)N)/c;
	}

	/* Make I0 and IN */
	I0[0] = sqrt(((double)N)*lambda_0)/c;
	IN[N] = sqrt(((double)N)*(*bubopts).lambda_N)/c;

	/* Solve the least squares problem */
	V = gsl_matrix_calloc(N+1,N+1);
	S = gsl_vector_calloc(N+1);
	work = gsl_vector_calloc(N+1);
	d_status = gsl_linalg_SV_decomp(A,V,S,work);
	if(d_status)
		status = EXIT_FAILURE;
	x = gsl_vector_calloc(N+1);
	s_status = gsl_linalg_SV_solve(A,V,S,b,x);
	if(s_status)
		status = EXIT_FAILURE;

	/* Read out the answer from x */
	memcpy(a,(*x).data,(N+1)*sizeof(double));

	/* Free memory */
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(x);

	FreeMatrixDouble(B_jN);
	free(X);
	free(D);
	free(f);
	free(p);
	gsl_matrix_free(A);
	gsl_vector_free(b);

	return status;
}

int bub1(int N,double m,double *a,double lambda_0,int K,int compat,struct bub_options *bubopts)
{
	int meshsize;
	int i,j,k;
	double *p,*Y,**X,**B_jN,*A_data,*b_data,**D,*I0,*IN;
	double *f,c;
	gsl_vector *b,*x,*tau,*norm,*S,*work;
	gsl_matrix *A,*V;
	int d_status, s_status, status = EXIT_SUCCESS;
	double best_error,max_error,h_MM;
	double *a_MM,*a_temp;
	int k1;

	/* Compute the a coefficients using the Miller-Madow method */
	a_MM = (double *)calloc(N+1,sizeof(double));
	for(j=0; j<=N; j++)
		a_MM[j] = XLOGX((double)j/(double)N) + (1-((double)j/(double)N))/(2*(double)N);

	meshsize = floor(((*bubopts).mesh_max-(*bubopts).mesh_min)/(*bubopts).mesh_inc)+1;

	/* Make mesh */
	p = (double *)malloc(meshsize*sizeof(double));
	for(i=0; i<meshsize; i++)
		p[i] = (*bubopts).mesh_min + i*(*bubopts).mesh_inc;
	
	/* Make weighting */
	f = (double *)malloc(meshsize*sizeof(double));
	make_weight(p,meshsize,m,(*bubopts).use_weighting,f,&c);

	B_jN = MatrixDouble(meshsize,N+1);
	make_binom(N,meshsize,p,B_jN);

	best_error = HUGE_VAL;

	for(k=0; k<MIN(K,N); k++)
	{
		k1 = k+1;

		a_temp = (double *)calloc(N+1,sizeof(double));
		memcpy(a_temp,a_MM,(N+1)*sizeof(double));

		/*** MAKE A ***/
		A = gsl_matrix_calloc(meshsize+k1+2,k1);
		A_data = (*A).data;
		X = VectorToMatrixDouble(&A_data[0],meshsize,k1);
		D = VectorToMatrixDouble(&A_data[meshsize*k1],k1,k1);
		I0 = &A_data[(meshsize+k1)*k1];
		IN = &A_data[(meshsize+k1+1)*k1];
		
		/*** MAKE B ***/
		b = gsl_vector_calloc(meshsize+k1+2);
		b_data = (*b).data;
		Y = &b_data[0];
		if(compat)      
			b_data[meshsize+k1+1] = (sqrt((double)N)/c)*a_MM[k];
		else
			b_data[meshsize+k1+1] = (sqrt((double)N)/c)*a_MM[k1];

		for(i=0; i<meshsize; i++)
		{
			h_MM = 0;
			for(j=k1; j<=N; j++)
				h_MM += B_jN[i][j]*a_MM[j];

			Y[i] = f[i]*(XLOGX(p[i]) - h_MM);
		}      
		
		for(i=0; i<meshsize; i++)
			for(j=0; j<k1; j++)
				X[i][j] = f[i]*B_jN[i][j];
		
		/* Make D */
		for(j=0; j<k1-1; j++)
		{
			D[j][j]=-1*sqrt((double)N)/c;
			D[j][j+1]=1*sqrt((double)N)/c;
		}
		
		/* Make I0 and IN */
		I0[0] = sqrt(((double)N)*lambda_0)/c;
		IN[k1-1] = sqrt(((double)N)*(*bubopts).lambda_N)/c;
		
		/* Solve the least squares problem */
		V = gsl_matrix_calloc(k1,k1);
		S = gsl_vector_calloc(k1);
		work = gsl_vector_calloc(k1);
		d_status = gsl_linalg_SV_decomp(A,V,S,work);
		if(d_status)
			status = EXIT_FAILURE;
		x = gsl_vector_calloc(k1);
		s_status = gsl_linalg_SV_solve(A,V,S,b,x);
		if(s_status)
			status = EXIT_FAILURE;
		
		/* Read out the answer from x */
		memcpy(a_temp,(*x).data,k1*sizeof(double));
		
#ifdef DEBUG
		printf("k1=%d\n",k1);
		for(j=0; j<=k; j++)
			printf("%f ",a_temp[j]);
		printf("\n");
#endif

		/* Find the error */
		max_error = compute_error(a_temp,N,m,compat,bubopts);

		if(max_error<best_error)
		{
			best_error = max_error;
			memcpy(a,a_temp,(N+1)*sizeof(double));
		}
				
		free(a_temp);

		/* Free memory */
		gsl_matrix_free(V);
		gsl_vector_free(S);
		gsl_vector_free(work);
		gsl_vector_free(x);
		free(X);
		free(D);
		gsl_matrix_free(A);
		gsl_vector_free(b);
	}
	
	free(f);
	free(p);
	FreeMatrixDouble(B_jN);
	free(a_MM);
}

void make_binom(int N,int meshsize,double *p,double **B_jN)
{
	int i,j;
	double *NCj;
	
	NCj = (double *)calloc(N+1,sizeof(double));

	for(j=0; j<N+1; j++)
		NCj[j] = gsl_sf_lngamma(N+1)-gsl_sf_lngamma(j+1)-gsl_sf_lngamma(N-j+1); 

	for(i=0; i<meshsize; i++)
		if(p[i]==0) 
		{
			B_jN[i][0]=1;
			for(j=1; j<N+1; j++)
				B_jN[i][j]=0;
		}
		else if(p[i]==1) 
		{
			for(j=0; j<N; j++)
				B_jN[i][j]=0;
			B_jN[i][N]=1;
		}
		else
		{
			for(j=0; j<N+1; j++)
				B_jN[i][j] = exp(NCj[j] + log(p[i])*j + log(1-p[i])*(N-j));
		}

	free(NCj);
}

void make_weight(double *p,int meshsize,double m,int use_weighting,double *f,double *c)
{
	int i;

	for(i=0; i<meshsize; i++)
		if(use_weighting && p[i]>1/m)
			f[i] = 1./p[i];
		else
			f[i] = m;

	if(use_weighting)
		*c = 2;
	else
		*c = 1;
}

double compute_error(double *a,int N,double m,int compat,struct bub_options *bubopts)
{
	int i,j;
	double max_bias,max_error,max_variance0,max_variance1;
	double *p,*f,**B_jN;
	double temp1,temp2,temp3,temp4,temp5;
	int meshsize;
	double c;

	meshsize = floor(((*bubopts).opt_mesh_max-(*bubopts).opt_mesh_min)/(*bubopts).opt_mesh_inc)+1;

	/* Make mesh */
	p = (double *)malloc(meshsize*sizeof(double));
	for(i=0; i<meshsize; i++)
		p[i] = (*bubopts).opt_mesh_min + i*(*bubopts).opt_mesh_inc;
	
	/* Make weighting */
	f = (double *)malloc(meshsize*sizeof(double));
	make_weight(p,meshsize,m,(*bubopts).opt_use_weighting,f,&c);

	/* Make B_jN */
	B_jN = MatrixDouble(meshsize,N+1);
	make_binom(N,meshsize,p,B_jN);

	/* Compute the worst-case bias */
	max_bias = 0;
	for(i=0; i<meshsize; i++)
	{
		temp1 = 0;
		for(j=0; j<N+1; j++)
			temp1+=B_jN[i][j]*a[j];
		temp2 = f[i]*fabs(temp1 - XLOGX(p[i]));
		if(temp2>max_bias)
			max_bias = temp2;
	}

	/* Steele loose bound on variance */
	max_variance1 = 0;
	for(i=0; i<meshsize; i++)
	{
		temp3 = 0;
		for(j=1; j<N+1; j++)
			temp3 += B_jN[i][j]*((double)j/(double)N)*pow(a[j]-a[j-1],2);
		temp4 = f[i]*temp3;
		if(temp4>max_variance1)
			max_variance1 = temp4;
	}
	if(compat)
		max_variance1 *= 4*N*c;
	else
		max_variance1 *= 2*N*c;

	/* McDiarmid bound on variance */
	max_variance0 = 0;
	for(j=1; j<N+1; j++)
	{
		temp5 = fabs(a[j]-a[j-1]);
		if(temp5>max_variance0)
			max_variance0=temp5;
	}
	max_variance0 = N*pow(max_variance0,2);

#ifdef DEBUG
	printf("maxbias=%f maxvariance0=%f maxvariance1=%f\n",max_bias,max_variance0,max_variance1);
#endif

	max_error = NAT2BIT(sqrt(pow(max_bias,2) + MIN(max_variance0,max_variance1)));
	free(p);
	free(f);
	FreeMatrixDouble(B_jN);

	return max_error;
}

