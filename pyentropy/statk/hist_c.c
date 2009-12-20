/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

#include "toolkit_c.h"

extern char ent_est_meth_list[ENT_EST_METHS][MAXCHARS];
extern char var_est_meth_list[GEN_VAR_EST_METHS+SPEC_VAR_EST_METHS][MAXCHARS];

/**************************************************/
/************* Memory allocation functions ********/
/**************************************************/

/*** hist1d ***/

struct hist1d *CAllocHist1D(int M,int *P_vec,int N,struct options_entropy *opts)
{
  int m;
  struct hist1d *hist_c;

  hist_c = (struct hist1d *)malloc(M*sizeof(struct hist1d));

  for(m=0;m<M;m++)
    {
      hist_c[m].wordlist = MatrixInt(P_vec[m],N);
      hist_c[m].wordcnt = (double *)calloc(P_vec[m],sizeof(double));
      hist_c[m].entropy = CAllocEst(opts);
    }

  return hist_c;
}

void CFreeHist1D(int M,struct hist1d *hist_c,struct options_entropy *opts)
{
  int m;

  for(m=0;m<M;m++)
    {
      FreeMatrixInt(hist_c[m].wordlist);
      free(hist_c[m].wordcnt);
      CFreeEst(hist_c[m].entropy,opts);
    }

  free(hist_c);
}

/*** hist2d ***/

struct hist2d *CAllocHist2D(int P_total,int N,struct options_entropy *opts)
{
  struct hist2d *hist_c;

  hist_c = (struct hist2d *)malloc(sizeof(struct hist2d));

  (*hist_c).joint = CAllocHist1D(1,&P_total,N,opts);
  (*hist_c).row = CAllocHist1D(1,&P_total,N,opts);
  (*hist_c).col = CAllocHist1D(1,&P_total,N,opts);
  (*hist_c).information = CAllocEst(opts);

  return hist_c;
}

void CFreeHist2D(struct hist2d *hist_c,struct options_entropy *opts)
{
  CFreeHist1D(1,(*hist_c).joint,opts);
  CFreeHist1D(1,(*hist_c).row,opts);
  CFreeHist1D(1,(*hist_c).col,opts);
  CFreeEst((*hist_c).information,opts);

  free(hist_c);
}

/*** histcond ***/

struct histcond *CAllocHistCond(int M,int P_total,int *P_vec,int N,struct options_entropy *opts)
{
  struct histcond *hist_c;

  hist_c = (struct histcond *)malloc(sizeof(struct histcond));

  (*hist_c).classcond = CAllocHist1DVec(M,P_vec,N,opts);
  (*hist_c).total = CAllocHist1D(1,&P_total,N,opts);
  (*hist_c).information = CAllocEst(opts);

  return hist_c;
}

void CFreeHistCond(int M,struct histcond *hist_c,struct options_entropy *opts)
{
  CFreeHist1DVec(M,(*hist_c).classcond,opts);
  CFreeHist1D(1,(*hist_c).total,opts);
  CFreeEst((*hist_c).information,opts);

  free(hist_c);
}

/*** hist1dvec ***/

struct hist1dvec *CAllocHist1DVec(int M,int *P_vec,int N,struct options_entropy *opts)
{
  struct hist1dvec *hist_c;

  hist_c = (struct hist1dvec *)malloc(sizeof(struct hist1dvec));

  (*hist_c).vec = CAllocHist1D(M,P_vec,N,opts);
  (*hist_c).entropy = CAllocEst(opts);

  return hist_c;
}

void CFreeHist1DVec(int M,struct hist1dvec *hist_c,struct options_entropy *opts)
{
  CFreeHist1D(M,(*hist_c).vec,opts);
  CFreeEst((*hist_c).entropy,opts);

  free(hist_c);
}

/*** est ***/

struct estimate *CAllocEst(struct options_entropy *opts)
{
  int e,v;
  struct estimate *in;

  in = (struct estimate *)malloc((*opts).E*sizeof(struct estimate));
  
  for(e=0;e<(*opts).E;e++)
    {
      strcpy(in[e].name,"temp");

      /* add the messages structure */
      in[e].messages = (struct message *)malloc(sizeof(struct message));
      in[e].messages->i = in[e].messages->j = in[e].messages->k = 0;
      in[e].messages->status = in[e].messages->warnings = in[e].messages->errors = NULL;

      /* add the extras structure */
      in[e].E = 0;
      in[e].extras = NULL;

      /* add the variance estimate structure */
      in[e].V = opts->V[e];
      in[e].ve = (struct nv_pair *)malloc(in[e].V*sizeof(struct nv_pair));
      for(v=0;v<in[e].V;v++)
	memcpy(in[e].ve[v].name,var_est_meth_list[(*opts).var_est_meth[e][v]-1],MAXCHARS*sizeof(char));
    }

  return in;
}

void CFreeEst(struct estimate *in,struct options_entropy *opts)
{
  int e;

  for(e=0;e<(*opts).E;e++)
    {
      if(in[e].messages->i)
        free(in[e].messages->status);
      if(in[e].messages->j)
        free(in[e].messages->warnings);
      if(in[e].messages->k)
        free(in[e].messages->errors);
      free(in[e].messages);

      if(in[e].E)
	free(in[e].extras);

      if(in[e].V)
	free(in[e].ve);
    }
  
  free(in);
}


/************************************************/
/****************** Estimate arithmetic *********/
/************************************************/

void AddEst(struct estimate *A,struct estimate *B,struct estimate *C,struct options_entropy *opts)
{
  int e,v;

  for(e=0;e<(*opts).E;e++)
    {
      C[e].value = A[e].value + B[e].value;

      /* scroll through the variance estimate requests */
      for(v=0;v<(*opts).V[e];v++)
	C[e].ve[v].value = A[e].ve[v].value + B[e].ve[v].value;
    }
}

void IncEst(struct estimate *A,struct estimate *B,struct options_entropy *opts)
{
  int e,v;

  for(e=0;e<(*opts).E;e++)
    {
      A[e].value += B[e].value;
      
      /* scroll through the variance estimate requests */
      for(v=0;v<(*opts).V[e];v++)
	A[e].ve[v].value += B[e].ve[v].value;
    }
}

void IncScalarEst(struct estimate *A,double B,struct options_entropy *opts)
{
  int e;

  for(e=0;e<(*opts).E;e++)
    A[e].value += B;
}

void SubtractEst(struct estimate *A,struct estimate *B,struct estimate *C,struct options_entropy *opts)
{
  int e,v;

  for(e=0;e<(*opts).E;e++)
    {
      C[e].value = A[e].value - B[e].value;

      /* scroll through the variance estimate requests */
      for(v=0;v<(*opts).V[e];v++)
	C[e].ve[v].value = A[e].ve[v].value + B[e].ve[v].value;
    }
}

void ScaleEst(struct estimate *A,double S,struct estimate *C,struct options_entropy *opts)
{
  int e,v;

  for(e=0;e<(*opts).E;e++)
    {
      C[e].value = S*(A[e].value);
      
      /* scroll through the variance estimate requests */
      for(v=0;v<(*opts).V[e];v++)
	C[e].ve[v].value = S*(A[e].ve[v].value);
    }
}

void ZeroEst(struct estimate *C,struct options_entropy *opts)
{
  int e,v;

  for(e=0;e<(*opts).E;e++)
    {
      C[e].value = 0;
      
      /* scroll through the variance estimate requests */
      for(v=0;v<(*opts).V[e];v++)
	C[e].ve[v].value = 0;
    }
}
