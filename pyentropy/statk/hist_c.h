/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */
struct hist1dvec{
  struct hist1d *vec;
  int M;
  int P;
  struct estimate *entropy;
};

struct hist2d{
  struct hist1d *joint;
  struct hist1d *row;
  struct hist1d *col;
  struct estimate *information;
};

struct histcond{
  struct hist1dvec *classcond;
  struct hist1d *total;
  struct estimate *information;
};

/* Histograms and estimates */
extern struct hist1d *CAllocHist1D(int M,int *P_vec,int N,struct options_entropy *opts);
extern void CFreeHist1D(int M,struct hist1d *hist_c,struct options_entropy *opts);
extern struct hist2d *CAllocHist2D(int P_total,int N,struct options_entropy *opts);
extern void CFreeHist2D(struct hist2d *hist_c,struct options_entropy *opts);
extern struct histcond *CAllocHistCond(int M,int P_total,int *P_vec,int N,struct options_entropy *opts);
extern void CFreeHistCond(int M,struct histcond *hist_c,struct options_entropy *opts);
extern struct hist1dvec *CAllocHist1DVec(int M,int *P_vec,int N,struct options_entropy *opts);
extern void CFreeHist1DVec(int M,struct hist1dvec *hist_c,struct options_entropy *opts);
extern struct estimate *CAllocEst(struct options_entropy *opts);
extern void CFreeEst(struct estimate *in, struct options_entropy *opts);

extern void AddEst(struct estimate *A,struct estimate *B,struct estimate *C,struct options_entropy *opts);
extern void IncEst(struct estimate *A,struct estimate *B,struct options_entropy *opts);
extern void IncScalarEst(struct estimate *A,double B,struct options_entropy *opts);
extern void SubtractEst(struct estimate *A,struct estimate *B,struct estimate *C,struct options_entropy *opts);
extern void ScaleEst(struct estimate *A,double S,struct estimate *C,struct options_entropy *opts);
extern void ZeroEst(struct estimate *C,struct options_entropy *opts);

extern int MatrixToHist2DComp(double **table,
			   int M, /* rows in table input */
			   int N, /* cols in table input */
			   struct hist2d *out, /* output structure */
			   struct options_entropy *opts /* entropy options */
			   );
extern int MarginalProc(int C_in,int N,int **list,double *cnt,int **wordlist,double *wordcnt);
extern int Info2DComp(struct hist2d *in,struct options_entropy *opts);
extern int InfoCondComp(struct histcond *in,struct options_entropy *opts);
extern int Entropy1DVecComp(struct hist1dvec *in,struct options_entropy *opts);
