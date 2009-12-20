/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */
/* Memory allocation functions */
extern double **VectorToMatrixDouble(double *in,int M,int N);
extern int **MatrixInt(int M,int N);
extern int ***Matrix3Int(int M,int N,int P);
extern double **MatrixDouble(int M,int N);
extern double ***Matrix3Double(int M,int N,int P);
extern double ****Matrix4Double(int M,int N,int P,int Q);
extern void FreeMatrixInt(int **in);
extern void FreeMatrix3Int(int ***in);
extern void FreeMatrixDouble(double **in);
extern void FreeMatrix3Double(double ***in);
extern void FreeMatrix4Double(double ****in);

extern int MaxMatrix3Int(int ***x,int M,int P,int N);
extern int MaxMatrixInt(int **x,int P,int N);
extern int MaxVectorInt(int *x,int P);
extern double MinVectorDouble(double *x,int P);

/* Functions for handling spike trains */
extern void GetLims(double *times,double t_start,double t_end,int Q,int *q_start,int *q_end);

