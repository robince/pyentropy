/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */
#include "toolkit_c.h"

int **MatrixInt(int M,int N)
{
  int **out;
  int m;

  out = (int **)calloc(M,sizeof(int *));
  out[0] = (int *)calloc(M*N,sizeof(int));
  for(m=1;m<M;m++)
    out[m] = out[m-1]+N;

  return out;
}

double **MatrixDouble(int M,int N)
{
  double **out;
  int m;

  out = (double **)calloc(M,sizeof(double *));
  out[0] = (double *)calloc(M*N,sizeof(double));
  for(m=1;m<M;m++)
    out[m] = out[m-1]+N;

  return out;
}

double **VectorToMatrixDouble(double *in,int M,int N)
{
  double **out;
  int m;

  out = (double **)calloc(M,sizeof(double *));
  out[0] = in;
  for(m=1;m<M;m++)
    out[m] = out[m-1]+N;

  return out;
}

int ***Matrix3Int(int M,int N,int P)
{
  int ***out;
  int i,j;

  /* These are pointers to the start of M+1 rows */
  out = (int ***)calloc(M,sizeof(int **));
  
  /* These are pointers to the start of M+1xN+1 pokes */
  out[0] = (int **)calloc(M*N,sizeof(int *));
  
  /* This is M+1xN+1xP pointers to ints */
  out[0][0] = (int *)calloc(M*N*P,sizeof(int));
  
  /* This sets pointers to the start of the rows */
  for(i=1;i<M;i++)
    {
      out[i]=out[i-1]+N;
      out[i][0]=out[i-1][0]+N*P;
    }
  
  /* This sets pointers to the start of the pokes */
  for(i=0;i<M;i++)
    for(j=1;j<N;j++)
      out[i][j]=out[i][j-1]+P;

  return out;
}

double ***Matrix3Double(int M,int N,int P)
{
  double ***out;
  int i,j;

  /* These are pointers to the start of M+1 rows */
  out = (double ***)calloc(M,sizeof(double **));
  
  /* These are pointers to the start of M+1xN+1 pokes */
  out[0] = (double **)calloc(M*N,sizeof(double *));
  
  /* This is M+1xN+1xP pointers to doubles */
  out[0][0] = (double *)calloc(M*N*P,sizeof(double));
  
  /* This sets pointers to the start of the rows */
  for(i=1;i<M;i++)
    {
      out[i]=out[i-1]+N;
      out[i][0]=out[i-1][0]+N*P;
    }
  
  /* This sets pointers to the start of the pokes */
  for(i=0;i<M;i++)
    for(j=1;j<N;j++)
      out[i][j]=out[i][j-1]+P;

  return out;
}

double ****Matrix4Double(int M,int N,int P,int Q)
{
  double ****out;
  int i,j,k;

  /* These are pointers to the start of M+1 rows */
  out = (double ****)calloc(M,sizeof(double ***));
  
  /* These are pointers to the start of M+1xN+1 pokes */
  out[0] = (double ***)calloc(M*N,sizeof(double **));
  
  /* This is M+1xN+1xP pointers to doubles */
  out[0][0] = (double **)calloc(M*N*P,sizeof(double *));
  
  /* blah */
  out[0][0][0] = (double *)calloc(M*N*P*Q,sizeof(double));
  
  /* This sets pointers to the start of the rows */
  for(i=1;i<M;i++)
    {
      out[i]=out[i-1]+N;
      out[i][0]=out[i-1][0]+N*P;
      out[i][0][0]=out[i-1][0][0]+N*P*Q;
    }
  
  /* This sets pointers to the start of the pokes */
  for(i=0;i<M;i++)
    for(j=1;j<N;j++)
      {
	out[i][j]=out[i][j-1]+P;
	out[i][j][0]=out[i][j-1][0]+P*Q;
      }

  for(i=0;i<M;i++)
    for(j=0;j<N;j++)
      for(k=1;k<P;k++)
	out[i][j][k]=out[i][j][k-1]+Q;

  return out;
}

void FreeMatrixInt(int **in)
{
  free(in[0]);
  free(in);
}

void FreeMatrixDouble(double **in)
{
  free(in[0]);
  free(in);
}

void FreeMatrix3Double(double ***in)
{
  free(in[0][0]);
  free(in[0]);
  free(in);
}

void FreeMatrix4Double(double ****in)
{
  free(in[0][0][0]);
  free(in[0][0]);
  free(in[0]);
  free(in);
}

void FreeMatrix3Int(int ***in)
{
  free(in[0][0]);
  free(in[0]);
  free(in);
}

/* This function returns the indices of the first spike and the last
   spike given the start time and the end time */
void GetLims(double *times,double t_start,double t_end,int Q,int *q_start,int *q_end)
{
  int start_flag,end_flag;
  int q;
  double cur_time;

  start_flag = 0;
  end_flag = 0;
  
  if(Q==0)
    {
      *q_start = 0;
      *q_end = -1;
    }
  else
    {
      /* Get the index of the first spike and the last spike */
      for(q=0;q<Q;q++)
	{
	  cur_time = times[q];
	  /* If we hit a spike after the end time before
	     hitting a spike after the start time,
	     then there are no spikes in the range. */
	  if((cur_time>=t_end) & (start_flag==0))
	    {
	      *q_start = 0;
	      *q_end = -1;
	      start_flag = 1;
	      end_flag = 1;
	    }
	  /* If we hit a spike after the start time
	     (but before the end time)
	     and we haven't hit the start spike yet,
	     then this is the start spike. */
	  else if((cur_time>=t_start) & (start_flag==0))
	    {
	      *q_start = q;
	      start_flag = 1;
	    }
	  /* If we hit a spike after the end time,
	     and we haven't hit the end spike yet,
	     then this is the end spike. */
	  else if((cur_time>=t_end) & (end_flag==0))
	    {
	      *q_end = q-1;
	      end_flag = 1;
	    }
	  /* Otherwise, the spike is neither the start spike
	     or the end spike. */
	}

      /* If we get to the end,
	 and we haven't raised the start flag,
	 then there are no spikes in the range */
      if(start_flag==0)
	{
	  *q_start = 0;
	  *q_end = -1;
	}
      /* If we get to the end,
	 and we haven't raised the end flag,
	 then raise the end flag. */
      else if(end_flag==0)
	*q_end = Q-1;
    }
}

int MaxMatrix3Int(int ***x,int M,int P,int N)
{
  int n,p,m;
  int y;

  y = x[0][0][0];

  for (m=0;m<M;m++)
    for (p=0;p<P;p++)
      for (n=0;n<N;n++)
	if(x[m][p][n]>y)
	  y=x[m][p][n];
  return y;
}

int MaxMatrixInt(int **x,int P,int N)
{
  int n,p;
  int y;

  y = x[0][0];

  for (p=0;p<P;p++)
    for (n=0;n<N;n++)
      if(x[p][n]>y)
	y=x[p][n];
  return y;
}

int MaxVectorInt(int *x,int P)
{
  int p;
  int y;

  y = x[0];

  for (p=0;p<P;p++)
    if(x[p]>y)
      y=x[p];
  return y;
}

double MinVectorDouble(double *x,int P)
{
  int p;
  double y;

  y = x[0];

  for (p=0;p<P;p++)
    if(x[p]<y)
      y=x[p];
  return y;
}
