/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

/** @file
 * @brief Sorting routines for C-style arrays.
 * This file contains C code used to sort C-style arrays of numbers.
 */

#include "toolkit_c.h"

/* Defines for dgsort() */
#define SWAP(a,b) {itemp=a; a=b; b=itemp;}
#define THRESH 7
#define NSTACK 50

/**
 * @brief Checks whether an array of doubles is sorted.
 * Returns an integer to indicate whether an array of doubles is (1) or is not
 * (0) sorted from greatest to least.
 */
int IsSortedDouble(int M,double *in)
{
  int m,out;

  out = 1;
  for(m=1;m<M;m++)
    if(in[m]<in[m-1])
      {
	out = 0;
	break;
      }
  
  return out;
}

int UniqueInt(int M,int *list,int *uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int *sorted;
  int *sort_idx;
  int u;

  sorted = (int *)malloc(M*sizeof(int));
  sort_idx = (int *)malloc(M*sizeof(int));

  SortInt(M,list,sorted,sort_idx);

  u=CountSortedInt(M,sorted,sort_idx,uni_list,uni_i,uni_j,cnt);
  
  free(sorted);
  free(sort_idx);
  
  return u;
}

int CountSortedInt(int M,int *sorted,int *sort_idx,int *uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int u,m;
  
  /* handle the first item specially */
  uni_list[0] = sorted[0];
  uni_i[0] = sort_idx[0];
  cnt[0]++;
  uni_j[sort_idx[0]] = 0;

  /* scroll down the list */
  u=0;
  for(m=1;m<M;m++)
    {
      if(sorted[m]!=sorted[m-1])
	{
	  u++;
	  uni_list[u] = sorted[m];
	  uni_i[u] = sort_idx[m];
	}
      cnt[u]++;
      uni_j[sort_idx[m]] = u;
    }
  u++;

  return u;
}

/**
 * @brief Returns the number of unique rows in an array.
 * Uses a sorting algorithm to determine the number of unique rows in the array
 * list, and places the list in uni_list. Also sets the sort indexes uni_i and
 * uni_j, and the row counts in cnt.
 * @param[in] M The number of rows in sorted.
 * @param[in] N The number of columns in sorted.
 * @param[in] list The array in which to count the unique rows.
 * @param[out] uni_list An array containing the unique rows of list.
 * @param[out] uni_i An index to relate uni_list to list (i.e.,
 * uni_list[i]==unsorted[uni_i[i]]).
 * @param[out] uni_j An index to relate list to uni_list (i.e.,
 * unsorted[i]==uni_list[uni_j[i]]).
 * @param[out] cnt An array of the row counts.
 */
int UniqueRowsInt(int M,int N,int **list,int **uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int **sorted;
  int *sort_idx;
  int u;

  /* First, sort the numbers */

  sorted = MatrixInt(M,N);
  sort_idx = (int *)malloc(M*sizeof(int));

  SortRowsInt(M,N,list,sorted,sort_idx);

  u=CountSortedRowsInt(M,N,sorted,sort_idx,uni_list,uni_i,uni_j,cnt);
  
  FreeMatrixInt(sorted);
  free(sort_idx);
  
  return u;
}

/**
 * @brief Count the unique rows in a sorted array.
 * Returns the number of unique rows in the array sorted, and places the list
 * of unique rows in uni_list. Also sets the sort indexes uni_i and uni_j, and
 * the row counts in cnt.
 * @param[in] M The number of rows in sorted.
 * @param[in] N The number of columns in sorted.
 * @param[in] sorted The sorted array in which to count the unique rows.
 * @param[in] sort_idx The sorting index that specifies from which row of the
 * unsorted array each row of sorted derives (i.e.,
 * sorted[i]==unsorted[sort_idx[i]]).
 * @param[out] uni_list An array containing the unique rows of sorted.
 * @param[out] uni_i An index to relate uni_list to the unsorted array (i.e.,
 * uni_list[i]==unsorted[uni_i[i]]).
 * @param[out] uni_j An index to relate the unsorted array to uni_list (i.e.,
 * unsorted[i]==uni_list[uni_j[i]]).
 * @param[out] cnt An array of the row counts.
 */
int CountSortedRowsInt(int M,int N,int **sorted,int *sort_idx,int **uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int m,u,n1,n;

  /* handle the first item specially */
  for(n1=0;n1<N;n1++)  
    uni_list[0][n1] = sorted[0][n1];
  uni_i[0] = sort_idx[0];
  cnt[0]++;
  uni_j[sort_idx[0]] = 0;

  /* scroll down the list */
  u=0;
  for(m=1;m<M;m++)
    {
      /* for each letter */
      for(n=0;n<N;n++)
	if(sorted[m][n]!=sorted[m-1][n])
	  break;

      /* is this item different from the previous item? */
      if(n<N)
	{
	  u++;
	  for(n1=0;n1<N;n1++)
	    uni_list[u][n1]=sorted[m][n1];
	  uni_i[u] = sort_idx[m];
	}
      cnt[u]++;
      uni_j[sort_idx[m]] = u;
    }
  u++;
  
  return u;
}

/* The calling sequence is a bit different from that of DirectCountComp */
void SortRowsInt(int P,int N,int **list,int **sorted,int *p_list)
{
  int p,n;
  int *col,*col_sort;
  int *idx,*p_list_temp;

  idx=(int *)malloc(P*sizeof(int));
  col = (int *)malloc(P*sizeof(int));
  col_sort = (int *)malloc(P*sizeof(int));
  p_list_temp = (int *)malloc(P*sizeof(int));

  /* initialize n */
  for(p=0;p<P;p++)
    p_list[p]=p;

  /* for each column, starting from the back */
  for (n=N-1;n>=0;n--)
    {
      /* make column */
      for(p=0;p<P;p++)
	col[p]=list[p_list[p]][n];
 
      /* We must have a stable sort:
	 In the case of a tie, the original order must be preserved! */
      SortInt(P,col,col_sort,idx);
 
      for(p=0;p<P;p++)
	p_list_temp[p] = p_list[idx[p]];
      
      for(p=0;p<P;p++)
	p_list[p]=p_list_temp[p];
    }
  
  for (p=0;p<P;p++)
    for (n=0;n<N;n++)
      sorted[p][n]=list[p_list[p]][n];
  
  free(idx);
  free(col);
  free(col_sort);
  free(p_list_temp);
}

void SortUnstableInt(int N, int *B, int *C, int *idx)
{
  int l,r;
  int i,j,k;
  int vidx,itemp;
  int stack_idx,*stack;
  int v;
  int *A;

  A = (int *)malloc(N*sizeof(int));

  /* Initialize the index vector */
  /* FOR A STABLE SORT, USE THE INDEX AS A TIE-BREAKER! */
  /* Modify the key vector so that each element is distinct */
  /* and the original index is the tie-breaker */
  for (j=0;j<N;j++)
    {
      idx[j]=j;
      A[j]=B[j];
    }

  /* Initialize the stack */
  stack = (int *)malloc(NSTACK*sizeof(int));
  stack_idx=0;
  stack[stack_idx++]=0; /* initialize the left pointer to zero */
  stack[stack_idx++]=N-1; /* initialize the right pointer to N-1 */
  
  while(stack_idx>0) /* loop until done */
    {
      /* Pop the stack and begin a new round */
      r = stack[--stack_idx];
      l = stack[--stack_idx];

      /* If the number of remaining elements is above the threshold, 
	 we do a quicksort */
      if(r-l>=THRESH)
	{
	  /*******************************************************/
	  /* Median-of-three modification */
	  /*******************************************************/

	  /* Choose the median of the left (l), center (k), and right (r)
	     elements as the partitioning element */
	  k=(l+r)>>1; /* find the center element */
	  SWAP(idx[k],idx[r]);
	  if(A[idx[l]] > A[idx[r]]) SWAP(idx[l],idx[r]);
	  if(A[idx[l+1]] > A[idx[r]]) SWAP(idx[l+1],idx[r]);
	  if(A[idx[l]] > A[idx[l+1]]) SWAP(idx[l],idx[l+1]);

	  /*******************************************************/
	  /* Partition the subfile */
	  /*******************************************************/

	  i=l+1; /* Position of partitioning element */
	  j=r;
	  vidx=idx[i]; /* Index of partitioning element */
	  v=A[vidx]; /* Value of partitioning element */
	  
	  /* While the pointers haven't yet crossed */
	  while(i <= j)
	    {
	      while (A[idx[++i]] < v); /* Scan up from l+1 to find element > a */
	      while (A[idx[--j]] > v); /* Scan down from r to find element < a */
	      if(i <= j)
		SWAP(idx[i],idx[j]); /* Pointer didn't cross, exchange elements */
	    }
	  /* At the end, i=j+1 */

	  /* Insert partitioning element */
	  idx[l+1]=idx[j]; 
	  idx[j]=vidx;

	  /*******************************************************/
	  /* Push the pointers to the larger subfile on the stack,
	     process the smaller subfile now */
	  /*******************************************************/

	  if(stack_idx+4>NSTACK)
	    printf("Imminent stack overflow!\n");

	  /* Push the pointers for the right subfile on the stack,
	     sort the left subfile immediately */
	  if(r-i+1 >= j-l)
	    {
	      stack[stack_idx++] = i;
	      stack[stack_idx++] = r;
	      stack[stack_idx++] = l;
	      stack[stack_idx++] = j-1;
	    }
	  /* Push the pointers for the left subfile on the stack,
	     sort the right subfile immediately */
	  else
	    {
	      stack[stack_idx++] = l;
	      stack[stack_idx++] = j-1;
	      stack[stack_idx++] = i;
	      stack[stack_idx++] = r;
	    }
	}

      /*******************************************************/
      /* If the number of remaining elements is below the threshold,
	 do an insertion sort */
      /*******************************************************/
      else
	{
	  /* For each element from l to r */
	  for(j=l+1;j<=r;j++)
	    {
	      /* Pick out the current element v */
	      vidx=idx[j];
	      v=A[vidx];
	      
	      /* Find the place to insert v */
	      /* Start from the element to the left of v and proceed leftward */
	      for (i=j-1;i>=l;i--)
		{
		  /* If the current element is less than or equal to v,
		     then v is in the proper place */
		  if (A[idx[i]]<=v)
		    break; 
		  
		  /* Otherwise, shift all of the elements over and keep looking */
		  idx[i+1]=idx[i];
		}
      
	      /* Once we've found the proper place for v, set its index */ 
	      /* Since i is the index of the element to the left of v, i+1 is the index of v */
	      /* If we complete the loop, then v goes at the front of the list */
	      /* At that point, i=l-1, which makes i+1=l */
	      idx[i+1]=vidx;
	    }
	}
    }

  /* Sort the list according to the indices */
  for (j=0;j<N;j++)
    C[j]=B[idx[j]];

  free(A);
  free(stack);
}

void SortInt(int N, int *A, int *C, int *idx)
{
  int u,n,j,U;
  int *run_length,*run_start,*run_idx,*run_sort;
  int run_flag;
  int *idx_temp;

  run_length = (int *)calloc(N,sizeof(int));
  run_start = (int *)malloc(N*sizeof(int));

  /* First do an unstable sort */
  SortUnstableInt(N,A,C,idx);

  /* Next, scan the sorted list for runs of equal values */
  u=0;
  run_flag=0;
  for(n=1;n<N;n++)
    /* If this element is equal to the previous, we are in a run */
    if(C[n]==C[n-1])

      /* If we are in the middle of a run */
      if(run_flag)
	run_length[u-1]++; /* Increment the length of the run */

      /* Else, a run is just beginning */
      else
	{
	  /* Get the starting point of the run */
	  u++;
	  run_start[u-1] = n-1;
	  run_length[u-1] = 2;
	  run_flag=1;
	}
    /* Otherwise there is no run or the run has just ended */
    else
      run_flag=0;

  U=u;

  idx_temp = (int *)malloc(N*sizeof(int));
  memcpy(idx_temp,idx,N*sizeof(int));

  /* Sort the indices in each run */
  for(u=0;u<U;u++)
    {
      run_idx = (int *)malloc(run_length[u]*sizeof(int));
      run_sort = (int *)malloc(run_length[u]*sizeof(int));

      SortInt(run_length[u],&idx[run_start[u]],run_sort,run_idx);
 
     /* Reorder the indices */
      for(j=0;j<run_length[u];j++)
	idx[run_start[u]+j] = idx_temp[run_start[u]+run_idx[j]]; 
    
      free(run_idx);
      free(run_sort);
    }
  
  free(idx_temp);
  free(run_start);
  free(run_length);
}

/**
 * @brief Returns the number of unique doubles in a list.
 * Uses a sorting algorithm to determine the number of unique doubles in the
 * list list, and places the list in uni_list. Also sets the sort indexes uni_i
 * and uni_j, and the counts in cnt.
 * @param[in] M The number of rows in list.
 * @param[in] list The list in which to count the unique doubles.
 * @param[out] uni_list A list containing the unique doubles of list.
 * @param[out] uni_i An index to relate uni_list to list (i.e.,
 * uni_list[i]==unsorted[uni_i[i]]).
 * @param[out] uni_j An index to relate list to uni_list (i.e.,
 * unsorted[i]==uni_list[uni_j[i]]).
 * @param[out] cnt An array of the counts.
 */
int UniqueDouble(int M,double *list,double *uni_list,int *uni_i,int *uni_j,int *cnt)
{
  double *sorted;
  int *sort_idx;
  int u;

  /* First, sort the numbers */

  sorted = (double *)malloc(M*sizeof(double));
  sort_idx = (int *)malloc(M*sizeof(int));

  /* A stable sort is not required here */
  SortUnstableDouble(M,list,sorted,sort_idx);

  u=CountSortedDouble(M,sorted,sort_idx,uni_list,uni_i,uni_j,cnt);

  free(sorted);
  free(sort_idx);
  
  return u;
}

int CountSortedDouble(int M,double *sorted,int *sort_idx,double *uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int u,m;

  /* handle the first item specially */
  uni_list[0] = sorted[0];
  uni_i[0] = sort_idx[0];
  cnt[0]++;
  uni_j[sort_idx[0]] = 0;

  /* scroll down the list */
  u=0;
  for(m=1;m<M;m++)
    {
      if(sorted[m]!=sorted[m-1])
	{
	  u++;
	  uni_list[u] = sorted[m];
	  uni_i[u] = sort_idx[m];
	}
      cnt[u]++;
      uni_j[sort_idx[m]] = u;
    }
  u++;
  
  return u;
}

int UniqueRowsDouble(int M,int N,double **list,double **uni_list,int *uni_i,int *uni_j,int *cnt)
{
  double **sorted;
  int *sort_idx;
  int u;

  sorted = MatrixDouble(M,N);
  sort_idx = (int *)malloc(M*sizeof(int));

  SortRowsDouble(M,N,list,sorted,sort_idx);

  u=CountSortedRowsDouble(M,N,sorted,sort_idx,uni_list,uni_i,uni_j,cnt);

  FreeMatrixDouble(sorted);
  free(sort_idx);
  
  return u;
}

/* The calling sequence is a bit different from that of DirectCountComp */
void SortRowsDouble(int P,int N,double **list,double **sorted,int *p_list)
{
  int p,n;
  double *col,*col_sort;
  int *idx,*p_list_temp;

  idx=(int *)malloc(P*sizeof(int));
  col = (double *)malloc(P*sizeof(double));
  col_sort = (double *)malloc(P*sizeof(double));
  p_list_temp = (int *)malloc(P*sizeof(int));

  /* initialize n */
  for(p=0;p<P;p++)
    p_list[p]=p;

  /* for each column, starting from the back */
  for (n=N-1;n>=0;n--)
    {
      /* make column */
      for(p=0;p<P;p++)
	col[p]=list[p_list[p]][n];
 
      /* We must have a stable sort:
	 In the case of a tie, the original order must be preserved! */
      SortDouble(P,col,col_sort,idx);
 
      for(p=0;p<P;p++)
	p_list_temp[p] = p_list[idx[p]];
      
      for(p=0;p<P;p++)
	p_list[p]=p_list_temp[p];
    }
  
  for (p=0;p<P;p++)
    for (n=0;n<N;n++)
      sorted[p][n]=list[p_list[p]][n];
  
  free(idx);
  free(col);
  free(col_sort);
  free(p_list_temp);
}

int CountSortedRowsDouble(int M,int N,double **sorted,int *sort_idx,double **uni_list,int *uni_i,int *uni_j,int *cnt)
{
  int m,u,n1,n;

  /* handle the first item specially */
  for(n1=0;n1<N;n1++)  
    uni_list[0][n1] = sorted[0][n1];
  uni_i[0] = sort_idx[0];
  cnt[0]++;
  uni_j[sort_idx[0]] = 0;

  /* scroll down the list */
  u=0;
  for(m=1;m<M;m++)
    {
      /* for each letter */
      for(n=0;n<N;n++)
	if(sorted[m][n]!=sorted[m-1][n])
	  break;

      /* is this item different from the previous item? */
      if(n<N)
	{
	  u++;
	  for(n1=0;n1<N;n1++)
	    uni_list[u][n1]=sorted[m][n1];
	  uni_i[u] = sort_idx[m];
	}
      cnt[u]++;
      uni_j[sort_idx[m]] = u;
    }
  u++;
  
  return u;
}

/* Unstable quicksort algorithm */
/* Based on Sedgewick, R. "Implementing Quicksort Programs",
   Helpful explanations in "Numerical Recipes in C" */ 
void SortUnstableDouble(int N, double *A, double *C, int *idx)
{
  int l,r;
  int i,j,k;
  int vidx,itemp;
  int stack_idx,*stack;
  double v;

  /* Initialize the index vector */
  for (j=0;j<N;j++)
    idx[j]=j;

  /* Initialize the stack */
  stack = (int *)malloc(NSTACK*sizeof(int));
  stack_idx=0;
  stack[stack_idx++]=0; /* initialize the left pointer to zero */
  stack[stack_idx++]=N-1; /* initialize the right pointer to N-1 */
  
  while(stack_idx>0) /* loop until done */
    {
      /* Pop the stack and begin a new round */
      r = stack[--stack_idx];
      l = stack[--stack_idx];

      /* If the number of remaining elements is above the threshold, 
	 we do a quicksort */
      if(r-l>=THRESH)
	{
	  /*******************************************************/
	  /* Median-of-three modification */
	  /*******************************************************/

	  /* Choose the median of the left (l), center (k), and right (r)
	     elements as the partitioning element */
	  k=(l+r)>>1; /* find the center element */
	  SWAP(idx[k],idx[r]);
	  if(A[idx[l]] > A[idx[r]]) SWAP(idx[l],idx[r]);
	  if(A[idx[l+1]] > A[idx[r]]) SWAP(idx[l+1],idx[r]);
	  if(A[idx[l]] > A[idx[l+1]]) SWAP(idx[l],idx[l+1]);

	  /*******************************************************/
	  /* Partition the subfile */
	  /*******************************************************/

	  i=l+1; /* Position of partitioning element */
	  j=r;
	  vidx=idx[i]; /* Index of partitioning element */
	  v=A[vidx]; /* Value of partitioning element */
	  
	  /* While the pointers haven't yet crossed */
	  while(i <= j)
	    {
	      while (A[idx[++i]] < v); /* Scan up from l+1 to find element > a */
	      while (A[idx[--j]] > v); /* Scan down from r to find element < a */
	      if(i <= j)
		SWAP(idx[i],idx[j]); /* Pointer didn't cross, exchange elements */
	    }
	  /* At the end, i=j+1 */

	  /* Insert partitioning element */
	  idx[l+1]=idx[j]; 
	  idx[j]=vidx;

	  /*******************************************************/
	  /* Push the pointers to the larger subfile on the stack,
	     process the smaller subfile now */
	  /*******************************************************/

	  if(stack_idx+4>NSTACK)
	    printf("Imminent stack overflow!\n");

	  /* Push the pointers for the right subfile on the stack,
	     sort the left subfile immediately */
	  if(r-i+1 >= j-l)
	    {
	      stack[stack_idx++] = i;
	      stack[stack_idx++] = r;
	      stack[stack_idx++] = l;
	      stack[stack_idx++] = j-1;
	    }
	  /* Push the pointers for the left subfile on the stack,
	     sort the right subfile immediately */
	  else
	    {
	      stack[stack_idx++] = l;
	      stack[stack_idx++] = j-1;
	      stack[stack_idx++] = i;
	      stack[stack_idx++] = r;
	    }
	}

      /*******************************************************/
      /* If the number of remaining elements is below the threshold,
	 do an insertion sort */
      /*******************************************************/
      else
	{
	  /* For each element from l to r */
	  for(j=l+1;j<=r;j++)
	    {
	      /* Pick out the current element v */
	      vidx=idx[j];
	      v=A[vidx];
	      
	      /* Find the place to insert v */
	      /* Start from the element to the left of v and proceed leftward */
	      for (i=j-1;i>=l;i--)
		{
		  /* If the current element is less than or equal to v,
		     then v is in the proper place */
		  if (A[idx[i]]<=v)
		    break; 
		  
		  /* Otherwise, shift all of the elements over and keep looking */
		  idx[i+1]=idx[i];
		}
      
	      /* Once we've found the proper place for v, set its index */ 
	      /* Since i is the index of the element to the left of v, i+1 is the index of v */
	      /* If we complete the loop, then v goes at the front of the list */
	      /* At that point, i=l-1, which makes i+1=l */
	      idx[i+1]=vidx;
	    }
	}
    }

  /* Sort the list according to the indices */
  for (j=0;j<N;j++)
    C[j]=A[idx[j]];

  free(stack);
}

void SortDouble(int N, double *A, double *C, int *idx)
{
  int u,n,j,U;
  int *run_length,*run_start,*run_idx,*run_sort;
  int run_flag;
  int *idx_temp;

  run_length = (int *)calloc(N,sizeof(int));
  run_start = (int *)malloc(N*sizeof(int));

  /* First do an unstable sort */
  SortUnstableDouble(N,A,C,idx);

  /* Next, scan the sorted list for runs of equal values */
  u=0;
  run_flag=0;
  for(n=1;n<N;n++)
    /* If this element is equal to the previous, we are in a run */
    if(C[n]==C[n-1])

      /* If we are in the middle of a run */
      if(run_flag)
	run_length[u-1]++; /* Increment the length of the run */

      /* Else, a run is just beginning */
      else
	{
	  /* Get the starting point of the run */
	  u++;
	  run_start[u-1] = n-1;
	  run_length[u-1] = 2;
	  run_flag=1;
	}
    /* Otherwise there is no run or the run has just ended */
    else
      run_flag=0;

  U=u;

  idx_temp = (int *)malloc(N*sizeof(int));
  memcpy(idx_temp,idx,N*sizeof(int));

  /* Sort the indices in each run */
  for(u=0;u<U;u++)
    {
      run_idx = (int *)malloc(run_length[u]*sizeof(int));
      run_sort = (int *)malloc(run_length[u]*sizeof(int));

      SortInt(run_length[u],&idx[run_start[u]],run_sort,run_idx);
 
     /* Reorder the indices */
      for(j=0;j<run_length[u];j++)
	idx[run_start[u]+j] = idx_temp[run_start[u]+run_idx[j]]; 
    
      free(run_idx);
      free(run_sort);
    }
  
  free(idx_temp);
  free(run_start);
  free(run_length);
}

void insertion_sort(double *A,int *idx,int l,int r)
{
  int i,j;
  int vidx;
  double v;

  /* For each element from l to r */
  for(j=l+1;j<=r;j++)
    {
      /* Pick out the current element v */
      vidx=idx[j];
      v=A[vidx];
      
      /* Find the place to insert v */
      /* Start from the element to the left of v and proceed leftward */
      for (i=j-1;i>=l;i--)
	{
	  /* If the current element is less than or equal to v,
	     then v is in the proper place */
	  if (A[idx[i]]<=v)
	    break; 
	  
	  /* Otherwise, shift all of the elements over and keep looking */
	  idx[i+1]=idx[i];
	}
      
      /* Once we've found the proper place for v, set its index */ 
      /* Since i is the index of the element to the left of v, i+1 is the index of v */
      /* If we complete the loop, then v goes at the front of the list */
      /* At that point, i=l-1, which makes i+1=l */
      idx[i+1]=vidx;
    }
}
