/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */
/* Sorting functions */
extern int IsSortedDouble(int M,double *in);
extern void SortUnstableInt(int N, int *B, int *C, int *idx);
extern void SortInt(int N, int *A, int *C, int *idx);
extern int CountSortedInt(int M,int *sorted,int *sort_idx,int *uni_list,int *uni_i,int *uni_j,int *cnt);
extern int UniqueInt(int M,int *list,int *uni_list,int *uni_i,int *uni_j,int *cnt);

extern void SortRowsInt(int P,int N,int **list,int **sorted,int *p_list);
extern int CountSortedRowsInt(int M,int N,int **sorted,int *sort_idx,int **uni_list,int *uni_i,int *uni_j,int *cnt);
extern int UniqueRowsInt(int M,int N,int **list,int **uni_list,int *uni_i,int *uni_j,int *cnt);

extern void SortUnstableDouble(int N, double *A, double *C, int *idx);
extern void SortDouble(int N, double *A, double *C, int *idx);
extern int CountSortedDouble(int M,double *sorted,int *sort_idx,double *uni_list,int *uni_i,int *uni_j,int *cnt);
extern int UniqueDouble(int M,double *list,double *uni_list,int *uni_i,int *uni_j,int *cnt);

extern void SortRowsDouble(int P,int N,double **list,double **sorted,int *p_list);
extern int CountSortedRowsDouble(int M,int N,double **sorted,int *sort_idx,double **uni_list,int *uni_i,int *uni_j,int *cnt);
extern int UniqueRowsDouble(int M,int N,double **list,double **uni_list,int *uni_i,int *uni_j,int *cnt);
