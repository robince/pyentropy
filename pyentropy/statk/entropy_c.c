/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

#include "toolkit_c.h"

/* Note: specific variance functions must be listed first with the entropy name followed by an underscore */
char ent_est_meth_list[ENT_EST_METHS][MAXCHARS] = {"plugin", "tpmc", "jack", "ma", "bub", "chaoshen", "ww", "nsb"};
char var_est_meth_list[GEN_VAR_EST_METHS+SPEC_VAR_EST_METHS][MAXCHARS] = {"nsb_var", "jack", "boot"};


/**
 * @brief Performs entropy and variance estimation
 * @param M Number of hist1d structures
 * @param in Pointer to hist1d structures
 * @param opts Pointer to an entropy options structure
 * @return Return status
 * Loops over M hist1d structures, using the specified options to calculate
 * entropy and variance estimates for each structure.
 */

/*int Entropy1DComp(int M,struct hist1d *in,struct options_entropy *opts)*/
/*{*/
	/*int m, e, v; [> loop variables <]*/
	/*int e_status, v_status, status = EXIT_SUCCESS; [> function status variables <]*/

	/* function pointers */
	/*int (*entropy_fun[ENT_EST_METHS+1])(struct hist1d *,struct options_entropy *,struct estimate *);*/
	/*int (*specific_variance_fun[SPEC_VAR_EST_METHS+1])(struct hist1d *,struct options_entropy *,struct nv_pair *);*/
	/*int (*general_variance_fun[GEN_VAR_EST_METHS])(struct hist1d *, int (*entropy_fun)(struct hist1d *, struct options_entropy *, struct estimate *), struct options_entropy *, struct nv_pair *);*/

	/*entropy_fun[0] = entropy_null;*/
	/*entropy_fun[1] = entropy_plugin;*/
	/*entropy_fun[2] = entropy_tpmc;*/
	/*entropy_fun[3] = entropy_jack;*/
	/*entropy_fun[4] = entropy_ma;*/
	/*entropy_fun[5] = entropy_bub;*/
	/*entropy_fun[6] = entropy_chaoshen;*/
	/*entropy_fun[7] = entropy_ww;*/
	/*entropy_fun[8] = entropy_nsb;*/

	/*specific_variance_fun[0] = variance_null;*/
	/*specific_variance_fun[1] = variance_nsb;*/

	/*general_variance_fun[0] = variance_jack;*/
	/*general_variance_fun[1] = variance_boot;*/

	/* loop over number of hist1d structures */
	/*for(m=0; m<M; m++)*/
	/*{*/
		/* loop over number of requested entropy methods */
		/*for(e=0; e<opts->E; e++)*/
		/*{*/
/*#ifdef DEBUG*/
			/*printf("\n(*opts).ent_est_meth[%d]=%d\n",e,(*opts).ent_est_meth[e]);*/
/*#endif*/
			/*if((*opts).ent_est_meth[e]>0) */
			/*{*/
/*#ifdef DEBUG*/
				/*printf("Applying entropy method ent_est_meth_list[%d]=\"%s\".\n",(*opts).ent_est_meth[e]-1,ent_est_meth_list[(*opts).ent_est_meth[e]-1]);*/
/*#endif*/
				/*e_status = entropy_fun[(*opts).ent_est_meth[e]](&in[m],opts,&(in[m].entropy[e]));*/
				/*if(e_status!=EXIT_SUCCESS)*/
					/*status = e_status;*/

				/* loop over number of variance methods */
				/*for(v=0; v<opts->V[e]; v++)*/
				/*{*/
/*#ifdef DEBUG*/
					/*printf("(*opts).var_est_meth[%d][%d]=%d\n",e,v,(*opts).var_est_meth[e][v]);*/
/*#endif*/
					/*if((*opts).var_est_meth[e][v]<=SPEC_VAR_EST_METHS)*/
					/*{*/
/*#ifdef DEBUG*/
						/*if((*opts).var_est_meth[e][v]==0)*/
							/*printf("Unrecognized specific variance method.\n");*/
						/*else if((*opts).ent_est_meth[e]==0)*/
							/*printf("Unrecognized entropy method.\n");*/
						/*else*/
							/*printf("Applying specific method var_est_meth_list[%d]=\"%s\".\n",(*opts).var_est_meth[e][v]-1,var_est_meth_list[(*opts).var_est_meth[e][v]-1]);*/
/*#endif*/
						/*v_status = specific_variance_fun[(*opts).var_est_meth[e][v]](&in[m],opts,&(in[m].entropy[e].ve[v]));*/
						/*if(v_status!=EXIT_SUCCESS)*/
							/*status = v_status;*/
					/*}*/
					/*else*/
					/*{*/
/*#ifdef DEBUG*/
						/*if((*opts).var_est_meth[e][v]==0)*/
							/*printf("Unrecognized general variance method.\n");*/
						/*else if((*opts).ent_est_meth[e]==0)*/
							/*printf("Unrecognized entropy method.\n");*/
						/*else*/
							/*printf("Applying general method var_est_meth_list[%d]=\"%s\" using entropy method ent_est_meth_list[%d]=\"%s\".\n",(*opts).var_est_meth[e][v]-1,var_est_meth_list[(*opts).var_est_meth[e][v]-1],(*opts).ent_est_meth[e]-1,ent_est_meth_list[(*opts).ent_est_meth[e]-1]);*/
/*#endif*/
						/*v_status = general_variance_fun[(*opts).var_est_meth[e][v]-SPEC_VAR_EST_METHS-1](&in[m],entropy_fun[(*opts).ent_est_meth[e]],opts,&(in[m].entropy[e].ve[v]));*/
						/*if(v_status!=EXIT_SUCCESS)*/
							/*status = v_status;*/
					/*}*/
				/*}*/
			/*}*/
		/*}*/
	/*}*/

	/*return status;*/
/*}*/

/**
 * @brief Returns the "plugin" or naive entropy estimate
 * @param in Pointer to a hist1d structure
 * Uses a hist1d structure to calculate an estimate of entropy via the naive
 * or "plugin" method.
 */
double EntropyPlugin(struct hist1d *in)
{
	double H;
	double N,*cnt;
	int i,m;

	N = (double) (*in).P;
	m = (*in).C;
	cnt = (*in).wordcnt;

	H=0;
	for(i=0; i<m; i++)
		H += XLOG2X(cnt[i]/N);
	
	return H;
}

/**
 * @brief Return the empirical maximum possible words.
 * @param in Pointer to a hist1d structure
 * @param minAgainstP Flag to indicate whether the result should be minimized against hist1d.P
 * @return Maximum number of possible words
 * Uses a hist1d structure to derive an empirical estimate of the maximum
 * number of possible words, which is the largest value in hist1d.wordlist
 * (biggest letter) plus 1, raised to the power hist1d.N (number of letters in
 * a word). If minAgainstP evaluates to true, then this result is minimized
 * against hist1d.P.
 */
double max_possible_words(struct hist1d *in,int minAgainstP)
{
	int max_bin,c,n;
	double m;

	max_bin = 0;
	for(c=0; c<(*in).C; c++)
		for(n=0; n<(*in).N; n++)
			if((*in).wordlist[c][n]>max_bin)
				max_bin = (*in).wordlist[c][n];
	max_bin++;
	if(minAgainstP)
		m = MIN(pow((double)max_bin,(double)(*in).N),(double)(*in).P);
	else
		m = pow((double)max_bin,(double)(*in).N);

#ifdef DEBUG
	printf("(*in).P=%d (*in).C=%d max_bin=%d (*in).N=%d m=%f\n",(*in).P,(*in).C,max_bin,(*in).N,m);
#endif

	return m;
}

/**
 * Place holder for an unrecognized entropy estimation technique
 */
int entropy_null(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy)
{
	printf("STAToolkit:entropy1d:unrecognizedOption: Unrecognized entropy estimation technique\n");
	
	return EXIT_FAILURE;
}

/**
 * Place holder for an unrecognized variance estimation technique
 */
int variance_null(struct hist1d *in,struct options_entropy *opts,struct nv_pair *variance)
{
	printf("STAToolkit:entropy1d:unrecognizedOption: Unrecognized variance estimation technique\n");
	
	return EXIT_FAILURE;
}

/**
 * This function is used by jni toolkit.
 * @param index Index of the element in ent_est_meth_list
 * @return Element at the specified index location
 */
char *getEstMethList(int index)
{
	return ent_est_meth_list[index];
}

