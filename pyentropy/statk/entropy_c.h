/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "statk_license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

#define ENT_EST_METHS 8
#define SPEC_VAR_EST_METHS 1
#define GEN_VAR_EST_METHS 2
#define DEFAULT_UNOCCUPIED_BINS_STRATEGY -1
#define DEFAULT_TPMC_POSSIBLE_WORDS_STRATEGY 0
#define DEFAULT_BUB_POSSIBLE_WORDS_STRATEGY 2
#define DEFAULT_BUB_LAMBDA_0 0
#define DEFAULT_BUB_K 11
#define DEFAULT_BUB_COMPAT 0
#define DEFAULT_WW_POSSIBLE_WORDS_STRATEGY 0
#define DEFAULT_WW_BETA 1
#define DEFAULT_BOOT_NUM_SAMPLES 100
#define DEFAULT_BOOT_RANDOM_SEED 1
#define DEFAULT_POSSIBLE_WORDS -1.0
#define DEFAULT_NSB_PRECISION 1e-6

struct options_entropy{
  int **var_est_meth; int var_est_meth_flag; int *V; /**< Number of variance estimation methods on each entropy estimation method. */
  int *ent_est_meth; int E; /**< Number of entropy estimation methods. */
  char **ent_est_meth_name; /**< used by jni toolkit */
  int useall; int useall_flag;
  int tpmc_possible_words_strategy; int tpmc_possible_words_strategy_flag;
  int bub_possible_words_strategy; int bub_possible_words_strategy_flag;
  double bub_lambda_0; int bub_lambda_0_flag;
  int bub_K; int bub_K_flag;
  int bub_compat; int bub_compat_flag;
  int ww_possible_words_strategy; int ww_possible_words_strategy_flag;
  double ww_beta; int ww_beta_flag;
  int boot_num_samples; int boot_num_samples_flag;
  int boot_random_seed; int boot_random_seed_flag;
  double possible_words; int possible_words_flag; /**< Possible words in BUB, NSB, TPMC, and WW calculations. */
  double nsb_precision; int nsb_precision_flag; /**< Relative precision on NSB calculations. */
};

struct message{
  int i, j, k; /* Number of status, warning, and error messages */
  char **status;
  char **warnings;
  char **errors;
};

struct nv_pair{
  char name[MAXCHARS];
  double value;
};

struct estimate{
  char name[MAXCHARS];
  double value;
  struct message *messages;
  int E; /* Number of extras */
  struct nv_pair *extras;
  int V; /* Number of variance estimation methods */
  struct nv_pair *ve;
};

struct hist1d{
  int P;       /* Number of words used to generate the estimate */
  int C;       /* Number of unique words */
  int N;       /* Number of subwords */
  int **wordlist; /* List of words that appear (C long) */
  double *wordcnt; /* Number of times each word occurs (C long) */
  struct estimate *entropy;
};

extern int Entropy1DComp(int M,struct hist1d *in,struct options_entropy *opts);
extern double EntropyPlugin(struct hist1d *in);
extern double max_possible_words(struct hist1d *in,int minAgainstP);

extern int entropy_null(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int variance_null(struct hist1d *in,struct options_entropy *opts,struct nv_pair *variance);

extern int entropy_plugin(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_tpmc(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_jack(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_ma(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_bub(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_chaoshen(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_ww(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);
extern int entropy_nsb(struct hist1d *in,struct options_entropy *opts,struct estimate *entropy);

extern int variance_jack(struct hist1d *in,int (*entropy_fun)(),struct options_entropy *opts,struct nv_pair *variance);
extern int variance_boot(struct hist1d *in,int (*entropy_fun)(),struct options_entropy *opts,struct nv_pair *variance);
extern int variance_nsb(struct hist1d *in,struct options_entropy *opts,struct nv_pair *variance);

extern char* getEstMethList(); /**< used by jni toolkit */
