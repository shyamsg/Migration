/*
 *  migrate.h
 *  Migrate
 *
 *  Created by Shyam Gopalakrishnan on 6/1/12.
 *  Copyright 2012 __Uchicago__. All rights reserved.
 *
 */

#ifndef __MIGRATE_H__
#define __MIGRATE_H__

#include "cmath"
#include "iostream"
#include "cstring"

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "nlopt.h"

#include "ext/hash_map"
#include "vector"
#include "iterator"
#include "algorithm"

using namespace std;
using namespace __gnu_cxx;

#define EPS_GRAD       1e-12
#define NRESTARTS      40
#define FTOL           1e-15
#define MAXEVAL        500000
#define LOW_COAL_RATE  1e-6

typedef unsigned int uint;
typedef struct {
  double t;
  gsl_vector_view obs_coal_rates;
  gsl_matrix * P0;
  vector<vector<int> > * popdict;
  uint count;
  bool logVal;
} cfnm_data;

void gsl_matrix_print(gsl_matrix * M);
int invert_matrix(const gsl_matrix *A, gsl_matrix *Ai)  throw (const char *);
gsl_matrix * comp_pw_coal_cont(gsl_matrix * m, gsl_vector * Ne_inv);
gsl_matrix * expM(gsl_matrix * Q);
gsl_matrix * conv_scrambling_matrix(gsl_matrix * P);
vector<vector <int> > pop_to_col(vector<vector<int> > &popD, int nd_old);
gsl_matrix * converge_pops(vector<vector<int> > & popDict, gsl_matrix * scramb);
gsl_matrix * compute_pw_coal_rates(vector<vector<double> > & Nes, \
				   vector<vector<double> > & ms,	\
				   vector<double> ts,			\
				   vector<vector<vector<int> > > & popmaps);
vector<vector<int> > construct_poparr(vector<vector<int> > popdict);
vector<vector<int> > find_pop_merges(gsl_vector * Ninv, vector<double> mtemp,\
				     double t, gsl_matrix * P0,		\
				     double merge_threshold, bool useMigration=false);
gsl_vector * average_coal_rates(gsl_vector_view origrates, \
				  vector<vector<int> > & popdict);
vector<vector<int> > make_merged_pd(vector<vector<vector<int> > > & pdlist);
double compute_dist_and_grad(uint n, const double * x, double * grad, void * data);
double compute_2norm_mig(cfnm_data * d, gsl_matrix * m, gsl_vector * Ne_inv);
vector< vector<double> > comp_params(gsl_matrix * obs_rates, vector <double> t, \
				     vector<vector<vector<int > > > &pdlist, \
				     bool logVal, double merge_threshold, \
				     bool useMigration=false);
#endif

