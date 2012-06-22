#include "migrate.h"

#define DIM_MISMATCH 127
#define MAXPOPS 50

using namespace std;
using namespace __gnu_cxx;

/* General Method to print matrices*/
void gsl_matrix_print(gsl_matrix * M)
{
  cout << "[ ";
  for (uint i=0; i<M->size1; i++) {
    cout << "[ ";
    copy(M->data+(i*M->size2), M->data+((i+1)*M->size2), ostream_iterator<double>(cout, " "));
    cout << " ]" << endl;
  }
  cout << "]" << endl;
}

/* matrix inversion routing using LU-factorization */
int invert_matrix(const gsl_matrix *A, gsl_matrix *Ai) throw (const char *)
{
  /* check matrix sizes */
  if(! A->size1 == A->size2 && Ai->size1 == Ai->size2 && A->size1 == Ai->size1)
    return (-1);
  /* create a permutation matrix for the LU-factorization */
  gsl_permutation * P = gsl_permutation_alloc(A->size1);
  int s; 			/* a dummy variable */
  /* allocate a buffer matrix */
  gsl_matrix *AA = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(AA, A);	/* copy A to AA */
  /* perform LU-factorization */
  gsl_linalg_LU_decomp(AA, P, &s);
  double det = gsl_linalg_LU_det(AA, s);
  if (fabs(det) < 1e-35) {
    throw "Determinant is low.";
  }
  /* backsubstitute to get the inverse */
  gsl_linalg_LU_invert(AA, P, Ai);
  gsl_matrix_free(AA);
  gsl_permutation_free(P);
  return 0;
}

/**************************************************************
Function to evaluate the coal. intensities for the given 
migration and effective population sizes, in the time window 
of t generations. Here we use compute an infinitesimal rate 
matrix Q from the migration rates and population sizes and use 
the exponentiation of this rate matrix to compute the transition
probabilities.
***************************************************************/
gsl_matrix * comp_pw_coal_cont(gsl_matrix * m, gsl_vector * Ne_inv)
{
  unsigned int numdemes = Ne_inv->size;
  if (m->size1 != numdemes || m->size2 != numdemes) {
    printf("Error number %d: Migration matrix and population size vector indicate different number of demes.", DIM_MISMATCH);
    exit(DIM_MISMATCH);
  }
  int nr = (numdemes*(numdemes+1))/2 + 1;
  gsl_matrix * Q = gsl_matrix_calloc(nr,nr);
  int rnum = -1;
  gsl_matrix_int * demePairs = gsl_matrix_int_alloc(numdemes, numdemes);
  for (unsigned int p1=0; p1<numdemes; p1++){
    for (unsigned int p2=p1; p2<numdemes; p2++){
      rnum++;
      gsl_matrix_int_set(demePairs, p1, p2, rnum);
      gsl_matrix_int_set(demePairs, p2, p1, rnum);
    }
  }
  for (unsigned int p1=0; p1<numdemes; p1++) {
    for (unsigned int p2=p1; p2<numdemes; p2++) {
      if (p1 == p2) {
	gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1,p2), nr-1, 0.5*gsl_vector_get(Ne_inv, p1));
      }
      for (unsigned int p3=0; p3<numdemes; p3++) {
	for (unsigned int p4=p3; p4<numdemes; p4++) {
	  if (p1 == p3 && p2 == p4) {
	    continue;
	  }
	  if (p1 == p2) { //both lines in same deme
	    if (p3 == p4) {
	      continue; // needs 2 migrations
	    } else if (p3 == p1) {
	      gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			     gsl_matrix_int_get(demePairs, p3, p4),		\
			     2*gsl_matrix_get(m, p2, p4));
	    } else if (p4 == p2) {
	      gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			     gsl_matrix_int_get(demePairs, p3, p4),		\
			     2*gsl_matrix_get(m, p1, p3));
	    }
	  } else { // two lines in 2 demes
	    if (p3 == p1 || p3 == p2) {
	      if (p1 == p3){
		gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			       gsl_matrix_int_get(demePairs, p3, p4),	\
			       gsl_matrix_get(m ,p2, p4));
	      } else {
		gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			       gsl_matrix_int_get(demePairs, p3, p4),	\
			       gsl_matrix_get(m ,p1, p4));
	      }
	    } else if (p4 == p2 || p4 == p1) {
	      if (p1 == p4){
		gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			       gsl_matrix_int_get(demePairs, p3, p4),	\
			       gsl_matrix_get(m ,p2, p3));
	      } else {
		gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),	\
			       gsl_matrix_int_get(demePairs, p3, p4),	\
			       gsl_matrix_get(m ,p1, p3));
	      }
	    }
	  }
	}
      }
    }
  }    
  for (unsigned int p1=0; p1<numdemes; p1++) {
    for (unsigned int p2=p1; p2<numdemes; p2++) {
      double temp = -gsl_matrix_get(Q, gsl_matrix_int_get(demePairs, p1, p2), \
				    nr-1);
      for (unsigned int p3=0; p3<numdemes; p3++) {
	for (unsigned int p4=p3; p4<numdemes; p4++) {
	  temp -= gsl_matrix_get(Q, gsl_matrix_int_get(demePairs, p1, p2), \
				 gsl_matrix_int_get(demePairs, p3, p4));
	}
      }
      gsl_matrix_set(Q, gsl_matrix_int_get(demePairs, p1, p2),		\
		     gsl_matrix_int_get(demePairs, p1, p2), temp);
    }
  }
  gsl_matrix_int_free(demePairs);
  return Q;
}

/***********************************************************
Computes the exponential of the matrix using eig decomposition.
***********************************************************/
gsl_matrix * expM(gsl_matrix * Q)
{
  gsl_vector_complex *eval = gsl_vector_complex_alloc(Q->size1);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(Q->size1, Q->size1); 
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(Q->size1); 
  gsl_eigen_nonsymmv(Q, eval, evec, w);
  gsl_eigen_nonsymmv_free(w); 
  gsl_matrix *U = gsl_matrix_calloc(Q->size1, Q->size1);
  for (unsigned int row=0; row<Q->size1; row++) {
    for (unsigned int col=0; col<Q->size1; col++) {
      gsl_complex z = gsl_matrix_complex_get(evec, row, col);
      gsl_matrix_set(U, row, col, GSL_REAL(z));
    }
  }
  gsl_matrix *L = gsl_matrix_calloc(Q->size1, Q->size1);
  for (unsigned int row=0; row<Q->size1; row++) {
    gsl_complex z = gsl_vector_complex_get(eval, row);
    gsl_matrix_set(L, row, row, exp(GSL_REAL(z)));
  }
  gsl_matrix *eQ = gsl_matrix_calloc(Q->size1, Q->size1);
  gsl_matrix *temp = gsl_matrix_calloc(Q->size1, Q->size1);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, U, L, 0.0, temp);
  gsl_matrix *Uinv = gsl_matrix_calloc(Q->size1, Q->size1);
  invert_matrix(U, Uinv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp, Uinv, 0.0, eQ);
  gsl_matrix_free(U);
  gsl_matrix_free(Uinv);
  gsl_matrix_free(L);
  gsl_matrix_free(temp);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(eval);
  return eQ;
}

/***********************************************************
This takes a Probability matrix and converts it to a scrambling matrix, to
get the probability of lines moving from the current demes to any config of
demes at the beginning of the time slice.
***********************************************************/
gsl_matrix * conv_scrambling_matrix(gsl_matrix * P) 
{
  int nr = P->size1;
  int nc = P->size2;
  gsl_matrix * Pnew = gsl_matrix_alloc(nr, nc);
  gsl_matrix_memcpy(Pnew, P);
  for (int row=0; row<nr-1; row++) {
    gsl_matrix_set(Pnew, row, nc-1, 0.0);
  }
  for (int row=0; row<nr; row++) {
    double sum = 0.0;
    for (int col=0; col<nc; col++) {
      sum += gsl_matrix_get(Pnew, row, col);
    }
    gsl_vector_view currRow = gsl_matrix_row(Pnew, row);
    gsl_vector_scale(&currRow.vector, 1.0/sum);
  }
  return Pnew;
}

/***********************************************************
Given the population indices, it tells us the correspoiding 
columns where lines come from any of the pops in the list.
***********************************************************/
vector<vector <int> > pop_to_col(vector<vector<int> > &popD, int nd_old)
{
  uint nd = popD.size();
  vector<vector<int> > ptc;
  for (uint i=0; i<nd; i++) {
    for (uint j=i; j<nd; j++) { 
      vector<int> cols;
      for (unsigned int k=0; k<popD[i].size(); k++){
	for (unsigned int l=0; l<popD[j].size(); l++) {
	  if (popD[i][k] > popD[j][l]) continue;
	  int tt = ((popD[i][k]*(2*nd_old-popD[i][k]+1))/2 + popD[j][l] - popD[i][k]);
	  cols.push_back(tt);
	}
      }
      ptc.push_back(cols);
    }
  }
  return ptc;
}

/***********************************************************
This function changes the scrambling matrix appropriately upon
the merging of populations.
***********************************************************/
gsl_matrix * converge_pops(vector<vector<int> > & popDict, gsl_matrix * scramb)
{
  uint npop = scramb->size2;
  npop = int(sqrt(8*npop - 7)/2);
  uint nr = scramb->size1;
  uint nc = popDict.size();
  nc = (nc*(nc+1))/2+1;
  gsl_matrix * temp = gsl_matrix_calloc(nr,nc);
  for (uint row=0; row<nr; row++) {
    gsl_matrix_set(temp, row, nc-1, gsl_matrix_get(scramb, row, scramb->size2-1));
  }
  vector<vector<int> > ptc = pop_to_col(popDict, npop);
  for (uint col=0; col<nc-1; col++) {
    for (uint row=0; row<nr; row++) {
      for (uint scol=0; scol<ptc[col].size(); scol++) {
	gsl_matrix_set(temp, row, col, (gsl_matrix_get(temp, row, col) + \
					gsl_matrix_get(scramb, row, ptc[col][scol])) \
		       );
      }
    }
  }
  return temp;
}

/************************************************************
Given a list of population sizes and the migration matrix,
one for each time slice, compute the expected coalescent rates
for all possible pairs of lines. Uses the comp_pw_coal_cont 
function to do this for each time slice.
The lists ms and Nes must be ordered from most recent to most 
ancient. Here ts is the list of time slice lengths - in the same
order. popmaps is a list which tells which pops merge back 
in time at the beginning of this time slice. Each entry
of popmaps is a list of arrays or tuples with all pops that 
merge into one
************************************************************/
gsl_matrix * compute_pw_coal_rates(vector<vector<double> > & Nes, \
                                   vector<vector<double> > & ms, \
                                   vector<double> ts, \
                                   vector<vector<vector<int> > > & popmaps)
{
  /* 
     if the number of keys in the popmap is not equal to 
     the number of pops (ms and Ns), do the converge thing
     to bring pops together --> change at the beginning of each
     iteration, so the pop map tells the method what to do 
     --->BEFORE<--- computing rates in each time slice
  */
    
  uint numslices = ts.size();
  uint numpops = Nes[0].size();
  uint nr = (numpops*(numpops+1))/2 + 1;
  gsl_matrix * exp_rates = gsl_matrix_calloc(nr-1, numslices);
  gsl_matrix * P0 = gsl_matrix_calloc(nr, nr);
  for (uint rc=0; rc<nr; rc++) {
    gsl_matrix_set(P0, rc, rc, 1.0);
  }
  for (uint i=0; i<numslices; i++) {
    //    cout << i << " " << numslices << " " << Nes[i].size() << endl;
    gsl_vector * Ne_inv = gsl_vector_calloc(Nes[i].size());
    for (uint el=0; el<Ne_inv->size; el++) {
      Ne_inv->data[el] = 1.0/Nes[i][el];
    }
    vector<double> mtemp = ms[i];
    numpops = Ne_inv->size;
    gsl_matrix * m = gsl_matrix_calloc(numpops,numpops);
    uint cnt = 0;
    // if nothing to do, popmaps[i] must be empty
    if (popmaps[i].size() == numpops){
      P0 = converge_pops(popmaps[i], P0);
    } else if (popmaps[i].size() > 0) {
      cout << "PopMaps say: " << popmaps[i].size() << \
	" populations whereas population sizes say: " << numpops << \
	"populations. Too much uncertainty. Does not compute :)." << endl;
      throw "Population map and population size parameters do not match";
    }
    for (uint ii=0; ii<numpops; ii++) {
      for (uint jj=ii+1; jj<numpops; jj++) {
	m->data[ii*m->size2+jj] = m->data[jj*m->size2+ii] = mtemp[cnt];
	cnt += 1;
      }
    }
    gsl_matrix * Q = comp_pw_coal_cont(m, Ne_inv);
    gsl_vector_free(Ne_inv);
    gsl_matrix_free(m);
    gsl_matrix_scale(Q, ts[i]);
    //    cout << "Q:"; gsl_matrix_print(Q);
    gsl_matrix * eQ = expM(Q);
    //    cout << "eQ:"; gsl_matrix_print(eQ);
    gsl_matrix_free(Q);
    gsl_matrix * P = gsl_matrix_alloc(P0->size1, eQ->size2);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, P0, eQ, 0.0, P);
    //    cout << "P:"; gsl_matrix_print(P);
    //    cout << "P0:"; gsl_matrix_print(P0);
    for (uint row=0; row<exp_rates->size1; row++) {
      exp_rates->data[row*exp_rates->size2+i] = P->data[row*P->size2+P->size2-1];
      if (exp_rates->data[row*exp_rates->size2+i] < 0) {
	exp_rates->data[row*exp_rates->size2+i] = 0.0;
      }
    }
    gsl_matrix * temp = conv_scrambling_matrix(eQ);
    gsl_matrix * temp2 = gsl_matrix_alloc(P0->size1, temp->size2);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, P0, temp, 0.0, temp2);
    gsl_matrix_free(P0);
    P0 = gsl_matrix_alloc(temp2->size1, temp->size2);
    gsl_matrix_memcpy(P0, temp2);
    gsl_matrix_free(temp);
    gsl_matrix_free(P);
    gsl_matrix_free(temp2);
  }
  return exp_rates;
}

/***********************************************************
Convert the array from comp_N_m into the desired format for 
converge pops.
***********************************************************/
vector<vector<int> > construct_poparr(vector<vector<int> > popdict)
{
  hash_map<int, int> popToIndex;
  uint nd = popdict.size();
  uint currIndex = 0;
  for (uint ii=0; ii<nd; ii++) {
    if (popdict[ii].size() == 0) {
      popToIndex[ii] = currIndex;
      currIndex += 1;
    } else if (popToIndex.find(ii) == popToIndex.end()) {
      popToIndex[ii] = currIndex;
      for (uint jj=0; jj<popdict[ii].size(); jj++) {
	  popToIndex[popdict[ii][jj]] = popToIndex[ii];
      }
      currIndex += 1;
    } else if (popToIndex.find(ii) != popToIndex.end()) {
      for (uint jj=0; jj<popdict[ii].size(); jj++) {
	popToIndex[popdict[ii][jj]] = popToIndex[ii];
      }
    }
  }

  vector<vector<int> > popmap;
  for (uint ii=0; ii<currIndex; ii++) {
    vector<int> temp;
    for (hash_map<int, int>::iterator hit=popToIndex.begin();	\
	 hit!=popToIndex.end(); hit++) {
      if (hit->second == int(ii)) {
	temp.push_back(hit->first);
      }
    }
    popmap.push_back(temp);
  }
  return popmap;
}

/***********************************************************
This function takes the optimal paramters found and figures 
out if the populations need to be merged.

 ***********************************************************/
vector<vector<int> > find_pop_merges(gsl_vector * Ninv, vector<double> mtemp,\
                                     double t, gsl_matrix * P0, \
                                     double merge_threshold, bool useMigration) 
{
  uint numdemes = Ninv->size;
  vector<vector<int> > popdict = vector<vector<int> >(numdemes);
  if (!useMigration) {
    cout << "Using coalescent rates." << endl;
    gsl_matrix * m = gsl_matrix_calloc(numdemes,numdemes);
    uint cnt = 0;
    for (uint ii=0; ii<numdemes; ii++) {
      for (uint jj=ii+1; jj<numdemes; jj++) {
	m->data[ii*m->size2+jj] = m->data[jj*m->size2+ii] = mtemp[cnt];
	cnt += 1;
      }
    }
    gsl_matrix * Q = comp_pw_coal_cont(m, Ninv);
    gsl_matrix_scale(Q, t);
    gsl_matrix *eQ = expM(Q);
    gsl_matrix_free(m);
    gsl_matrix_free(Q);
    //gsl_vector_view P = gsl_matrix_subcolumn(eQ, eQ->size2-1, 0, eQ->size1-1);
    gsl_vector_view Prow = gsl_matrix_column(eQ, eQ->size2-1);
    gsl_vector_view P = gsl_vector_subvector(&Prow.vector, 0, eQ->size1-1);
    for (uint i=0; i<numdemes; i++) {
      for (uint j=i+1; j<numdemes; j++) {
	gsl_vector * dists = gsl_vector_alloc(3);;
	dists->data[0] = P.vector.data[((2*numdemes-i+1)*i)/2];
	dists->data[1] = P.vector.data[((2*numdemes-i+1)*i)/2+(j-i)];
	dists->data[2] = P.vector.data[((2*numdemes-j+1)*j)/2];
	double meanRates = gsl_stats_mean(dists->data, 1, dists->size);
	//	medianRates = np.median(dists);
        double rangeRates = gsl_vector_max(dists) - gsl_vector_min(dists);
	//      print i,j, medianRates, meanRates
	if (rangeRates/meanRates < merge_threshold) {
	  //	  print 'Hello:', i, j, np.real(dists);
	  popdict[i].push_back(j);
	  popdict[j].push_back(i);
	}
	gsl_vector_free(dists);
      }
    }
    gsl_matrix_free(eQ);
  } else {
    //    cout << "mtemp " << mtemp.size() << endl;
    //    copy(mtemp.begin(), mtemp.end(), ostream_iterator<double>(cout, " "));
    //    cout << endl;
    uint cnt = 0;
    for (uint ii=0; ii<numdemes; ii++) {
      for (uint jj=ii+1; jj<numdemes; jj++) {
	if (mtemp[cnt] > merge_threshold) {
	  popdict[ii].push_back(jj);
	  popdict[jj].push_back(ii);
	}
	cnt += 1;
      }
    }
  }
  return construct_poparr(popdict);
}

/***********************************************************
This function changes the rates to average them
when the populations merge. popdictlist is a list
of population dictionaries. rates is the column vector
of coalescent rates.
When applying to observed rates, the popdictlist must just 
contain the latest list. While applying to the estimated 
rates from the current N and m, the whole list must be given.
************************************************************/
gsl_vector * average_coal_rates(gsl_vector_view origrates, \
                                vector<vector<int> > & popdict)
{
  gsl_vector *newRates;
  if (popdict.size() == 0) {
    newRates = gsl_vector_alloc(origrates.vector.size);
    gsl_vector_memcpy(newRates, &origrates.vector);
    return newRates;
  }
  uint numdemes = 0;
  for (uint i=0; i<popdict.size(); i++) {
    numdemes += popdict[i].size();
  }
  vector<vector<int> > ptc = pop_to_col(popdict, numdemes);
  newRates = gsl_vector_alloc(ptc.size());
  for (uint row=0; row<ptc.size(); row++) {
    double totRates = 0.0;
    for (uint row2=0; row2<ptc[row].size(); row2++) {
      totRates += gsl_vector_get(&origrates.vector, ptc[row][row2]);
    }
    totRates /= ptc[row].size();
    newRates->data[row] = totRates;
  }
  return newRates;
}

/************************************************************
Merges the sequence of population dictionaries to get one 
big one.
************************************************************/
vector<vector<int> > make_merged_pd(vector<vector<vector<int> > > & pdlist)
{
  vector<vector<int> > pdnew = vector<vector<int> >();
  if (pdlist.size() == 0) {
    return pdnew;
  } else if (pdlist.size() == 1) {
    pdnew.insert(pdnew.end(), pdlist[0].begin(), pdlist[0].end());
    return pdnew;
  }
  vector<vector<int> >::iterator vvit;
  vector<int>::iterator vit;
  pdnew.insert(pdnew.end(), pdlist[pdlist.size()-1].begin(), pdlist[pdlist.size()-1].end());
  for (int i=pdlist.size()-1; i>0; i--) {
    vector<vector<int> > pdprev = vector<vector<int> >(pdlist[i-1]);
    vector<vector<int> > pdtemp;
    for (vvit=pdnew.begin(); vvit<pdnew.end();vvit++) {
      vector<int> temp;  
      for (vit=vvit->begin(); vit<vvit->end(); vit++) {
        temp.insert(temp.end(), pdprev[*vit].begin(), pdprev[*vit].end());
      }
      sort(temp.begin(), temp.end());
      pdtemp.push_back(temp);
    }
    pdnew.clear();
    pdnew.insert(pdnew.end(), pdtemp.begin(), pdtemp.end());
  }
  return pdnew;
}
/***********************************************************
Function takes the migration matrix, the popsize inverses, 
the cfnm_data (auxillary data).
************************************************************/

double compute_2norm_mig(cfnm_data * d, gsl_matrix * m, gsl_vector * Ne_inv)
{
  gsl_matrix * Q = comp_pw_coal_cont(m, Ne_inv);
  gsl_matrix_scale(Q, d->t);
  gsl_matrix * Pcurr = expM(Q);
  gsl_matrix_free(Q);
  gsl_matrix * temp = gsl_matrix_alloc(d->P0->size1, Pcurr->size2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, d->P0, Pcurr, 0.0, temp);
  gsl_matrix_free(Pcurr);
  //  gsl_vector_view estRates = gsl_matrix_subcolumn(temp, temp->size2-1, 0, 
  //						  temp->size1-1);
  gsl_vector_view estRow = gsl_matrix_column(temp, temp->size2-1);
  gsl_vector_view estRates = gsl_vector_subvector(&estRow.vector, 0, temp->size1-1);
  gsl_matrix_free(temp);
  if (d->obs_coal_rates.vector.size != estRates.vector.size) {
    throw "Dimensions of observed and computed vectors do not match.";
  }
  //  cout << "testing" << endl;
  //  cout << d->popdict->size() << endl;
  gsl_vector * avobsrates = average_coal_rates(d->obs_coal_rates, *(d->popdict));
  gsl_vector * avestrates = average_coal_rates(estRates, *(d->popdict));
  //  cout << "t2" << endl;
  gsl_vector_sub(avobsrates, avestrates);
  gsl_vector_free(avestrates);
  double fnorm = gsl_blas_dnrm2(avobsrates);
  gsl_vector_free(avobsrates);
  return fnorm;
}

/***********************************************************
Function to compute the order 2 norm distance between the 
observed coalescent intensities and coalescent intensities 
computed using the given N and m (contained in x).
The input is given as x - a vector with N and m values
P0 is the initial scrambling matrix. Note that we work with 1/N and
m values, as preconditioning to make the problem more spherical.
***********************************************************/
double compute_dist_and_grad(unsigned int n, const double * x, double * grad, void * data)
{
  cfnm_data * d = (cfnm_data *) data;
  d->count++;
#ifdef DEBUG
  if(d->count%100000 == 0) cout << "Done with " << d->count << " evals." << endl;
#endif
  uint k = int((sqrt(1+8*n) - 1)/2.0);
  gsl_vector * Ne_inv = gsl_vector_alloc(k);
  for (uint i=0; i<k; i++) {
    Ne_inv->data[i] = x[i];
  }
  gsl_matrix * m = gsl_matrix_calloc(k,k);
  uint cnt = 0;
  // make the m matrix from ms
  for (uint i=0; i<k; i++) {
    for (uint j=i+1; j<k; j++) {
      m->data[i*k+j] = m->data[j*k+i] = x[k+cnt];
      cnt += 1;
    }
  }
  double fnorm = compute_2norm_mig(d, m, Ne_inv);
  //#    grad = grad_Frob_cdiff(x, t, obs_coal_rates, P0, popdict, epsilon=1e-5)
  //#    print grad
  //#    return tempo, grad
  if (grad) {
    for (uint nd=0; nd<k; nd++) {
      Ne_inv->data[nd] += EPS_GRAD;
      grad[nd] = (fnorm - compute_2norm_mig(d, m, Ne_inv))/EPS_GRAD;
      Ne_inv->data[nd] -= EPS_GRAD;
    }
    cnt = 0;
    for (uint p1=0; p1<k; p1++) {
      for (uint p2=p1+1; p2<k; p2++) {
	m->data[p1*k+p2] += EPS_GRAD; 
	m->data[p2*k+p1] += EPS_GRAD;
	grad[k+cnt] = (fnorm - compute_2norm_mig(d, m, Ne_inv))/EPS_GRAD;
	m->data[p1*k+p2] -= EPS_GRAD;
	m->data[p2*k+p1] -= EPS_GRAD;
      }
    }
  }
  gsl_matrix_free(m);
  gsl_vector_free(Ne_inv);
  return fnorm;
}

/***********************************************************
This function estimates the N and m parameters for the various time slices. The 
time slices are given in a vector form 't'. t specifies the length of the time 
slice, not the time from present to end of time slice (not cumulative but atomic).
Also, the obs_rates are given, for each time slice are given in columns of the 
obs_rates matrix. Both obs_rates and time slice lengths are given from present to 
past.
***********************************************************/
vector< vector<double> > comp_params(gsl_matrix * obs_rates, vector <double> t, \
                                     vector<vector<vector<int > > > &pdlist, \
                                     double merge_threshold, bool useMigration)
{
  /* SETUP FOR UNIFORM RANDOM GENERATION */
  const gsl_rng_type * T; 
  gsl_rng * r; 
  gsl_rng_env_setup(); 
  T = gsl_rng_default; 
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, time(NULL));
  /* END RANDGEN SETUP */

  double minf;
  uint numslices = t.size();
  uint nr = obs_rates->size1 + 1;
  uint numdemes = int((sqrt(8*nr - 7))/2);
  cfnm_data * d = (cfnm_data * ) malloc(sizeof(cfnm_data));
  d->P0 = gsl_matrix_alloc(nr, nr);
  vector<vector<int> > tempPopdict;
  tempPopdict.reserve(MAXPOPS);
  d->popdict = &tempPopdict;
  gsl_matrix_set_identity(d->P0);
  vector<vector<double> > xopts;
  pdlist.clear();
  nlopt_opt opt;
  for (uint ns=0; ns<numslices; ns++) {
    cout << endl << "starting slice number " << ns << endl;
    gsl_vector * Ninv;
    vector<double> mtemp;
    bool reestimate = 0;
    do {
      uint nparams = uint((numdemes*(numdemes+1))/2.0);
      // Set aug params for function
      d->t = t[ns];
      vector<vector<int> > pdmerged = make_merged_pd(pdlist);
      tempPopdict.clear();
      tempPopdict.insert(tempPopdict.end(), pdmerged.begin(), pdmerged.end());
      d->obs_coal_rates = gsl_matrix_column(obs_rates, ns);

      double lb[nparams];
      double ub[nparams];
      for (uint np=0; np<numdemes; np++) {
	lb[np] = 1e-15;
	ub[np] = 1e-1;
      }
      gsl_vector * temprates = average_coal_rates(d->obs_coal_rates, pdmerged);
      uint cnt = numdemes;
      for (uint ii=0; ii<numdemes; ii++) {
	for (uint jj=ii+1; jj<numdemes; jj++) {
	  lb[cnt] = 1e-15; //(temprates->data[((2*numdemes-ii+1)*ii)/2+(jj-ii)] < LOW_COAL_RATE) ? 0 : 1e-15;
	  ub[cnt] = (temprates->data[((2*numdemes-ii+1)*ii)/2+(jj-ii)] < LOW_COAL_RATE)\
	    ? 1e-15 : 1e-1;
	  cnt++;
	}
      }
      gsl_vector_free(temprates);

      // Set up minimizer
      /* Various possible algorithms,
	 NLOPT_GN_DIRECT
	 NLOPT_GN_DIRECT_L
	 NLOPT_GLOBAL_DIRECT_L_RAND
	 NLOPT_GLOBAL_DIRECT_NOSCAL, NLOPT_GLOBAL_DIRECT_L_NOSCAL
	 NLOPT_GN_ORIG_DIRECT, NLOPT_GN_ORIG_DIRECT_L
	 NLOPT_GN_CRS2_LM
	 NLOPT_G_MLSL_LDS, NLOPT_G_MLSL
	 NLOPT_GD_STOGO, or NLOPT_GD_STOGO_RAND
	 NLOPT_GN_ISRES
	 NLOPT_LN_COBYLA
	 NLOPT_LN_BOBYQA
	 NLOPT_LN_NEWUOA
	 NLOPT_LN_NEWUOA_BOUND
	 NLOPT_LN_PRAXIS
	 NLOPT_LN_NELDERMEAD
	 NLOPT_LN_SBPLX
	 NLOPT_LD_MMA
	 NLOPT_LD_SLSQP
	 NLOPT_LD_LBFGS
	 NLOPT_LD_TNEWTON_PRECOND_RESTART, NLOPT_LD_TNEWTON_PRECOND
	 NLOPT_LD_TNEWTON_RESTART, NLOPT_LD_TNEWTON
	 NLOPT_LD_VAR2, NLOPT_LD_VAR1
	 NLOPT_AUGLAG, NLOPT_AUGLAG_EQ
       */
      opt = nlopt_create(NLOPT_LN_SBPLX, nparams);
      nlopt_set_lower_bounds(opt, lb);
      nlopt_set_upper_bounds(opt, ub);
      nlopt_set_min_objective(opt, compute_dist_and_grad, d);
      //      nlopt_set_maxeval(opt, MAXEVAL);
      nlopt_set_ftol_abs(opt, 1e-15);

#ifdef LOCAL
      nlopt_opt opt_local;
      opt_local = nlopt_create(NLOPT_LN_NELDERMEAD, nparams);
      nlopt_set_lower_bounds(opt_local, lb);
      nlopt_set_upper_bounds(opt_local, ub);
      nlopt_set_min_objective(opt_local, compute_dist_and_grad, d);
      //      nlopt_set_maxeval(opt_local, MAXEVAL/10);
      nlopt_set_ftol_abs(opt_local, 1e-15);
#endif

      double * x = (double *) malloc(sizeof(double)*nparams);
      double bestfun = 1e200;
      double * bestxopt = (double *) malloc(sizeof(double)*nparams);

      for (uint rest=NRESTARTS; rest>0; rest--) {
        bool reset = false;
	d->count = 0;
	//random starting point
	if (bestfun < FTOL) break;
	for (uint kind=0; kind<numdemes; kind++) {
	  x[kind] = gsl_ran_flat(r, 1e-5, 1e-3);
	}
	for (uint kind=numdemes; kind<nparams; kind++) {
	  x[kind] = gsl_ran_flat(r, 1e-8, 1e-3);
	}
	int retcode = 0;
	try {
	  retcode = nlopt_optimize(opt, x, &minf);
	} catch (const char * detError) {
#ifdef DEBUG
	  cout << detError << endl;
#endif
	  rest++;
	  reset = true;
	}
	if (reset) continue;
	if (retcode < 0) {
	  printf("nlopt first step failed - %d!\n", retcode);
	}
	else {
#ifdef DEBUG
	  printf("Completed first step optimization in %d fevals.\n", d->count);
#endif
	  if (minf < bestfun) {
	    bestfun = minf;
	    memcpy(bestxopt, x, sizeof(double)*nparams);
	  }
	}
#ifdef LOCAL
	d->count = 0;
	try {
	  retcode = nlopt_optimize(opt_local, x, &minf);
	} catch (const char * detError) {
#ifdef DEBUG
	  cout << detError << endl;
#endif
	  rest++; 
	  reset = true;
	}
	if (reset) continue;
	if (retcode < 0) {
	  printf("nlopt second step failed - %d!\n", retcode);
	} else if (minf < bestfun) {
	  bestfun = minf;
	  memcpy(bestxopt, x, sizeof(double)*nparams);
	}
#ifdef DEBUG
	printf("Completed second step optimization with %d fevals.\n", d->count);
#endif
#endif
      }
      free(x);
      nlopt_destroy(opt);
#ifdef LOCAL
      nlopt_destroy(opt_local);
#endif
      //check for population mergers
      Ninv = gsl_vector_alloc(numdemes);
      memcpy(Ninv->data, bestxopt, sizeof(double)*numdemes);
      mtemp.clear();
      mtemp.insert(mtemp.end(), bestxopt+numdemes, bestxopt+nparams);
      vector<vector<int > > popdict = find_pop_merges(Ninv, mtemp, d->t, d->P0, \
						      merge_threshold, useMigration);
      if (popdict.size() < numdemes) {
	gsl_matrix * temp = converge_pops(popdict, d->P0);
	gsl_matrix_free(d->P0);
	d->P0 = gsl_matrix_alloc(temp->size1, temp->size2);
	gsl_matrix_memcpy(d->P0, temp);
	gsl_matrix_free(temp);
	reestimate = true;
	bestfun = 1e200;
	pdlist.push_back(popdict);
	numdemes = popdict.size();
#ifdef DEBUG
	cout << "\tre-estimating due to population merging.\nCurrent estimate: ";
        copy(bestxopt, bestxopt+nparams, ostream_iterator<double>(cout, " "));
	cout << endl;
#endif
      } else {
	vector<double> currxopt;
	for (uint kind=0; kind<numdemes; kind++) {
	  currxopt.push_back(1./bestxopt[kind]);
	}
	for (uint testind=numdemes; testind < nparams; testind++) {
	  if (*(bestxopt+testind) < 1e-10) {
	    currxopt.push_back(0.0);
	  } else {
	    currxopt.push_back(*(bestxopt+testind));
	  }
	}
	    
	reestimate=false;
	xopts.push_back(currxopt);
      }
#ifdef DEBUG
      cout << "Best function estimate " << bestfun << endl;
      if (bestfun < 1) {
	cout << "Nparams: " << nparams << "\t";
        copy(bestxopt, bestxopt+nparams, ostream_iterator<double>(cout, " "));
	cout << endl;       
      }
#endif
      free(bestxopt);
    } while(reestimate);
    gsl_matrix * m = gsl_matrix_calloc(numdemes, numdemes);
    uint cnt=0;
    for (uint row=0; row<numdemes; row++) {
      for (uint col=row+1; col<numdemes; col++) {
	m->data[row*numdemes+col] = m->data[col*numdemes+row] = mtemp[cnt];
	cnt++;
      }
    }
    gsl_matrix *Q = comp_pw_coal_cont(m, Ninv);
    gsl_matrix_scale(Q, d->t);
    gsl_matrix *P = expM(Q);
    gsl_matrix_free(Q);
    gsl_matrix * temp = conv_scrambling_matrix(P);
    gsl_matrix_free(P);
    gsl_matrix * temp2 = gsl_matrix_alloc(d->P0->size1, temp->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, d->P0, temp, 0.0, temp2);
    gsl_matrix_free(temp);
    gsl_matrix_free(d->P0);
    d->P0 = gsl_matrix_alloc(temp2->size1, temp2->size2);
    gsl_matrix_memcpy(d->P0, temp2);
    gsl_matrix_free(temp2);
    gsl_vector_free(Ninv);
    mtemp.clear();
  }
  return xopts;
}
