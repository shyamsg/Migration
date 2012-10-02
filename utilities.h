#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include "migrate.h"
#include <stdio.h>

/***************************************
This function takes in the npops, pops
of interest and sets the appropriate vectors 
in d
 ***************************************/
double * (vector<double> parms, cfmn_data *d, unsigned int npops, \
      vector<unsigned int> popsOfInt) 
{
  // This is the vector with the values of parms that belong to the 
  // popsOfInt.
  double * chosen;
  // otherParms is the vector with values of parms for the other pops
  d->otherParms.clear();
  // vector of index of other parms (ones not optimized)
  d->indexOthers.clear();
  // vector of index of optmized parms
  d->indexOpt.clear();
  // Compute the indices of the optimized pops
  int cnt = -1;
  chosen = (double *) malloc(sizeof(double)*(popsOfInt.size()*(popsOfInt.size()+1))/2);
  for (unsigned int t = 0; t < popsOfInt.size(); t++) {
    for (unsigned int u = t; u < popsOfInt.size(); u++) {
      d->indexOpt.insert((popsOfInt[t]*(2*npops-popsOfInt[t]-1))/2 + popsOfInt[u]);
      chosen[++cnt] = parms[(popsOfInt[t]*(2*npops-popsOfInt[t]-1))/2 + popsOfInt[u]];
    }
  }
  for (unsigned int t = 0; t < (npops*(npops+1))/2; t++) {
    if (d->indexOpt.find(t) == d->indexOpt.end()) {
      d->indexOthers.insert(t);
      d->otherParms.push_back(parms[t]);
    }
  }
  return chosen;
}

/*
    next_comb(vector<int> comb, int k, int n)
        Generates the next combination of n elements as k after comb

    comb => the previous combination ( use (0, 1, 2, ..., k) for first)
    k => the size of the subsets to generate
    n => the size of the original set

    Returns: 1 if a valid combination was found
        0, otherwise
*/
bool next_comb(vector<int> comb, int k, int n) {
    int i = k - 1;
    ++comb[i];
    while ((i >= 0) && (comb[i] >= n - k + 1 + i)) {
        --i;
        ++comb[i];
    }

    if (comb[0] > n - k) /* Combination (n-k, n-k+1, ..., n) reached */
        return false; /* No more combinations can be generated */

    /* comb now looks like (..., x, n, n, n, ..., n).
    Turn it into (..., x, x + 1, x + 2, ...) */
    for (i = i + 1; i < k; ++i)
        comb[i] = comb[i - 1] + 1;

    return true;
}

int main(int argc, char *argv[]) {
    int n = 10; /* The size of the set; for {1, 2, 3, 4} it's 4 */
    int k = 2; /* The size of the subsets; for {1, 2}, {1, 3}, ... it's 2 */
    vector<int> comb = vector<int>(k); /* comb[i] is the index of the i-th element in the
            combination */

    /* Setup comb for the initial combination */
    int i;
    for (i = 0; i < k; ++i)
        comb[i] = i;

    return 0;
}

#endif
