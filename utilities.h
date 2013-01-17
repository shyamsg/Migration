#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include "migrate.h"
#include <stdio.h>

/***************************************
This function takes in the npops, pops
of interest and sets the appropriate vectors 
in d
 ***************************************/
void genCombo(double * parms, cfnm_data * d, unsigned int npops, int * popsOfInt, double *chosen, unsigned int numPopsInt);
/*
    next_comb(vector<int> comb, int k, int n)
        Generates the next combination of n elements as k after comb

    comb => the previous combination ( use (0, 1, 2, ..., k) for first)
    k => the size of the subsets to generate
    n => the size of the original set

    Returns: 1 if a valid combination was found
        0, otherwise
*/
bool next_comb(int *comb, int k, int n);
/*
int main(int argc, char *argv[]) {
    int n = 10; / * The size of the set; for {1, 2, 3, 4} it's 4 * /
    int k = 2; / * The size of the subsets; for {1, 2}, {1, 3}, ... it's 2 * /
    vector<int> comb = vector<int>(k); / * comb[i] is the index of the i-th element in the
            combination * /

    / * Setup comb for the initial combination * /
    int i;
    for (i = 0; i < k; ++i)
        comb[i] = i;

    return 0;
}
*/
#endif
