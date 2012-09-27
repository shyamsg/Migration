#ifndef __UTILITIES_H__
#define __UTILITIES_H__

/***************************************
This function takes in the npops, pops
of interest and sets the appropriate vectors 
in d
 ***************************************/
vector<double> (vector<double> parms, cfmn_data *d, unsigned int npops, \
      vector<unsigned int> popsOfInt) 
{
  vector<double> chosen;
  d->otherParms.clear();
  d->indexOthers.clear();
  d->indexOpt.clear();

}


#endif
