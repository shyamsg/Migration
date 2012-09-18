#include <iostream>
#include "migrate.h"
#include <time.h>

void parseCmdLine(){};

int main(int argc, char **argv)
{  
  gsl_matrix * conv;
  uint id = (unsigned int)atoi(argv[1]);
  vector<vector<double> > Ns = vector<vector<double> >(id);
  vector<vector<double> > ms = vector<vector<double> >(id);
  vector<double> ts = vector<double>();
  vector<vector<vector<int> > > popmap = vector<vector<vector<int> > >(id);
  if (id > 9) {
    Ns[id-10].push_back(10000.0), Ns[id-10].push_back(15000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 10000.0), Ns[id-10].push_back( 20000.0), Ns[id-10].push_back( 20000.0);
    ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back( 1.0e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back(0.0e+00), ms[id-10].push_back( 0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(1.5e-05), ms[id-10].push_back( 1.0e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back(1.0e-05), ms[id-10].push_back( 0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back( 1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back( 0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back( 1.5e-05), ms[id-10].push_back(1.5e-05), ms[id-10].push_back(0.0e+00), ms[id-10].push_back( 0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(1.5e-05), ms[id-10].push_back( 0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back(0.0e+00), ms[id-10].push_back( 1.2e-05), ms[id-10].push_back(1.2e-05), ms[id-10].push_back(1.2e-05), ms[id-10].push_back( 1.2e-05), ms[id-10].push_back(1.2e-05), ms[id-10].push_back(1.2e-05);    
    ts.push_back(1500);
  }
  if (id > 8) {
    Ns[id-9].push_back(10000.0), Ns[id-9].push_back( 15000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 10000.0), Ns[id-9].push_back( 20000.0); 
    ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(1e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.5e-05), ms[id-9].push_back(0.0), ms[id-9].push_back(0.0), ms[id-9].push_back(1.2e-05), ms[id-9].push_back(1.2e-05), ms[id-9].push_back(1.2e-05);
    ts.push_back(1800);
    if (id-9 > 0) {
      popmap[id-9] = vector<vector<int> >(9); popmap[id-9][0].push_back(0), popmap[id-9][1].push_back(1), popmap[id-9][2].push_back(2), popmap[id-9][3].push_back(3), popmap[id-9][4].push_back(4), popmap[id-9][5].push_back(5), popmap[id-9][6].push_back(6), popmap[id-9][7].push_back(7), popmap[id-9][8].push_back(8), popmap[id-9][8].push_back(9);
    }
  }
  if (id > 7) {
    Ns[id-8].push_back(10000.0), Ns[id-8].push_back( 15000.0), Ns[id-8].push_back( 10000.0), Ns[id-8].push_back( 10000.0), Ns[id-8].push_back( 10000.0), Ns[id-8].push_back( 10000.0), Ns[id-8].push_back( 10000.0), Ns[id-8].push_back( 10000.0); 
    ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(0.0), ms[id-8].push_back(0.0), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(1e-05), ms[id-8].push_back(0.0), ms[id-8].push_back(0.0), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(0.0), ms[id-8].push_back(0.0), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(0.0), ms[id-8].push_back(0.0), ms[id-8].push_back(1.5e-05), ms[id-8].push_back(0.0), ms[id-8].push_back(0.0), ms[id-8].push_back(1.2e-05), ms[id-8].push_back(1.2e-05), ms[id-8].push_back(1.2e-05);
    ts.push_back(1800);
    if (id-8 > 0) {
      popmap[id-8] = vector<vector<int> >(8); popmap[id-8][0].push_back(0), popmap[id-8][1].push_back(1), popmap[id-8][2].push_back(2), popmap[id-8][2].push_back(3), popmap[id-8][3].push_back(4), popmap[id-8][4].push_back(5), popmap[id-8][5].push_back(6), popmap[id-8][6].push_back(7), popmap[id-8][7].push_back(8);
    }
  }
  if (id > 6) {
    Ns[id-7].push_back(10000.0), Ns[id-7].push_back( 10000.0), Ns[id-7].push_back( 10000.0), Ns[id-7].push_back( 10000.0), Ns[id-7].push_back( 10000.0), Ns[id-7].push_back( 12000.0), Ns[id-7].push_back( 12000.0); 
    ms[id-7].push_back(1.5e-05), ms[id-7].push_back(1.5e-05), ms[id-7].push_back(1e-05), ms[id-7].push_back(1e-05), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(1.5e-05), ms[id-7].push_back(1e-05), ms[id-7].push_back(1e-05), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(1.5e-05), ms[id-7].push_back(1.5e-05), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(1.5e-05), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(0.0), ms[id-7].push_back(1.2e-05);
    ts.push_back(2000);
    if (id-7 > 0) {
      popmap[id-7] = vector<vector<int> >(7); popmap[id-7][0].push_back(0), popmap[id-7][1].push_back(1), popmap[id-7][2].push_back(2), popmap[id-7][3].push_back(3), popmap[id-7][4].push_back(4), popmap[id-7][4].push_back(5), popmap[id-7][5].push_back(6), popmap[id-7][6].push_back(7);
    }
  }
  if (id > 5) {
    Ns[id-6].push_back(10000.0), Ns[id-6].push_back( 10000.0), Ns[id-6].push_back( 10000.0), Ns[id-6].push_back( 10000.0), Ns[id-6].push_back( 12000.0), Ns[id-6].push_back( 12000.0); 
    ms[id-6].push_back(3e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(0), ms[id-6].push_back(1e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(0), ms[id-6].push_back(1e-05), ms[id-6].push_back(1e-05), ms[id-6].push_back(0), ms[id-6].push_back(1e-05), ms[id-6].push_back(0), ms[id-6].push_back(0);
    ts.push_back(2500);
    if (id-6 > 0) {
      popmap[id-6] = vector<vector<int> >(6); popmap[id-6][0].push_back(0), popmap[id-6][1].push_back(1), popmap[id-6][1].push_back(2), popmap[id-6][2].push_back(3), popmap[id-6][3].push_back(4), popmap[id-6][4].push_back(5), popmap[id-6][5].push_back(6);
    }
  }
  if (id > 4) {
    Ns[id-5].push_back(10000.0), Ns[id-5].push_back( 10000.0), Ns[id-5].push_back( 10000.0), Ns[id-5].push_back( 12000.0), Ns[id-5].push_back( 12000.0);
    ms[id-5].push_back(2e-05), ms[id-5].push_back(1e-05), ms[id-5].push_back(1e-05), ms[id-5].push_back(0), ms[id-5].push_back(1e-05), ms[id-5].push_back(1e-05), ms[id-5].push_back(0), ms[id-5].push_back(1e-05), ms[id-5].push_back(0), ms[id-5].push_back(2e-05);
    ts.push_back(2500);
    if (id-5 > 0) {
      popmap[id-5] = vector<vector<int> >(5); popmap[id-5][0].push_back(0), popmap[id-5][0].push_back(1), popmap[id-5][1].push_back(2), popmap[id-5][2].push_back(3), popmap[id-5][3].push_back(4), popmap[id-5][4].push_back(5);
    }
  }
  if (id > 3) {
    Ns[id-4].push_back(10000.0), Ns[id-4].push_back( 10000.0), Ns[id-4].push_back( 12000.0), Ns[id-4].push_back( 12000.0);
    ms[id-4].push_back(1e-05), ms[id-4].push_back(1e-05), ms[id-4].push_back(1e-05), ms[id-4].push_back(1e-05), ms[id-4].push_back(1e-06), ms[id-4].push_back(4e-05);
    ts.push_back(2500);  
    if (id-4 > 0) {
      popmap[id-4] = vector<vector<int> >(4); popmap[id-4][0].push_back(0), popmap[id-4][0].push_back(1), popmap[id-4][1].push_back(2), popmap[id-4][2].push_back(3), popmap[id-4][3].push_back(4);
    }
  }
  if (id > 2) {
    Ns[id-3].push_back(1e4);Ns[id-3].push_back(1e4);Ns[id-3].push_back(1.2e4);
    ms[id-3].push_back(2e-5);ms[id-3].push_back(1e-5);ms[id-3].push_back(1e-5);
    ts.push_back(2500);
    if (id-3 > 0) {
      popmap[id-3] = vector<vector<int> >(3); popmap[id-3][0].push_back(0), popmap[id-3][1].push_back(1), popmap[id-3][2].push_back(2), popmap[id-3][2].push_back(3);
    }
  }
  if (id > 1) {
    Ns[id-2].push_back(1e4);Ns[id-2].push_back(1.2e4);
    ms[id-2].push_back(2e-5);
    ts.push_back(2500);
    if (id-2 > 0) {
      popmap[id-2] = vector<vector<int> >(2); popmap[id-2][0].push_back(0), popmap[id-2][0].push_back(1), popmap[id-2][1].push_back(2);
    }
  }
  if (id > 0) {
    Ns[id-1].push_back(1e4);
    ts.push_back(5000);
    if (id-1 > 0) {
      popmap[id-1] = vector<vector<int> >(1); popmap[id-1][0].push_back(0), popmap[id-1][0].push_back(1);
    }
  }
  time_t start, end;
  time(&start);
  conv = compute_pw_coal_rates(Ns, ms, ts, popmap);
  gsl_matrix_print(conv);
  //  cout << "The computed rates have dimension " << conv->size1 << "x" << conv->size2 << endl;
  vector<vector<vector<int> > > testlist;
  vector<vector<double> > estParms = comp_params(conv, ts, testlist, 1e-3, true);
  //  cout << estParms.size() << endl;
  for (uint i=0; i<estParms.size(); i++) {
    cout << "Estimates for slice number " << i << endl;
    copy(estParms[i].begin(), estParms[i].end(), ostream_iterator<double>(cout, " "));
    cout << endl;
  }
  gsl_matrix_free(conv);
  /*
  conv = gsl_matrix_alloc(2,2);
  conv->data[0] = 1.3;
  conv->data[1] = 1.7;
  conv->data[2] = 2.1;
  conv->data[3] = 1.4;

  gsl_matrix_print(conv);
  gsl_matrix * test = gsl_matrix_alloc(2,2);
  invert_matrix(conv, test);
  gsl_matrix_print(test);
  */
  time(&end);
  double diff = difftime(end, start);
  printf("Time Elapsed: %.2lf\n", diff);
  return 0;
}
