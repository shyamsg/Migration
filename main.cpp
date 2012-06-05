#include <iostream>
#include "migrate.h"

void parseCmdLine(){};

int main()
{
  vector<vector<int> > pd= vector<vector<int> >(1);
  pd[0] = vector<int>(); pd[0].push_back(0); pd[0].push_back(1);
  //  pd[1] = vector<int>(); pd[1].push_back(1); pd[1].push_back(2);
  
  vector<vector<int> > ptc = pop_to_col(pd, 2);
  for(uint i=0; i<ptc.size(); i++) {
    printf("%d:\t", i);
    copy(ptc[i].begin(), ptc[i].end(), ostream_iterator<int>(cout, " "));
    cout << endl;
  }
  
  gsl_matrix * testSc = gsl_matrix_calloc(4,4);
  testSc->data[0] = 0.31;
  testSc->data[1] = 0.3;
  testSc->data[2] = 0.29;
  testSc->data[3] = 0.1;
  testSc->data[4] = 0.15;
  testSc->data[5] = 0.24;
  testSc->data[6] = 0.28;
  testSc->data[7] = 0.33;
  testSc->data[8] = 0.42;
  testSc->data[9] = 0.47;
  testSc->data[10] = 0.1;
  testSc->data[11] = 0.01;
  testSc->data[15] = 1.0;
  
  gsl_matrix * conv = converge_pops(pd, testSc);
  for (uint i=0; i<conv->size1; i++) {
    cout << "[ ";
    for (uint j=0; j<conv->size2; j++) {
      cout << conv->data[i*conv->size2+j] << " ";
    }
    cout << "]" << endl;
  }
  gsl_matrix_free(conv);
  conv = expM(testSc);
  gsl_matrix_print(conv);
  vector<vector<double> > Ns = vector<vector<double> >(10);
  Ns[0].push_back(10000.0), Ns[0].push_back(15000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 10000.0), Ns[0].push_back( 20000.0), Ns[0].push_back( 20000.0);
  Ns[1].push_back(10000.0), Ns[1].push_back( 15000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 10000.0), Ns[1].push_back( 20000.0); 
  Ns[2].push_back(10000.0), Ns[2].push_back( 15000.0), Ns[2].push_back( 10000.0), Ns[2].push_back( 10000.0), Ns[2].push_back( 10000.0), Ns[2].push_back( 10000.0), Ns[2].push_back( 10000.0), Ns[2].push_back( 10000.0); 
  Ns[3].push_back(10000.0), Ns[3].push_back( 10000.0), Ns[3].push_back( 10000.0), Ns[3].push_back( 10000.0), Ns[3].push_back( 10000.0), Ns[3].push_back( 12000.0), Ns[3].push_back( 12000.0); 
  Ns[4].push_back(10000.0), Ns[4].push_back( 10000.0), Ns[4].push_back( 10000.0), Ns[4].push_back( 10000.0), Ns[4].push_back( 12000.0), Ns[4].push_back( 12000.0); 
  Ns[5].push_back(10000.0), Ns[5].push_back( 10000.0), Ns[5].push_back( 10000.0), Ns[5].push_back( 12000.0), Ns[5].push_back( 12000.0);
  Ns[6].push_back(10000.0), Ns[6].push_back( 10000.0), Ns[6].push_back( 12000.0), Ns[6].push_back( 12000.0);
  Ns[7].push_back(1e4);Ns[7].push_back(1e4);Ns[7].push_back(1.2e4);
  Ns[8].push_back(1e4);Ns[8].push_back(1.2e4);
  Ns[9].push_back(1e4);
  vector<vector<double> > ms = vector<vector<double> >(10);
  ms[0].push_back(1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back(1.0e-05), ms[0].push_back(1.0e-05), ms[0].push_back(1.0e-05), ms[0].push_back(0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back(1.0e-05), ms[0].push_back( 1.0e-05), ms[0].push_back(1.0e-05), ms[0].push_back(0.0e+00), ms[0].push_back( 0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(1.5e-05), ms[0].push_back( 1.0e-05), ms[0].push_back(1.0e-05), ms[0].push_back(1.0e-05), ms[0].push_back( 0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back( 1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back( 0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back( 1.5e-05), ms[0].push_back(1.5e-05), ms[0].push_back(0.0e+00), ms[0].push_back( 0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(1.5e-05), ms[0].push_back( 0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back(0.0e+00), ms[0].push_back( 1.2e-05), ms[0].push_back(1.2e-05), ms[0].push_back(1.2e-05), ms[0].push_back( 1.2e-05), ms[0].push_back(1.2e-05), ms[0].push_back(1.2e-05);
  ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.5e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(1e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.5e-05), ms[1].push_back(1.5e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.5e-05), ms[1].push_back(0.0), ms[1].push_back(0.0), ms[1].push_back(1.2e-05), ms[1].push_back(1.2e-05), ms[1].push_back(1.2e-05);
  ms[2].push_back(1.5e-05), ms[2].push_back(1.5e-05), ms[2].push_back(1e-05), ms[2].push_back(1e-05), ms[2].push_back(1e-05), ms[2].push_back(0.0), ms[2].push_back(0.0), ms[2].push_back(1.5e-05), ms[2].push_back(1e-05), ms[2].push_back(1e-05), ms[2].push_back(1e-05), ms[2].push_back(0.0), ms[2].push_back(0.0), ms[2].push_back(1.5e-05), ms[2].push_back(1.5e-05), ms[2].push_back(1.5e-05), ms[2].push_back(0.0), ms[2].push_back(0.0), ms[2].push_back(1.5e-05), ms[2].push_back(1.5e-05), ms[2].push_back(0.0), ms[2].push_back(0.0), ms[2].push_back(1.5e-05), ms[2].push_back(0.0), ms[2].push_back(0.0), ms[2].push_back(1.2e-05), ms[2].push_back(1.2e-05), ms[2].push_back(1.2e-05);
  ms[3].push_back(1.5e-05), ms[3].push_back(1.5e-05), ms[3].push_back(1e-05), ms[3].push_back(1e-05), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(1.5e-05), ms[3].push_back(1e-05), ms[3].push_back(1e-05), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(1.5e-05), ms[3].push_back(1.5e-05), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(1.5e-05), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(0.0), ms[3].push_back(1.2e-05);
  ms[4].push_back(3e-05), ms[4].push_back(1e-05), ms[4].push_back(1e-05), ms[4].push_back(1e-05), ms[4].push_back(0), ms[4].push_back(1e-05), ms[4].push_back(1e-05), ms[4].push_back(1e-05), ms[4].push_back(0), ms[4].push_back(1e-05), ms[4].push_back(1e-05), ms[4].push_back(0), ms[4].push_back(1e-05), ms[4].push_back(0), ms[4].push_back(0);
  ms[5].push_back(2e-05), ms[5].push_back(1e-05), ms[5].push_back(1e-05), ms[5].push_back(0), ms[5].push_back(1e-05), ms[5].push_back(1e-05), ms[5].push_back(0), ms[5].push_back(1e-05), ms[5].push_back(0), ms[5].push_back(2e-05);
  ms[6].push_back(1e-05), ms[6].push_back(1e-05), ms[6].push_back(1e-05), ms[6].push_back(1e-05), ms[6].push_back(1e-06), ms[6].push_back(4e-05);
  ms[7].push_back(2e-5);ms[7].push_back(1e-5);ms[7].push_back(1e-5);
  ms[8].push_back(2e-5);
  vector<double> ts = vector<double>();
  ts.push_back(1500), ts.push_back(1800), ts.push_back(1800), ts.push_back(2000), ts.push_back(2500), ts.push_back(2500), ts.push_back(2500), ts.push_back(2500), ts.push_back(2500), ts.push_back(5000);
  vector<vector<vector<int> > > popmap = vector<vector<vector<int> > >(10);
  //popmap[0] = vector<vector<int> >(); popmap[0][0] = vector<int>(); //empty vector for 1st popmap entry
  popmap[1] = vector<vector<int> >(9); popmap[1][0].push_back(0), popmap[1][1].push_back(1), popmap[1][2].push_back(2), popmap[1][3].push_back(3), popmap[1][4].push_back(4), popmap[1][5].push_back(5), popmap[1][6].push_back(6), popmap[1][7].push_back(7), popmap[1][8].push_back(8), popmap[1][8].push_back(9);
  popmap[2] = vector<vector<int> >(8); popmap[2][0].push_back(0), popmap[2][1].push_back(1), popmap[2][2].push_back(2), popmap[2][2].push_back(3), popmap[2][3].push_back(4), popmap[2][4].push_back(5), popmap[2][5].push_back(6), popmap[2][6].push_back(7), popmap[2][7].push_back(8);
  popmap[3] = vector<vector<int> >(7); popmap[3][0].push_back(0), popmap[3][1].push_back(1), popmap[3][2].push_back(2), popmap[3][3].push_back(3), popmap[3][4].push_back(4), popmap[3][4].push_back(5), popmap[3][5].push_back(6), popmap[3][6].push_back(7);
  popmap[4] = vector<vector<int> >(6); popmap[4][0].push_back(0), popmap[4][1].push_back(1), popmap[4][1].push_back(2), popmap[4][2].push_back(3), popmap[4][3].push_back(4), popmap[4][4].push_back(5), popmap[4][5].push_back(6);
  popmap[5] = vector<vector<int> >(5); popmap[5][0].push_back(0), popmap[5][0].push_back(1), popmap[5][1].push_back(2), popmap[5][2].push_back(3), popmap[5][3].push_back(4), popmap[5][4].push_back(5);
  popmap[6] = vector<vector<int> >(4); popmap[6][0].push_back(0), popmap[6][0].push_back(1), popmap[6][1].push_back(2), popmap[6][2].push_back(3), popmap[6][3].push_back(4);
  popmap[7] = vector<vector<int> >(3); popmap[7][0].push_back(0), popmap[7][1].push_back(1), popmap[7][2].push_back(2), popmap[7][2].push_back(3);
  popmap[8] = vector<vector<int> >(2); popmap[8][0].push_back(0), popmap[8][0].push_back(1), popmap[8][1].push_back(2);
  popmap[9] = vector<vector<int> >(1); popmap[9][0].push_back(0), popmap[9][0].push_back(1);
  
  gsl_matrix_free(conv);
  conv = compute_pw_coal_rates(Ns, ms, ts, popmap);
  gsl_matrix_print(conv);
  cout << "The computed rates have dimension " << conv->size1 << "x" << conv->size2 << endl;
  vector<vector<vector<int> > > testlist;
  vector<vector<double> > estParms = comp_params(conv, ts, testlist, 1e-2, true);
  cout << estParms.size() << endl;
  for (uint i=0; i<estParms.size(); i++) {
    cout << "Estimates for slice number " << i << endl;
    for (uint j=0; j<estParms[i].size(); j++) {
      cout << estParms[i][j] << " ";
    }
    cout << endl;
  }
  
  return 0;
}
