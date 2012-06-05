#include <iostream>
#include "migrate.h"

void parseCmdLine()

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
  vector<vector<double> > Ns = vector<vector<double> >(3);
  Ns[0].push_back(1e4);Ns[0].push_back(1e4);Ns[0].push_back(1.2e4);
  Ns[1].push_back(1e4);Ns[1].push_back(1.2e4);
  Ns[2].push_back(1e4);
  vector<vector<double> > ms = vector<vector<double> >(3);
  ms[0].push_back(2e-5);ms[0].push_back(1e-5);ms[0].push_back(1e-5);
  ms[1].push_back(2e-5);
  vector<double> ts = vector<double>();
  ts.push_back(2.5e3); ts.push_back(2.5e3); ts.push_back(5e3);
  vector<vector<vector<int> > > popmap = vector<vector<vector<int> > >(3);
  //popmap[0] = vector<vector<int> >(); popmap[0][0] = vector<int>(); //empty vector for 1st popmap entry
  popmap[1] = vector<vector<int> >(2); popmap[1][0].push_back(0); popmap[1][0].push_back(1); popmap[1][1].push_back(2);
  popmap[2] = vector<vector<int> >(1); popmap[2][0].push_back(0); popmap[2][0].push_back(1);
  
  gsl_matrix_free(conv);
  conv = compute_pw_coal_rates(Ns, ms, ts, popmap);
  gsl_matrix_print(conv);
  vector<vector<vector<int> > > yahoo = vector<vector<vector<int> > >(2);
  yahoo[0] = vector<vector<int> >(2); yahoo[0][0].push_back(0); yahoo[0][0].push_back(1); yahoo[0][1].push_back(2);
  yahoo[1] = vector<vector<int> >(1); yahoo[1][0].push_back(0); yahoo[1][0].push_back(1);	
  vector<vector<int> > test = make_merged_pd(yahoo);
  for (uint i=0; i<test.size(); i++) {
    cout << "[ ";
    for (uint j=0; j<test[i].size(); j++) {
      cout << test[i][j] << " ";
    }
    cout << "]" << endl;
  }
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
