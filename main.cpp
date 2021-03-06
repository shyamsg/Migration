#include <iostream>
#include "migrate.h"
#include <time.h>
#include <tclap/CmdLine.h>
#include <fstream>

#define MERGE_THRESHOLD 1e-2
using namespace std;

string ratefile;
string timefile;
string outfile = "testme";
int npop;
int skip;
bool logVal;
unsigned int poi;

/***********************************************************
This function reads the rates and times files and constructs
a matrix and vector out of these respectively. It also handles
the skip argument and omits the first 'skip' timeslices and
rates.
************************************************************/
gsl_matrix * readRatesAndTimes(vector<double> & times)
{
  gsl_matrix * rates = NULL;
  ifstream timestream(timefile.c_str());
  float timeGen;
  while(!timestream.eof()) {
    timeGen = -1;
    timestream >> timeGen;
    if (timeGen != -1)
      times.push_back(timeGen);
  }
  timestream.close();
  ifstream ratestream(ratefile.c_str());
  int comps = (npop*(npop+1))/2;
  double rate;
  rates = gsl_matrix_alloc(comps, times.size()-skip);
  int rowindex = -1;
  while(!ratestream.eof()) {
    rowindex += 1;
    for (int i = 0; i < (int)times.size(); i++) {
      rate = -1.0;
      ratestream >> rate;
      if (rate == -1) break;
      if ((i - skip) > -1) {
	gsl_matrix_set(rates, rowindex, i-skip, rate);
      }
    }
  }
  times.erase(times.begin(), times.begin()+skip);
  cout << times.size() << " " << skip << endl;
  return rates;
}

/***********************************************************
This function is used to write the output from the migration
estimation algorithm. It writes out 4 things, the slice number,
the time in gens, the current version of the population merger
state and the estimates for the remaining populations.

************************************************************/
void writeOutput(vector<double> & times, vector<vector<vector<int> > > & pdout, \
		 vector<vector<double> > & estParms)
{
  string outEsts = outfile + ".ests";
  string outPd = outfile +".pd";
  ofstream out(outEsts.c_str());
  ofstream pdo(outPd.c_str());
#ifdef DEBUG
  cout << "Writing output " << times.size() << "\t" << pdout.size() << endl;
#endif
  for (uint i=0; i < times.size(); i++) {
    out << i << "\t" << times[i] << "\t";
    copy(estParms[i].begin(), estParms[i].end(), ostream_iterator<double>(out, ","));
    out << "\n";
  }
  for (uint i=0; i < pdout.size(); i++) {
    for (uint j=0; j < pdout[i].size(); j++) {
      copy(pdout[i][j].begin(), pdout[i][j].end(), ostream_iterator<int>(pdo, ","));
      pdo << "|";
    }
    pdo << "\n";
  }
  pdo.close();
  out.close();
}

/***********************************************************
This is a command line parser function. It parses cla using
the TCLAP classes. 
***********************************************************/
void parseCmdLine(int argc, char **argv)
{
  try {  
    TCLAP::CmdLine cmd("Migrate rate estimator", ' ', "0.1");
    // Various arguments that this module takes are
    // rate file.
    // time file
    // output file.
    // skip number of initial slices
    // pattern -- add later

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<int> skipArg("s", "skip", "Initial time slices to skip", false, 0, " integer > 0 ");
    cmd.add(skipArg);
    TCLAP::ValueArg<unsigned int> npArg("n", "npoi", "Number of pops to solve at once", false, 0, " integer >= 2 ");
    cmd.add(npArg);
    TCLAP::SwitchArg logArg("l", "logCost", "Log the cost function", false);
    cmd.add(logArg);
    TCLAP::ValueArg<int> popArg("p", "pop", "Number of populations", true, 2, "integer > 0 ");
    cmd.add(popArg);
    TCLAP::ValueArg<std::string> rateArg("r","rate","Coalescent rates file", true, "rates.txt", "filename");
    cmd.add(rateArg);
    TCLAP::ValueArg<std::string> timeArg("t","time","Time slice file", true, "times.txt", "filename");
    cmd.add(timeArg);
    TCLAP::ValueArg<std::string> outArg("o","out","Output file prefix", false, "test", "filename");
    cmd.add(outArg);

    // Parse the argv array.
    cmd.parse( argc, argv );
    // Get the value parsed by each arg. 
    ratefile = rateArg.getValue();
    timefile = timeArg.getValue();
    outfile = outArg.getValue();
    skip = skipArg.getValue();
    npop = popArg.getValue();
    logVal = logArg.getValue();
    poi = npArg.getValue();
    if (poi == 0) poi = npop;
  } catch (TCLAP::ArgException &e)  { // catch any exceptions 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
  }
}

int main(int argc, char **argv)
{
  vector<double> times;
  gsl_matrix * conv;
#ifndef TEST
  parseCmdLine(argc, argv);
  conv = readRatesAndTimes(times);
  gsl_matrix_print(conv);
#else
  /*************DEBUG_CODE***************/
  uint id = (unsigned int)atoi(argv[1]);
  logVal = (argc > 2) ? true: false;
  poi= 2;
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
  conv = compute_pw_coal_rates(Ns, ms, ts, popmap);
  gsl_matrix_print(conv);
  times = ts;
  cout << "The computed rates have dimension " << conv->size1 << "x" << conv->size2 << endl;
  /*************DEBUG_CODE***************/
#endif
  time_t start, end;
  time(&start);
  vector<vector<vector<int> > > pdout;
  vector<vector<double> > estParms = comp_params(conv, times, pdout, logVal, poi, MERGE_THRESHOLD, true);
  //  cout << estParms.size() << endl;
  for (uint i=0; i<estParms.size(); i++) {
    cout << "Estimates for slice number " << i << endl;
    copy(estParms[i].begin(), estParms[i].end(), ostream_iterator<double>(cout, " "));
    cout << endl;
  }
  gsl_matrix_free(conv);
  writeOutput(times, pdout, estParms);
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
