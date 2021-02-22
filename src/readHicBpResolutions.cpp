#include <fstream>
#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// returns whether or not this is valid HiC file
bool readMagicString(ifstream& fin) {
  string str;
  getline(fin, str, '\0' );
  return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

// [[Rcpp::export]]
NumericVector readHicBpResolutions(string hicFilename)
{
  ifstream fin(hicFilename, ios::in | ios::binary);

  if (!readMagicString(fin)) {
    cerr << "Hi-C magic string is missing, does not appear to be a hic file" << endl;
    fin.close();
    return NumericVector(0);
  }

  int version;
  fin.read((char*)&version, sizeof(int));
  if (version < 6) {
    cerr << "Version " << version << " no longer supported" << endl;
    fin.close();
    return NumericVector(0);
  }
  long master;
  fin.read((char*)&master, sizeof(long));
  string genome;
  getline(fin, genome, '\0' );
  int nattributes;
  fin.read((char*)&nattributes, sizeof(int));
  // reading and ignoring attribute-value dictionary
  for (int i=0; i<nattributes; i++) {
    string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  int nChrs;
  fin.read((char*)&nChrs, sizeof(int));
  // chromosome map for finding matrix
  for (int i=0; i<nChrs; i++) {
    string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*)&length, sizeof(int));
  }
  int nBpResolutions;
  fin.read((char*)&nBpResolutions, sizeof(int));
  NumericVector bpResolutions(nBpResolutions);
  for (int i=0; i<nBpResolutions; i++) {
    int resBP;
    fin.read((char*)&resBP, sizeof(int));
    bpResolutions[i] = resBP;
  }

  fin.close();

  return bpResolutions;
}
