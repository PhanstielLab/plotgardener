#include <Rcpp.h>
using namespace Rcpp;

int rowNum(NumericVector v, NumericMatrix x, int maxRows, int rowCol, double wig) {

  // genomic coord start and stop, incorporating wiggle (in bp)
  int start = std::min(v[0], v[1]) - wig;
  int stop = std::max(v[0], v[1]) + wig;

  // check every row value
  for (int j = 1; j <= maxRows; j++){
    // initialize a row counter
    int rowCounter = 0;
    // go through dataframe to check for approriate row values
    int dfRows = x.nrow();
    for (int row = 0; row < dfRows; row++){
      NumericVector v = x(row,_);

      if (((v[rowCol] == j) &
          (((start >= v[0]) & (start <= v[1])) |
          ((stop >= v[0]) & (stop <= v[1])) |
          ((start <= v[0]) & (stop >= v[1]))))){

        rowCounter++;
      }
    }

    if (rowCounter == 0){

      return j;

    }

  }

  return 0;

}


// [[Rcpp::export]]
NumericMatrix checkRow(NumericMatrix x, int maxRows, int rowCol, double wig){

  // go through each row of the dataframe
  int dfRows = x.nrow();
  for (int i = 0; i < dfRows; ++i){
    // get the row
    NumericVector row = x(i,_);
    // determine its rownumber from rowNum function
    int rowNumber = rowNum(row, x, maxRows, rowCol, wig);
    // update whole dataframe with that value's rownumber
    x(i,rowCol) = rowNumber;
  }

  return x;

}
