#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rebinBigwig(DataFrame bigwig, DataFrame bins)
{
  NumericVector bw_start = bigwig["start"];
  NumericVector bw_end = bigwig["end"];
  NumericVector bw_score = bigwig["score"];
  NumericVector bins_start = bins["start"];
  NumericVector bins_end = bins["end"];
  NumericVector max_scores(bins.nrows());
  int bw_i = 0, b_i = 0;
  while(bw_i < bigwig.nrows() && b_i < bins.nrows())
  {
    if(bw_end[bw_i]<bins_start[b_i])
      bw_i++;
    else if(bins_end[b_i]<bw_start[bw_i])
      b_i++;
    else
    {
      max_scores[b_i] = fmax(max_scores[b_i], bw_score[bw_i]);
      if(bw_end[bw_i]<bins_end[b_i])
        bw_i++;
      else if(bins_end[b_i]<bw_end[bw_i])
        b_i++;
      else if(bw_end[bw_i]==bins_end[b_i])
      {
        bw_i++;
        b_i++;
      }
    }
  }
  return max_scores;
}
