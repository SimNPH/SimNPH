#include <Rcpp.h>
using namespace Rcpp;

// all those functions assume that Tint is ordered and
// that Tint and lambda have compatible sizes
// this is checked in the R function that creates the calling R function

// [[Rcpp::export]]
NumericVector hazFunCpp (const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);

  int i, j;

  for(i=0; i < n; i++){
    for(j=0; j < m; j++){
      if(v[i] >= Tint[m-1-j]){
        result[i] = lambda[m-1-j];
        break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector cumhazFunCpp (const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);

  int i, j;

  for(i=0; i < n; i++){
    for(j=1; j < m; j++){
      result[i] += (lambda[j-1] * std::max(std::min(Tint[j], v[i]) - Tint[j-1], 0.) );
    }
    result[i] += std::max(v[i] - Tint[m-1], 0.) * lambda[m-1];
  }

  return result;
}

// [[Rcpp::export]]
NumericVector cdfFunCpp (const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);

  int i, j;

  for(i=0; i < n; i++){
    for(j=1; j < m; j++){
      result[i] += (lambda[j-1] * std::max(std::min(Tint[j], v[i]) - Tint[j-1], 0.) );
    }
    result[i] += std::max(v[i] - Tint[m-1], 0.) * lambda[m-1];
  }

  result=1-exp(-result);
  return result;
}

// [[Rcpp::export]]
NumericVector pdfFunCpp (const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);

  int i,j;

  double haz, ee;

  for(i=0; i < n; i++){
    ee = 0;
    haz = 0;

    for(j=0; j < m; j++){
      if(v[i] >= Tint[m-1-j]){
        haz = lambda[m-1-j];
        break;
      }
    }

    for(j=0; j < m; j++){
      ee += (lambda[j-1] * std::max(std::min(Tint[j], v[i]) - Tint[j-1], 0.) );
    }
    ee += std::max(v[i] - Tint[m-1], 0.) * lambda[m-1];

    result[i] = exp(-ee) * haz;
  }

  return result;
}


// [[Rcpp::export]]
NumericVector survFunCpp (const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);

  int i, j;

  for(i=0; i < n; i++){
    for(j=1; j < m; j++){
      result[i] += (lambda[j-1] * std::max(std::min(Tint[j], v[i]) - Tint[j-1], 0.) );
    }
    result[i] += std::max(v[i] - Tint[m-1], 0.) * lambda[m-1];
  }

  result=exp(-result);
  return result;
}

/*
// [[Rcpp::export]]
NumericVector quantFunCpp(const NumericVector& Tint, const NumericVector& lambda, const NumericVector& v){
  const int n = v.size();
  const int m = Tint.size();
  NumericVector result(n);
  NumericVector Qint(M);

  int i, j;

  // calculate at which quantiles the intervall limits lie
  for(j=0; j < m; j++){

  }

  // calculate the quantiles at the probabilities v
  for(i=0; i < n; i++){
    // find in which intervall the probability falls
    for(j=0; j < m; j++){
      if(Qint[j] < v[i]){
        break;
      }
    }
    result[i] = (Tint[j-1] - log(v[i] / Qint[j]) / lambda[j])
  }

  return result;
}
*/
