#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double C_gini(NumericVector x){
  std::map<double, double> counts;
  R_xlen_t n = x.size();
  // NumericVector::iterator i;
  for (NumericVector::iterator i = x.begin(); i != x.end(); ++i) {
    counts[ *i ]++;
  }

  double sqsum = 0.0;
  for (std::map<double, double>::iterator it = counts.begin(); it != counts.end(); ++it)  {
    sqsum += 1.0 * (it->second) * (it->second);
  }
  return (1.0 - sqsum/(1.0*n*n));
}


// [[Rcpp::export]]
double C_ginicorr(NumericVector x, double k){
  double eps = 1e-5;
  if(std::abs(k - 1.0) < eps) return 1.0;

  return C_gini(x) / (1.0-1.0/k);
}

double C_sqrtginicorr(NumericVector x, double k){
  double eps = 1e-5;
  if(std::abs(k - 1.0) < eps) return 1.0;

  return (1 - sqrt(1 - C_gini(x))) / (1-1.0/k);
}

// [[Rcpp::export]]
double C_cfunc(NumericVector x, NumericVector y, double c, double k, bool sqrtflag){
  double varx, vary;
  if(sqrtflag){
    varx = C_sqrtginicorr(x, k);
    vary = C_sqrtginicorr(y, k);
  } else {
  varx = C_ginicorr(x, k);
  vary = C_ginicorr(y, k);
  }

  return(2*sqrt(varx * vary) + c)/(varx + vary + c);
}

// [[Rcpp::export]]
double C_meansfunc(NumericVector x, NumericVector y, double c){
  //R_xlen_t n = x.size();
  if (x.size() != y.size()) Rcpp::stop("X and Y must have the same length.");
  std::map<double, double> countsx;
  std::map<double, double> countsy;
  NumericVector::iterator x_i, y_i;
  for (x_i = x.begin(),  y_i = y.begin();
       x_i != x.end() && y_i != y.end(); ++x_i, ++y_i) {
    countsx[ *x_i ]++;
    countsy[ *y_i ]++;
  }

  double sqsum = 0.0;
  for (std::map<double,double>::iterator it = countsx.begin(); it != countsx.end(); ++it)  {
    sqsum += 1.0 * (it->second) * (it->second);
  }
  for (std::map<double,double>::iterator it = countsy.begin(); it != countsy.end(); ++it)  {
    sqsum += 1.0 * (it->second) * (it->second);
  }

  double xysum = 0.0;
  std::map<double,double>::iterator il = countsx.begin();
  std::map<double,double>::iterator ir = countsy.begin();
  while (il != countsx.end() && ir != countsy.end())
  {
    if (il->first < ir->first)
      ++il;
    else if (ir->first < il->first)
      ++ir;
    else
    {
      xysum += (il->second) * (ir->second);
      ++il;
      ++ir;
    }
  }


  return (2*(xysum) + c)/(sqsum + c);
}

// [[Rcpp::export]]
double C_Cohen(NumericVector x, NumericVector y){
  R_xlen_t n = x.size();
  if (x.size() != y.size()) Rcpp::stop("X and Y must have the same length.");
  NumericMatrix xy(n, 2);
  xy.column(0) = x;
  xy.column(1) = y;
  std::map<double, double> countsx;
  std::map<double, double> countsy;
  std::map<double, double> countsxy;

  countsx.clear();
  countsy.clear();

  NumericVector::iterator x_i, y_i;
  for (x_i = x.begin(),  y_i = y.begin();
       x_i != x.end() && y_i != y.end(); ++x_i, ++y_i) {
    countsx[ *x_i ]++;
    countsy[ *y_i ]++;
    if ( *x_i == *y_i ){
      countsxy[ *x_i ]++;
    }

  }

  double xxyysum = 0.0;
  std::map<double,double>::iterator il = countsx.begin();
  std::map<double,double>::iterator ir = countsy.begin();
  while (il != countsx.end() && ir != countsy.end())
  {
    if (il->first < ir->first)
      ++il;
    else if (ir->first < il->first)
      ++ir;
    else
    {
      xxyysum += 1.00 * (il->second) * (ir->second);
      ++il;
      ++ir;
    }
  }
  double xysum = 0.0;

  for (std::map<double, double>::iterator it = countsxy.begin(); it != countsxy.end(); ++it)  {
    xysum  += 1.0 * (it->second);
    }



  double pe = xxyysum / (1.0*n*n);
  double po = xysum / (1.0*n);

  if ((1.0-pe) < 1e-6 ) {
    return 1.0;
  }
  return (po-pe)/(1.0-pe);

}



// [[Rcpp::export]]
double C_AdjRand(NumericVector x, NumericVector y){
  double eps = 1e-3;
  R_xlen_t n = x.size();
  if (x.size() != y.size()) Rcpp::stop("X and Y must have the same length.");
  NumericMatrix xy(n, 2);
  xy.column(0) = x;
  xy.column(1) = y;
  std::map<double, double> countsx;
  std::map<double, double> countsy;
  std::map<std::vector<double>, double> count_rows;
  countsx.clear();
  countsy.clear();
  count_rows.clear();
  NumericVector::iterator x_i, y_i;
  R_xlen_t xy_i = 0;
  for (x_i = x.begin(),  y_i = y.begin(), xy_i = 0;
       x_i != x.end() && y_i != y.end(), xy_i != n; ++x_i, ++y_i, ++xy_i) {
    countsx[ *x_i ]++;
    countsy[ *y_i ]++;
    NumericVector a = xy.row(xy_i);
    std::vector<double> b = Rcpp::as< std::vector<double> >(a);

    // Add to map
    count_rows[ b ] += 1.0;
  }

  double ai = 0.0;
  double bi = 0.0;
  double nij = 0.0;

  for (std::map<double,double>::iterator it = countsx.begin(); it != countsx.end(); ++it)  {
    double tmp = it->second;
    ai += (tmp) * (tmp - 1.0)/2.0;
  }
  for (std::map<double,double>::iterator it = countsy.begin(); it != countsy.end(); ++it)  {
    double tmp = it->second;
    bi += (tmp) * (tmp - 1.0)/2.0;
  }
  for (std::map<std::vector<double>, double>::iterator it = count_rows.begin(); it != count_rows.end(); ++it)  {
    double tmp = it->second;
    nij += (tmp) * (tmp - 1.0)/2.0;
  }

  return (nij - ai * bi / (1.0 * n * (n-1.0)/2) + eps) / (.5 * (ai + bi) - ai * bi / (1.0 * n * (n-1.0)/2) + eps);

}
