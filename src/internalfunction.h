#ifndef INTERNALFUNCTION_H_
#define INTERNALFUNCTION_H_


 double C_gini(Rcpp::NumericVector x);
 double C_ginicorr(Rcpp::NumericVector x, double k);
 double C_cfunc(Rcpp::NumericVector x, Rcpp::NumericVector y, double c, double k, bool sqrtflag);
 double C_meansfunc(Rcpp::NumericVector x, Rcpp::NumericVector y, double c);
 double C_Cohen(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double C_AdjRand(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double C_Rand(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double C_sqrtginicorr(Rcpp::NumericVector x, double k);


#endif /* INTERNALFUNCTION_H_ */
