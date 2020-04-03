#ifndef INTERNALFUNCTION_H_
#define INTERNALFUNCTION_H_


 double c_gini(Rcpp::NumericVector x);
 double c_ginicorr(Rcpp::NumericVector x, double k);
 double c_cfunc(Rcpp::NumericVector x, Rcpp::NumericVector y, double c, double k, bool sqrtflag);
 double c_meansfunc(Rcpp::NumericVector x, Rcpp::NumericVector y, double c);
 double c_cohen(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_adj_rand(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_rand(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_nmi(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_ami(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_sqrtginicorr(Rcpp::NumericVector x, double k);


#endif /* INTERNALFUNCTION_H_ */
