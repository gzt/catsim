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
 double c_jaccard(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_dice(Rcpp::NumericVector x, Rcpp::NumericVector y);
 double c_hamming(Rcpp::NumericVector x, Rcpp::NumericVector y);
 Rcpp::NumericVector c_catssim_2d(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y, Rcpp::IntegerVector window, std::string method, double c1, double c2, bool sqrtgini);

#endif /* INTERNALFUNCTION_H_ */
