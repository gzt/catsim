#ifndef INTERNALFUNCTION_H_
#define INTERNALFUNCTION_H_


 double c_gini(Rcpp::NumericVector x);
 double c_ginicorr(SEXP x, double k);
 double c_cfunc(SEXP x, SEXP y, double c, double k, bool sqrtflag);
 double c_meansfunc(SEXP x, SEXP y, double c);
 double c_cohen(SEXP x, SEXP y);
 double c_adj_rand(SEXP x, SEXP y);
 double c_rand(SEXP x, SEXP y);
 double c_nmi(SEXP x, SEXP y);
 double c_ami(SEXP x, SEXP y);
 double c_sqrtginicorr(SEXP x, double k);


#endif /* INTERNALFUNCTION_H_ */
