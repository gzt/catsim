#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double c_gini(NumericVector x) {
  std::map<double, double> counts;
  R_xlen_t n = x.size();
  // NumericVector::iterator i;
  for (NumericVector::iterator i = x.begin(); i != x.end(); ++i) {
    counts[*i]++;
  }

  double sqsum = 0.0;
  for (std::map<double, double>::iterator it = counts.begin();
       it != counts.end(); ++it) {
    sqsum += 1.0 * (it->second) * (it->second);
  }
  return (1.0 - sqsum / (1.0 * n * n));
}

// [[Rcpp::export]]
double c_ginicorr(NumericVector x, double k) {
  double eps = 1e-5;
  if (std::abs(k - 1.0) < eps) return 1.0;

  return c_gini(x) / (1.0 - 1.0 / k);
}

double c_sqrtginicorr(NumericVector x, double k) {
  double eps = 1e-5;
  if (std::abs(k - 1.0) < eps) return 1.0;

  return (1 - sqrt(1 - c_gini(x))) / (1 - 1.0 / k);
}

// [[Rcpp::export]]
double c_cfunc(NumericVector x, NumericVector y, double c, double k,
               bool sqrtflag) {
  double varx, vary;
  if (sqrtflag) {
    varx = c_sqrtginicorr(x, k);
    vary = c_sqrtginicorr(y, k);
  } else {
    varx = c_ginicorr(x, k);
    vary = c_ginicorr(y, k);
  }

  return (2 * sqrt(varx * vary) + c) / (varx + vary + c);
}

// [[Rcpp::export]]
double c_meansfunc(NumericVector x, NumericVector y, double c) {
  // R_xlen_t n = x.size();
  if (x.size() != y.size()) Rcpp::stop("X and Y must have the same length.");
  std::map<double, double> countsx;
  std::map<double, double> countsy;
  NumericVector::iterator x_i, y_i;
  for (x_i = x.begin(), y_i = y.begin(); x_i != x.end() && y_i != y.end();
       ++x_i, ++y_i) {
    countsx[*x_i]++;
    countsy[*y_i]++;
  }

  double sqsum = 0.0;
  for (std::map<double, double>::iterator it = countsx.begin();
       it != countsx.end(); ++it) {
    sqsum += 1.0 * (it->second) * (it->second);
  }
  for (std::map<double, double>::iterator it = countsy.begin();
       it != countsy.end(); ++it) {
    sqsum += 1.0 * (it->second) * (it->second);
  }

  double xysum = 0.0;
  std::map<double, double>::iterator il = countsx.begin();
  std::map<double, double>::iterator ir = countsy.begin();
  while (il != countsx.end() && ir != countsy.end()) {
    if (il->first < ir->first)
      ++il;
    else if (ir->first < il->first)
      ++ir;
    else {
      xysum += (il->second) * (ir->second);
      ++il;
      ++ir;
    }
  }

  return (2 * (xysum) + c) / (sqsum + c);
}

// [[Rcpp::export]]
double c_cohen(NumericVector x, NumericVector y) {
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
  for (x_i = x.begin(), y_i = y.begin(); x_i != x.end() && y_i != y.end();
       ++x_i, ++y_i) {
    countsx[*x_i]++;
    countsy[*y_i]++;
    if (*x_i == *y_i) {
      countsxy[*x_i]++;
    }
  }

  double xxyysum = 0.0;
  std::map<double, double>::iterator il = countsx.begin();
  std::map<double, double>::iterator ir = countsy.begin();
  while (il != countsx.end() && ir != countsy.end()) {
    if (il->first < ir->first)
      ++il;
    else if (ir->first < il->first)
      ++ir;
    else {
      xxyysum += 1.00 * (il->second) * (ir->second);
      ++il;
      ++ir;
    }
  }
  double xysum = 0.0;

  for (std::map<double, double>::iterator it = countsxy.begin();
       it != countsxy.end(); ++it) {
    xysum += 1.0 * (it->second);
  }

  double pe = xxyysum / (1.0 * n * n);
  double po = xysum / (1.0 * n);

  if ((1.0 - pe) < 1e-6) {
    return 1.0;
  }
  return (po - pe) / (1.0 - pe);
}

Rcpp::NumericVector c_randRaw(NumericVector x, NumericVector y) {
  Rcpp::NumericVector resultvector(3);
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
  for (x_i = x.begin(), y_i = y.begin(), xy_i = 0;
       x_i != x.end() && y_i != y.end() && xy_i != n; ++x_i, ++y_i, ++xy_i) {
    countsx[*x_i]++;
    countsy[*y_i]++;
    NumericVector a = xy.row(xy_i);
    std::vector<double> b = Rcpp::as<std::vector<double> >(a);

    // Add to map
    count_rows[b] += 1.0;
  }

  double ai = 0.0;
  double bi = 0.0;
  double nij = 0.0;

  for (std::map<double, double>::iterator it = countsx.begin();
       it != countsx.end(); ++it) {
    ai += (it->second) * ((it->second) - 1.0) / 2.0;
  }
  for (std::map<double, double>::iterator it = countsy.begin();
       it != countsy.end(); ++it) {
    bi += (it->second) * ((it->second) - 1.0) / 2.0;
  }
  for (std::map<std::vector<double>, double>::iterator it = count_rows.begin();
       it != count_rows.end(); ++it) {
    nij += (it->second) * ((it->second) - 1.0) / 2.0;
  }

  resultvector[0] = nij;
  resultvector[1] = ai;
  resultvector[2] = bi;

  return resultvector;
}

// [[Rcpp::export]]
double c_adj_rand(NumericVector x, NumericVector y) {
  double eps = 0.0;
  R_xlen_t n = x.size();

  double ai = 0.0;
  double bi = 0.0;
  double nij = 0.0;

  NumericVector resultvector(3);
  resultvector = c_randRaw(x, y);
  nij = resultvector[0];
  ai = resultvector[1];
  bi = resultvector[2];
  if ((.5 * (ai + bi) - ai * bi / (1.0 * n * (n - 1.0) / 2)) < 1e-9) eps = 1e-9;

  return (nij - ai * bi / (1.0 * n * (n - 1.0) / 2) + eps) /
         (.5 * (ai + bi) - ai * bi / (1.0 * n * (n - 1.0) / 2) + eps);
}

// [[Rcpp::export]]
double c_rand(NumericVector x, NumericVector y) {
  R_xlen_t n = x.size();

  double ai = 0.0;
  double bi = 0.0;
  double nij = 0.0;

  NumericVector resultvector(3);
  resultvector = c_randRaw(x, y);
  nij = resultvector[0];
  ai = resultvector[1];
  bi = resultvector[2];

  return (n * (n - 1.0) / 2.0 + 2.0 * nij - ai - bi) / (n * (n - 1.0) / 2.0);
}

// [[Rcpp::export]]
double c_nmi(NumericVector x, NumericVector y) {
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
  for (x_i = x.begin(), y_i = y.begin(), xy_i = 0;
       x_i != x.end() && y_i != y.end() && xy_i != n; ++x_i, ++y_i, ++xy_i) {
    countsx[*x_i]++;
    countsy[*y_i]++;
    NumericVector a = xy.row(xy_i);
    std::vector<double> b = Rcpp::as<std::vector<double> >(a);

    // Add to map
    count_rows[b] += 1.0;
  }

  double IXY = 0.0;
  double HX = 0.0;
  double HY = 0.0;

  for (std::map<double, double>::iterator it = countsx.begin();
       it != countsx.end(); ++it) {
    double tmp = it->second;
    HX += -tmp / n * log(tmp / n);
  }
  for (std::map<double, double>::iterator it = countsy.begin();
       it != countsy.end(); ++it) {
    double tmp = it->second;
    HY += -tmp / n * log(tmp / n);
  }
  for (std::map<std::vector<double>, double>::iterator it = count_rows.begin();
       it != count_rows.end(); ++it) {
    double tmp = it->second;
    std::vector<double> tmpvec = it->first;
    double xn = countsx[tmpvec[0]];
    double yn = countsy[tmpvec[1]];
    IXY += tmp / n * log(tmp * n / (xn * yn));
  }
  if (HX + HY < 1e-9) return 0.0;

  return 2 * IXY / (HX + HY);
}

double lfact(int x) { return lgamma(x + 1.0); }

double hypergeomfunc(double ai, double bj, R_xlen_t N) {
  double tmp = 0.0;
  int maxval = (ai < bj ? round(ai) : round(bj));  // min(ai, bj);
  int minval =
      (1 > (ai + bj - N) ? 1 : round(ai + bj - N));  // max(1,ai + bj - N);
  for (int nij = minval; nij <= maxval; ++nij) {
    if (ai > 0 && bj > 0) {
      tmp += nij / N * log(N * nij / (ai * bj)) *
             exp(lfact(ai) + lfact(bj) + lfact(N - ai) + lfact(N - bj) -
                 lfact(N) - lfact(nij) - lfact(ai - nij) - lfact(bj - nij) -
                 lfact(N - ai - bj + nij));
    }
  }
  return tmp;
}

// [[Rcpp::export]]
double c_ami(NumericVector x, NumericVector y) {
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
  for (x_i = x.begin(), y_i = y.begin(), xy_i = 0;
       x_i != x.end() && y_i != y.end() && xy_i != n; ++x_i, ++y_i, ++xy_i) {
    countsx[*x_i]++;
    countsy[*y_i]++;
    NumericVector a = xy.row(xy_i);
    std::vector<double> b = Rcpp::as<std::vector<double> >(a);

    // Add to map
    count_rows[b] += 1.0;
  }

  double IXY = 0.0;
  double HX = 0.0;
  double HY = 0.0;
  double EMI = 0.0;
  for (std::map<double, double>::iterator it = countsx.begin();
       it != countsx.end(); ++it) {
    double tmp = it->second;
    HX += -tmp / n * log(tmp / n);
  }
  for (std::map<double, double>::iterator it = countsy.begin();
       it != countsy.end(); ++it) {
    double tmp = it->second;
    HY += -tmp / n * log(tmp / n);
  }
  for (std::map<std::vector<double>, double>::iterator it = count_rows.begin();
       it != count_rows.end(); ++it) {
    double tmp = it->second;
    double xn = countsx[(it->first)[0]];
    double yn = countsy[(it->first)[1]];
    IXY += tmp / n * log(tmp * n / (xn * yn));
    EMI += hypergeomfunc(xn, yn, n);
  }
  return (IXY - EMI) / ((HX > HY ? HX : HY) - EMI);
}
