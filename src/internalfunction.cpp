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

// Helper functions for Jaccard, Dice, and Hamming
// [[Rcpp::export]]
double c_jaccard(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  R_xlen_t n = x.size();

  if (x.size() != y.size()) {
    Rcpp::stop("X and Y must have the same length.");
  }

  double sum_xy_or = 0.0;
  double sum_xy_and = 0.0;

  for (R_xlen_t i = 0; i < n; i++) {
    const double xi = x[i];
    const double yi = y[i];

    const bool x_na = Rcpp::traits::is_na<REALSXP>(xi);
    const bool y_na = Rcpp::traits::is_na<REALSXP>(yi);

    // R's `x | y`, with na.rm = TRUE:
    //
    // NA | TRUE  -> TRUE
    // NA | FALSE -> NA
    // NA | NA    -> NA
    //
    // Therefore OR is TRUE whenever either non-NA value is nonzero.
    bool or_true =
      (!x_na && xi != 0.0) ||
      (!y_na && yi != 0.0);

    if (or_true) {
      sum_xy_or += 1.0;
    }

    // R's `x & y`, with na.rm = TRUE:
    //
    // FALSE & NA -> FALSE
    // TRUE  & NA -> NA
    // NA    & NA -> NA
    //
    // Therefore AND is TRUE only when both values are known and nonzero.
    bool and_true =
      (!x_na && xi != 0.0) &&
      (!y_na && yi != 0.0);

    if (and_true) {
      sum_xy_and += 1.0;
    }
  }

  if (sum_xy_or == 0.0) {
    return NA_REAL;
  }

  return sum_xy_and / sum_xy_or;
}

// [[Rcpp::export]]
double c_dice(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  double jacc = c_jaccard(x, y);
  if (Rcpp::traits::is_nan<REALSXP>(jacc)) {
    return NA_REAL;
  }
  return 2.0 * jacc / (1.0 + jacc);
}

// [[Rcpp::export]]
double c_hamming(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  R_xlen_t n = x.size();
  if (x.size() != y.size()) Rcpp::stop("X and Y must have the same length.");

  double n_valid = 0.0;
  double n_agree = 0.0;

  for (R_xlen_t i = 0; i < n; i++) {
    if (Rcpp::traits::is_nan<REALSXP>(x[i]) || Rcpp::traits::is_nan<REALSXP>(y[i])) {
      continue;
    }
    n_valid += 1.0;
    if (x[i] == y[i]) {
      n_agree += 1.0;
    }
  }

  if (n_valid < 1e-9) {
    return NA_REAL;
  }

  return n_agree / n_valid;
}

// [[Rcpp::export]]
Rcpp::NumericVector c_catssim_2d(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y,
                                 Rcpp::IntegerVector window, std::string method,
                                 double c1, double c2, bool sqrtgini) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  int win_rows = window[0];
  int win_cols = window[1];

  int out_nrow = nrow - win_rows + 1;
  int out_ncol = ncol - win_cols + 1;

  Rcpp::NumericMatrix resultmatrix(out_nrow, out_ncol * 3);

  double (*method_func)(NumericVector, NumericVector);
  if (method == "Cohen" || method == "cohen" || method == "C" || method == "c" ||
      method == "kappa" || method == "Kappa") {
    method_func = c_cohen;
  } else if (method == "AdjRand" || method == "adjrand" || method == "Adj" ||
             method == "adj" || method == "a" || method == "A" ||
             method == "ARI" || method == "ari") {
    method_func = c_adj_rand;
  } else if (method == "Rand" || method == "rand" || method == "r" ||
             method == "R") {
    method_func = c_rand;
  } else if (method == "NMI" || method == "MI" || method == "mutual" ||
             method == "information" || method == "nmi" || method == "mi") {
    method_func = c_nmi;
  } else if (method == "AMI" || method == "ami") {
    method_func = c_ami;
  } else if (method == "Jaccard" || method == "jaccard" || method == "j" ||
             method == "J") {
    method_func = c_jaccard;
  } else if (method == "Dice" || method == "dice" || method == "D" ||
             method == "d") {
    method_func = c_dice;
  } else if (method == "Accuracy" || method == "accuracy" || method == "acc" ||
             method == "Hamming" || method == "hamming" || method == "H" ||
             method == "h") {
    method_func = c_hamming;
  } else {
    Rcpp::stop("Error: invalid method");
  }

  std::set<double> unique_vals;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      unique_vals.insert(x(i, j));
      unique_vals.insert(y(i, j));
    }
  }
  double k = unique_vals.size();

  for (int i = 0; i < out_nrow; i++) {
    for (int j = 0; j < out_ncol; j++) {
      NumericVector subx(win_rows * win_cols);
      NumericVector suby(win_rows * win_cols);

      int idx = 0;
      for (int wi = 0; wi < win_rows; wi++) {
        for (int wj = 0; wj < win_cols; wj++) {
          subx[idx] = x(i + wi, j + wj);
          suby[idx] = y(i + wi, j + wj);
          idx++;
        }
      }

      double comp1 = c_meansfunc(subx, suby, c1);
      double comp2 = c_cfunc(subx, suby, c2, k, sqrtgini);
      double comp3 = method_func(subx, suby);

      resultmatrix(i, j) = comp1;
      resultmatrix(i, j + out_ncol) = comp2;
      resultmatrix(i, j + 2 * out_ncol) = comp3;
    }
  }

  for (int i = 0; i < out_nrow; i++) {
    for (int j = 0; j < out_ncol; j++) {
      if (resultmatrix(i, j) < 0.0) resultmatrix(i, j) = 0.0;
      if (resultmatrix(i, j + out_ncol) < 0.0) resultmatrix(i, j + out_ncol) = 0.0;
      if (resultmatrix(i, j + 2 * out_ncol) < 0.0) resultmatrix(i, j + 2 * out_ncol) = 0.0;
    }
  }

  Rcpp::NumericVector col_means(3);
  for (int comp = 0; comp < 3; comp++) {
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < out_nrow; i++) {
      for (int j = 0; j < out_ncol; j++) {
        double val = resultmatrix(i, j + comp * out_ncol);
        if (!Rcpp::traits::is_nan<REALSXP>(val)) {
          sum += val;
          count++;
        }
      }
    }
    col_means[comp] = (count > 0) ? sum / count : NA_REAL;
  }

  return col_means;
}
