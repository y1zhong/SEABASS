#ifndef PGVEC_H
#define PGVEC_H

#include <Rcpp.h>
Rcpp::NumericVector rcpp_pgdraw(Rcpp::NumericVector b, Rcpp::NumericVector c);
double samplepg(double z);
double exprnd(double mu);
double aterm(int n, double x, double t);
double randinvg(double mu);
double truncgamma();
double tinvgauss(double z, double t);

#endif
