#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ddv(NumericMatrix xx, NumericVector v, NumericVector w) {

  Environment base("package:base");
  Function mat_Mult = base["%*%"];

  int dimxm = xx.ncol();

  NumericVector xv;
  NumericMatrix xm1;
  NumericMatrix xm2;
  NumericMatrix xm;

  NumericMatrix tot_xm(dimxm, dimxm);


  for(int i = 0; i < xx.nrow(); ++i) {
    xv = xx( i, _ );
    xv.attr("dim") = Dimension(dimxm, 1);
    xm1 = as<NumericMatrix>(xv);
    xm2 = transpose(xm1);
    xm = mat_Mult(xm1, xm2);
    xm = w(i) * xm / v(i); // added in pre-multiplication by w(i)

   for (int j = 0; j < dimxm; j++) {     // loop over rows
     for (int k = 0; k < dimxm; k++) {   // loop over columns
       tot_xm(j,k) = tot_xm(j,k) + xm(j,k); // elementwise addition
     }
   }

  }

  return tot_xm;

}

// [[Rcpp::export]]
NumericVector dv(NumericMatrix xx, NumericVector adj) {

  int dimxv = xx.ncol();

  NumericVector xv;
  NumericVector tot_xv(dimxv);

  for(int i = 0; i < xx.nrow(); ++i) {

    xv = xx( i, _ );
    xv = xv / adj(i);

    for (int j = 0; j < dimxv; j++) {     // loop over rows
      tot_xv(j) = tot_xv(j) + xv(j); // elementwise addition
    }
  }

  return tot_xv;

}

// [[Rcpp::export]]
NumericVector dv2(NumericMatrix xx, NumericVector adj, NumericVector w) {
  
  int dimxv = xx.ncol();

  NumericVector xv;
  NumericVector tot_xv(dimxv);

  for(int i = 0; i < xx.nrow(); ++i) {
  
    xv = xx( i, _ );
    xv = w(i) * xv / adj(i); // scale by weights
    
    for (int j = 0; j < dimxv; j++) {     // loop over rows
      tot_xv(j) = tot_xv(j) + xv(j); // elementwise addition
    }
  }
  
  return tot_xv;
  
}

// [[Rcpp::export]]
NumericVector dvm(NumericMatrix xx, NumericVector adj, NumericVector w) {

  int dimxv = xx.ncol();

  NumericVector xv;
  NumericVector tot_xv(dimxv);

  for(int i = 0; i < xx.nrow(); ++i) {

    xv = xx( i, _ );
    xv = xv * adj(i) * w(i); // scaled by weights 

    for (int j = 0; j < dimxv; j++) {     // loop over rows
      tot_xv(j) = tot_xv(j) + xv(j); // elementwise addition
    }
  }

  return tot_xv;

}



