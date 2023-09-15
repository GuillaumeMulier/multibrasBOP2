
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;


////////////////////////////////////////////////////////////////
/////// Helper functions
////////////////////////////////////////////////////////////////


// Set the seed like in R
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}


// Generate one sample out of multinomial distribution
// [[Rcpp::export]]
IntegerVector rmultinom_1(unsigned int size, NumericVector probs, unsigned int N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}


// Generate n samples out of multinomial distribution
// [[Rcpp::export]]
IntegerMatrix rmultinom_rcpp(unsigned int n, unsigned int size, NumericVector probs) {
  unsigned int N = probs.length();
  IntegerMatrix sim(N, n);
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int s = 0; s < size; s ++) {
      sim(_,i) = sim(_,i) + rmultinom_1(1, probs, N);
    }
  }
  return sim;
}


// Posterior probability for comparison of 2 beta distributions
// [[Rcpp::export]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
class PostProba: public Func {
private:
  double delta;
  double alpha0;
  double beta0;
  double alpha1;
  double beta1;
  int n0;
  int n1;
  int nevt0;
  int nevt1;

public:
  PostProba(double delta_, double alpha0_, double beta0_, double alpha1_, double beta1_, int n0_, int n1_, int nevt0_, int nevt1_) : delta(delta_), alpha0(alpha0_), beta0(beta0_), alpha1(alpha1_), beta1(beta1_), n0(n0_), n1(n1_), nevt0(nevt0_), nevt1(nevt1_) {}

  double operator()(const double& x) const {
    return R::pbeta(x + delta, alpha1 + nevt1, beta1 + n1 - nevt1, false, false) * R::dbeta(x, alpha0 + nevt0, beta0 + n0 - nevt0, false);
  }
};
