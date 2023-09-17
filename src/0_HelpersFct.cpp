
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


////////////////////////////////////////////////////////////////
/////// Generate patients
////////////////////////////////////////////////////////////////


// Generate matrix of patients' outcomes according to multinomial distribution (equivalent to method 2 of package multibrasBOP2)
// [[Rcpp::export]]
IntegerMatrix GenPts(int NSim,
                     IntegerVector AnaInter,
                     NumericMatrix Probs,
                     double Seed = 1024) {

  int NTtt = Probs.nrow();
  int NbAna = AnaInter.length();
  IntegerVector CumAnaInter = cumsum(AnaInter);
  IntegerMatrix MatPts (NSim * NTtt * NbAna, 10);
  IntegerVector ResOutcomes (2);
  set_seed(Seed);

  for (int t = 0; t < NTtt; t ++) { // Loop over the treatments

    for (int sim = 0; sim < NSim; sim ++) { // Loop over the different simulated trials for each treatment

      IntegerVector Outcomes (4);

      for (int ana = 0; ana < NbAna; ana ++) { // Lastly, loop over interim analyses

        Outcomes = Outcomes + rmultinom_rcpp(1, AnaInter(ana), Probs(t, _));
        ResOutcomes[0] = Outcomes[0] + Outcomes[1];
        ResOutcomes[1] = Outcomes[1] + Outcomes[3];
        MatPts(t * NSim * NbAna + sim * NbAna + ana, _) = IntegerVector::create(t, sim, ana, CumAnaInter(ana),
               Outcomes[0], Outcomes[1], Outcomes[2], Outcomes[3], ResOutcomes[0], ResOutcomes[1]);

      }

    }

  }

  return MatPts;

}
