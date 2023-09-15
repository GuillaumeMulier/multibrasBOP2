
#include <Rcpp.h>
using namespace Rcpp;


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
