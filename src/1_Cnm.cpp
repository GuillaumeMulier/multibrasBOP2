
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;

// Declare functions

void set_seed(double seed);

IntegerVector rmultinom_1(unsigned int size, NumericVector probs, unsigned int N);

IntegerMatrix rmultinom_rcpp(unsigned int n, unsigned int size, NumericVector probs);

IntegerMatrix GenPts(int NSim,
                     IntegerVector AnaInter,
                     NumericMatrix Probs,
                     double Seed = 1024);

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


////////////////////////////////////////////////////////////////
/////// Cnm threshold
////////////////////////////////////////////////////////////////


// Generate list of matrix of decisions for Cnm
// [[Rcpp::export]]
List GetDecision(IntegerVector AnasInters,
                 double LambdaM,
                 double GammaM,
                 NumericVector Prior,
                 NumericVector Delta) {

  double PriorEff = Prior[0] + Prior[1];
  double PriorNoTox = Prior[1] + Prior[3];
  int NbAna = AnasInters.length();
  double MaxPts = max(AnasInters);
  List ValeursEff(NbAna);
  List ValeursTox(NbAna);

  // Limit values for efficacy and toxicity
  for (int i = 0; i < NbAna; i++) {

    int npts = AnasInters[i];
    double nnpts = AnasInters[i];
    LogicalMatrix TempMatEff (npts + 1, npts + 1);
    LogicalMatrix TempMatNoTox (npts + 1, npts + 1);
    double seuil = LambdaM * pow(nnpts / MaxPts, GammaM);

    // Loop over possible values for responses in control and treatment group
    for (int nevt0 = 0; nevt0 < npts + 1; nevt0++) {

      for (int nevt = 0; nevt < npts + 1; nevt++) {

        double err_est;
        int err_code;
        PostProba DistEff(Delta[0], PriorEff, 1 - PriorEff, PriorEff, 1 - PriorEff, npts, npts, nevt0, nevt);
        double PostProbEff = integrate(DistEff, 0, 1 - Delta[0], err_est, err_code);
        TempMatEff(nevt0, nevt) = PostProbEff < seuil;
        PostProba DistTox(Delta[1], PriorNoTox, 1 - PriorNoTox, PriorNoTox, 1 - PriorNoTox, npts, npts, nevt0, nevt);
        double PostProbTox = integrate(DistTox, 0, 1 - Delta[1], err_est, err_code);
        TempMatNoTox(nevt0, nevt) = PostProbTox < seuil;

      }

    }

    ValeursEff[i] = TempMatEff;
    ValeursTox[i] = TempMatNoTox;

  }

  return List::create(ValeursEff, ValeursTox);

}


// Given a threshold and a set of patients, what are the operating characteristics of the threshold
// [[Rcpp::export]]
NumericMatrix GetOCCnm(IntegerMatrix MatPts,
                       NumericVector Prior,
                       NumericVector Phi,
                       IntegerVector AnaEff,
                       IntegerVector AnaTox,
                       int NbAna,
                       double LambdaM,
                       double GammaM,
                       List Decisions,
                       bool Control = false,
                       NumericVector Delta = 0) {

  double PriorEff = Prior[0] + Prior[1];
  double PriorNoTox = Prior[1] + Prior[3];
  int NSimu = max(MatPts(_, 1)) + 1;
  int NTtt = max(MatPts(_, 0)) + 1;
  NumericMatrix Resultats(4, NTtt);
  NumericVector ArretPrec (NTtt);
  NumericVector NPts (NTtt);
  NumericVector PositiveTrials (NTtt);
  NumericVector PromisingTtt (NSimu);
  double MaxPtsTtt = max(MatPts(_, 3));

  if (Control) {

    List ListeDeciEff = Decisions[0];
    List ListeDeciTox = Decisions[1];

    for (int s = 0; s < NSimu; s++) { // Loop over all simulated trials

      // We are inside a trial. Contrary to non controlled case, we are gonna analyse trial by trial to reduce the number of accesses to the control group
      LogicalVector ArretsTtt (NTtt);

      for (int ana = 0; ana < NbAna; ana++) { // Loop over interim analyses first to get the control group once

        // Data of control group
        IntegerVector TempTrialCont = MatPts(s * NbAna + ana, _);
        int eff0 = TempTrialCont[8];
        int notox0 = TempTrialCont[9];

        // Decision matrix
        LogicalMatrix DecisionEff = ListeDeciEff[ana];
        LogicalMatrix DecisionTox = ListeDeciTox[ana];

        for (int t = 1; t < NTtt; t++) { // Loop over treatments and compare them with control group

          if (ArretsTtt[t]) continue; // Go to next treatment if this arm is terminated

          // Data of experimental arm
          IntegerVector TempTrial = MatPts(t * NSimu * NbAna + s * NbAna + ana, _); // Gain time to subset like that by "guessing" the row number
          int n = TempTrial[3];
          int eff = TempTrial[8];
          int notox = TempTrial[9];

          // Determine if there is an efficacy/toxicity analysis
          IntegerVector npts = IntegerVector::create(n);
          LogicalVector ana_eff = !is_na(match(npts, AnaEff));
          LogicalVector ana_tox = !is_na(match(npts, AnaTox));

          // Results of the interim analysis
          bool arret = false;
          bool arret_eff = false;
          bool arret_tox = false;
          if (ana_eff[0]) {
            arret_eff = DecisionEff(eff0, eff);
            if (arret_eff) arret = true;
          }
          if (ana_tox[0]) {
            arret_tox = DecisionTox(notox0, notox);
            if (arret_tox) arret = true;
          }

          if ((ana < (NbAna - 1)) & arret) {
            ArretPrec[t]++;
            NPts[t] = NPts[t] + n;
            ArretsTtt[t] = true; // Mark that this treatment arm is terminated
          } else if (ana == (NbAna - 1)){
            NPts[t] = NPts[t] + n;
            if (!arret) {
              PositiveTrials[t]++;
              PromisingTtt[s]++;
            }
          }

        }

      }

    }

    // Compute per trial operating characteristics
    double PromisingGlob = mean(PromisingTtt > 0);
    for (int t = 1; t < NTtt; t++) {

      Resultats(0, t) = PromisingGlob;
      Resultats(1, t) = ArretPrec[t] / NSimu;
      Resultats(2, t) = NPts[t] / NSimu;
      Resultats(3, t) = PositiveTrials[t] / NSimu;

    }
    Resultats = Resultats(_, Range(1, NTtt - 1));

  } else {

    for (int s = 0; s < NSimu; s++) { // Loop over all simulated trials

      for (int t = 0; t < NTtt; t++) { // Loop over all treatments

        // We are here inside a trial' and a treatment's arm
        bool arret = false;
        bool arret_eff = false;
        bool arret_tox = false;

        for (int ana = 0; ana < NbAna; ana++) {

          // Data for trial analysis
          IntegerVector TempTrial = MatPts(t * NSimu * NbAna + s * NbAna + ana, _); // Gain time to subset like that by "guessing" the row number
          int n = TempTrial[3];
          double nn = TempTrial[3];
          int eff = TempTrial[8];
          int notox = TempTrial[9];
          double seuil = 1 - LambdaM * pow(nn / MaxPtsTtt, GammaM);
          double PostProbEff = R::pbeta(Phi[0], PriorEff + eff, 1 - PriorEff + n - eff, true, false);
          double PostProbTox = R::pbeta(Phi[1], PriorNoTox + notox, 1 - PriorNoTox + n - notox, true, false);

          // Determine if there is an efficacy/toxicity analysis
          IntegerVector npts = IntegerVector::create(n);
          LogicalVector ana_eff = !is_na(match(npts, AnaEff));
          LogicalVector ana_tox = !is_na(match(npts, AnaTox));

          // Results of the interim analysis
          if (ana_eff[0]) {
            arret_eff = PostProbEff > seuil;
            if (arret_eff) arret = true;
          }
          if (ana_tox[0]) {
            arret_tox = PostProbTox > seuil;
            if (arret_tox) arret = true;
          }
          if ((ana < (NbAna - 1)) & arret) {
            ArretPrec[t]++;
            NPts[t] = NPts[t] + n;
            break;
          } else if (ana == (NbAna - 1)){
            NPts[t] = NPts[t] + n;
            if (!arret) {
              PositiveTrials[t]++;
              PromisingTtt[s]++;
            }
          }

        }

      }

    }

    // Compute per trial operating characteristics
    double PromisingGlob = mean(PromisingTtt > 0);
    for (int t = 0; t < NTtt; t++) {

      Resultats(0, t) = PromisingGlob;
      Resultats(1, t) = ArretPrec[t] / NSimu;
      Resultats(2, t) = NPts[t] / NSimu;
      Resultats(3, t) = PositiveTrials[t] / NSimu;

    }

  }

  return Resultats;

}


// Main function to compute the hyperparameters of the threshold that are optimal in sense of FWER and power
// [[Rcpp::export]]
NumericVector DeterCnm(double Fwer,
                       int NBras,
                       int NSim,
                       IntegerVector AnaInter,
                       IntegerVector AnaEff,
                       IntegerVector AnaTox,
                       NumericVector Prior,
                       NumericVector PN,
                       NumericVector PA,
                       NumericVector Phi,
                       NumericVector Delta,
                       NumericVector LambdaMSeq,
                       NumericVector GammaMSeq,
                       bool Control = false,
                       double Seed = 1024) {

  // Global variables used after
  int LongueurLambda = LambdaMSeq.length();
  int LongueurGamma = GammaMSeq.length();
  int NbAna = AnaInter.length();
  IntegerVector AnaInterCum = cumsum(AnaInter);
  int Ligne = 0;
  NumericMatrix MatH0(NBras + 1, 4);
  NumericMatrix MatLFC(NBras + 1, 4);
  for (int i = 0; i < NBras + 1; i++) {
    if (i == 1) {
      MatH0(i, _) = PN;
      MatLFC(i, _) = PA;
    } else {
      MatH0(i, _) = PN;
      MatLFC(i, _) = PN;
    }
  }
  NumericMatrix MatriceCarac(LongueurGamma * 5, 4);
  int DebutGamma = 0;

  if (Control) { // With control group

    IntegerMatrix TabH0 = GenPts(NSim, AnaInter, MatH0, Seed);
    IntegerMatrix TabLFC = GenPts(NSim, AnaInter, MatLFC, Seed);

    for (int l = 0; l < LongueurLambda; l++) { // Loop over possible values of Lambda

      for (int g = DebutGamma; g < LongueurGamma; g++) { // Loop over possible values of Gamma

        if (DebutGamma == (LongueurGamma - 1)) break;
        List ListeDeci = GetDecision(AnaInterCum, LambdaMSeq[l], GammaMSeq[g], Prior, Delta);
        NumericMatrix ResH0 = GetOCCnm(TabH0, Prior, Phi, AnaEff, AnaTox, NbAna, LambdaMSeq[l], GammaMSeq[g], ListeDeci, true, Delta);
        if (ResH0(0, 0) > Fwer) break;
        NumericMatrix ResLFC = GetOCCnm(TabLFC, Prior, Phi, AnaEff, AnaTox, NbAna, LambdaMSeq[l], GammaMSeq[g], ListeDeci, true, Delta);
        MatriceCarac(Ligne, _) = NumericVector::create(LambdaMSeq[l], GammaMSeq[g], ResH0(0, 0), ResLFC(3, 0));
        Ligne++;
        DebutGamma = g;

      }

      DebutGamma = DebutGamma; // Update the value outside of the inside loop to have access in outside loop

    }

  } else { // Against reference value

    NumericMatrix MatH0Ref = MatH0(Range(1, NBras), _);
    NumericMatrix MatLFCRef = MatLFC(Range(1, NBras), _);
    IntegerMatrix TabH0 = GenPts(NSim, AnaInter, MatH0Ref, Seed);
    IntegerMatrix TabLFC = GenPts(NSim, AnaInter, MatLFCRef, Seed);
    List ListeDeci = GetDecision(1, 1, 1, Prior, Delta); // Dummy list because I don't know how to put default argument of list

    for (int l = 0; l < LongueurLambda; l++) { // Loop over possible values of Lambda

      for (int g = DebutGamma; g < LongueurGamma; g++) { // Loop over possible values of Gamma

        if (DebutGamma == (LongueurGamma - 1)) break;
        NumericMatrix ResH0 = GetOCCnm(TabH0, Prior, Phi, AnaEff, AnaTox, NbAna, LambdaMSeq[l], GammaMSeq[g], ListeDeci, false);
        if (ResH0(0, 0) > Fwer) break;
        NumericMatrix ResLFC = GetOCCnm(TabLFC, Prior, Phi, AnaEff, AnaTox, NbAna, LambdaMSeq[l], GammaMSeq[g], ListeDeci, false);
        MatriceCarac(Ligne, _) = NumericVector::create(LambdaMSeq[l], GammaMSeq[g], ResH0(0, 0), ResLFC(3, 0));
        Ligne++;
        DebutGamma = g;

      }

      DebutGamma = DebutGamma; // Update the value outside of the inside loop to have access in outside loop

    }

  }

  // Return the first couple with highest power
  NumericVector Puissance = MatriceCarac(_, 3);
  return MatriceCarac(which_max(Puissance), _);

}
