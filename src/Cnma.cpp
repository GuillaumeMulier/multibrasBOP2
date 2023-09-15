
#include <Rcpp.h>
using namespace Rcpp;

////////////////////////////////////////////////////////////////
/////// Cnma threshold
////////////////////////////////////////////////////////////////


// Generate list of matrix of decisions for Cnma
// [[Rcpp::export]]
List GetDecisionCnma(IntegerVector AnasInters,
                     int NBras,
                     double LambdaH,
                     double GammaH,
                     double LambdaMono,
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
    LogicalMatrix TempMatEff ((npts + 1) * NBras, npts + 1); // Extend the number of lines as if we binded rows for different number of active arms
    LogicalMatrix TempMatNoTox ((npts + 1) * NBras, npts + 1);

    for (int a = 1; a < (NBras + 1); a++) { // Loop over different cases of remaining active arms

      double aa = a;
      double eta = 1 + NBras - aa;
      double seuil = (1 - LambdaH / eta) * pow(nnpts / MaxPts, GammaH);
      if ((i + 1) == NbAna) {
        seuil = std::max(seuil, LambdaMono);
      }

      // Loop over possible values for responses in control and treatment group
      for (int nevt0 = 0; nevt0 < (npts + 1); nevt0++) {

        for (int nevt = 0; nevt < (npts + 1); nevt++) {

          double err_est;
          int err_code;
          PostProba DistEff(Delta[0], PriorEff, 1 - PriorEff, PriorEff, 1 - PriorEff, npts, npts, nevt0, nevt);
          double PostProbEff = integrate(DistEff, 0, 1 - Delta[0], err_est, err_code);
          TempMatEff(nevt0 + (a - 1) * (npts + 1), nevt) = PostProbEff < seuil;
          PostProba DistTox(Delta[1], PriorNoTox, 1 - PriorNoTox, PriorNoTox, 1 - PriorNoTox, npts, npts, nevt0, nevt);
          double PostProbTox = integrate(DistTox, 0, 1 - Delta[1], err_est, err_code);
          TempMatNoTox(nevt0 + (a - 1) * (npts + 1), nevt) = PostProbTox < seuil;

        }

      }

    }

    // This time, outside of inner loop because we want to have a list element per interim analysis
    ValeursEff[i] = TempMatEff;
    ValeursTox[i] = TempMatNoTox;

  }

  return List::create(ValeursEff, ValeursTox);

}


// Given a threshold and a set of patients, what are the operating characteristics of the threshold
// [[Rcpp::export]]
NumericMatrix GetOCCnma(IntegerMatrix MatPts,
                        NumericVector Prior,
                        NumericVector Phi,
                        IntegerVector AnaInter,
                        IntegerVector AnaEff,
                        IntegerVector AnaTox,
                        int NbAna,
                        double LambdaH,
                        double GammaH,
                        double LambdaMono,
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
  IntegerVector AnaInterCum = cumsum(AnaInter);

  if (Control) {

    List ListeDeciEff = Decisions[0];
    List ListeDeciTox = Decisions[1];

    for (int s = 0; s < NSimu; s++) { // Loop over all simulated trials

      // We are inside a trial. Contrary to non controlled case, we are gonna analyse trial by trial to reduce the number of accesses to the control group
      LogicalVector ArretsTtt (NTtt - 1);

      for (int ana = 0; ana < NbAna; ana++) { // Loop over interim analyses first to get the control group once

        // Data of control group
        IntegerVector TempTrialCont = MatPts(s * NbAna + ana, _);
        int eff0 = TempTrialCont[8];
        int notox0 = TempTrialCont[9];

        // Decision matrix
        LogicalMatrix DecisionEff = ListeDeciEff[ana];
        LogicalMatrix DecisionTox = ListeDeciTox[ana];
        int NActifs = sum(1 - ArretsTtt);
        int NPtsAna = AnaInterCum[ana];

        for (int t = 1; t < NTtt; t++) { // Loop over treatments and compare them with control group

          if (ArretsTtt[t - 1]) continue; // Go to next treatment if this arm is terminated

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
            arret_eff = DecisionEff(eff0 + (NActifs - 1) * (NPtsAna + 1), eff);
            if (arret_eff) arret = true;
          }
          if (ana_tox[0]) {
            arret_tox = DecisionTox(notox0 + (NActifs - 1) * (NPtsAna + 1), notox);
            if (arret_tox) arret = true;
          }

          if ((ana < (NbAna - 1)) & arret) {
            ArretPrec[t]++;
            NPts[t] = NPts[t] + n;
            ArretsTtt[t - 1] = true; // Mark that this treatment arm is terminated
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

      // Here we are inside a simulated trial
      LogicalVector BrasArretes(NTtt);

      for (int ana = 0; ana < NbAna; ana++) { // Loop over interim analyses

        double NBrasAct = sum(!BrasArretes);
        double Eta = NTtt + 1 - NBrasAct;

        for (int t = 0; t < NTtt; t++) { // Loop over all treatments

          if (BrasArretes[t]) continue; // If the treatment is already stopped, continue to next treatment and don't compute anything

          bool arret = false;
          bool arret_eff = false;
          bool arret_tox = false;

          // Data for trial analysis
          IntegerVector TempTrial = MatPts(t * NSimu * NbAna + s * NbAna + ana, _); // Gain time to subset like that by "guessing" the row number
          int n = TempTrial[3];
          double nn = TempTrial[3];
          int eff = TempTrial[8];
          int notox = TempTrial[9];
          double seuil = 1 - (1 - LambdaH / Eta) * pow(nn / MaxPtsTtt, GammaH);
          if ((nn == MaxPtsTtt) & (seuil > 1 - LambdaMono)) { // Correction with mono-arm threshold at final interim analysis if needed
            seuil = 1 - LambdaMono;
          }
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
            BrasArretes[t] = true;
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
NumericVector DeterCnma(double Fwer,
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
                        NumericVector LambdaHSeq,
                        NumericVector GammaHSeq,
                        NumericVector LambdaMonoSeq,
                        NumericVector GammaMonoSeq,
                        bool Control = false,
                        double Seed = 1024) {

  // Global variables used after
  int LongueurLambda = LambdaHSeq.length();
  int LongueurGamma = GammaHSeq.length();
  int LongueurLMono = LambdaMonoSeq.length();
  int LongueurGMono = GammaMonoSeq.length();
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
  NumericMatrix MatriceCarac(LongueurGamma * 5, 6);
  int DebutGamma = 0;
  NumericVector TempVec = rev(LambdaHSeq); // Reverse order for algorithm to work
  LambdaHSeq = TempVec;

  if (Control) { // With control group

    IntegerMatrix TabH0 = GenPts(NSim, AnaInter, MatH0, Seed);
    IntegerMatrix TabLFC = GenPts(NSim, AnaInter, MatLFC, Seed);
    double LMono;
    double GMono;

    if ((LongueurLMono == 1) & (LongueurGMono == 1)) {

      // Get the mono-arm threshold if declared
      LMono = LambdaMonoSeq[0];
      GMono = GammaMonoSeq[0];

    } else {

      // Optimize monoarm threshold if not supplied
      NumericVector ResMono = DeterCnm(Fwer, 1, NSim, AnaInter, AnaEff, AnaTox, Prior, PN, PA, Phi, Delta,
                                       LambdaMonoSeq, GammaMonoSeq, true, Seed);
      LMono = ResMono[0];
      GMono = ResMono[1];

    }

    for (int l = 0; l < LongueurLambda; l++) { // Loop over possible values of Lambda

      for (int g = DebutGamma; g < LongueurGamma; g++) { // Loop over possible values of Gamma

        if (DebutGamma == (LongueurGamma - 1)) break;
        List ListeDeci = GetDecisionCnma(AnaInterCum, NBras, LambdaHSeq[l], GammaHSeq[g], LMono, Prior, Delta);
        NumericMatrix ResH0 = GetOCCnma(TabH0, Prior, Phi, AnaInter, AnaEff, AnaTox, NbAna, LambdaHSeq[l], GammaHSeq[g], LMono, ListeDeci, true, Delta);
        if (ResH0(0, 0) > Fwer) break;
        NumericMatrix ResLFC = GetOCCnma(TabLFC, Prior, Phi, AnaInter, AnaEff, AnaTox, NbAna, LambdaHSeq[l], GammaHSeq[g], LMono, ListeDeci, true, Delta);
        MatriceCarac(Ligne, _) = NumericVector::create(LambdaHSeq[l], GammaHSeq[g], LMono, GMono, ResH0(0, 0), ResLFC(3, 0));
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
    double LMono;
    double GMono;

    if ((LongueurLMono == 1) & (LongueurGMono == 1)) {

      // Get the mono-arm threshold if declared
      LMono = LambdaMonoSeq[0];
      GMono = GammaMonoSeq[0];

    } else {

      // Optimize monoarm threshold if not supplied
      NumericVector ResMono = DeterCnm(Fwer, 1, NSim, AnaInter, AnaEff, AnaTox, Prior, PN, PA, Phi, Delta,
                                       LambdaMonoSeq, GammaMonoSeq, false, Seed);
      LMono = ResMono[0];
      GMono = ResMono[1];

    }

    // Optimize Holm threshold
    for (int l = 0; l < LongueurLambda; l++) { // Loop over possible values of Lambda

      for (int g = DebutGamma; g < LongueurGamma; g++) { // Loop over possible values of Gamma

        if (DebutGamma == (LongueurGamma - 1)) break;
        NumericMatrix ResH0 = GetOCCnma(TabH0, Prior, Phi, AnaInter, AnaEff, AnaTox, NbAna, LambdaHSeq[l], GammaHSeq[g], LMono, ListeDeci, false);
        if (ResH0(0, 0) > Fwer) break;
        NumericMatrix ResLFC = GetOCCnma(TabLFC, Prior, Phi, AnaInter, AnaEff, AnaTox, NbAna, LambdaHSeq[l], GammaHSeq[g], LMono, ListeDeci, false);
        MatriceCarac(Ligne, _) = NumericVector::create(LambdaHSeq[l], GammaHSeq[g], LMono, GMono, ResH0(0, 0), ResLFC(3, 0));
        Ligne++;
        DebutGamma = g;

      }

      DebutGamma = DebutGamma; // Update the value outside of the inside loop to have access in outside loop

    }

  }

  // Return the first couple with highest power
  NumericVector Puissance = MatriceCarac(_, 5);
  return MatriceCarac(which_max(Puissance), _);

}
