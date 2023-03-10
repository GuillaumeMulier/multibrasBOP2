% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_patients_TOP.R
\encoding{UTF-8}
\name{gen_patients_multinomTOP}
\alias{gen_patients_multinomTOP}
\title{Generate a list of trials with time to event for each patient}
\usage{
gen_patients_multinomTOP(
  n_sim,
  ana_inter,
  ana_inter_tox = NULL,
  interpatient,
  loi_gen_patients = "unif",
  max_teff,
  prop_demieff = 0.5,
  max_ttox,
  prop_demitox = 0.5,
  rand_ratio = NULL,
  multinom_cont = NULL,
  multinom_ttt = list(),
  mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
  seed = 1024
)
}
\arguments{
\item{n_sim}{The number of simulated trials.}

\item{ana_inter}{A vector of the supplementray number of patients at each interim analysis.}

\item{ana_inter_tox}{Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
number of patients, set it to NULL.}

\item{interpatient}{Mean time between 2 patients.}

\item{loi_gen_patients}{The law of generation of time to event (uniform, weibull, log-logistic).}

\item{max_teff}{The length of the observation window for efficacy.}

\item{prop_demieff}{The proportion of responses at \code{max_teff}/2.}

\item{max_ttox}{The length of the observation window for toxicity.}

\item{prop_demitox}{The proportion of toxicities at \code{max_ttox}/2.}

\item{rand_ratio}{A vector of the same length of \code{multinom_ttt} or length 1 giving the ratio of patients between control and treatment groups.}

\item{multinom_cont}{Vector of length 4 giving the law of control group (Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox). NULL if no control arm.}

\item{multinom_ttt}{List of vector of length 4 giving the law of treatment groups (Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox). 1 vector if 1 group is accepted.}

\item{mat_beta_xi}{The link matrix between outcomes of the multinomial distribution and endpoints.}

\item{seed}{The seed used to generate the responses/toxicities with the multinomial distribution and the times.}
}
\value{
A list of \code{n_sim} elements. Each element is a data.frame of 14 columns :
\itemize{
\item ttt: the treatment arm;
\item V1, V2, V3 and V4: the Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox endpoints;
\item nb_pat: the index of the patient in his arm;
\item analyse: the index of the interim analysis;
\item eff, tox: 1 if the patient experience response/toxicity, 0 otherwise;
\item temps_recrutement: time of the patient's start in the trial;
\item temps_eff, temps_tox: time of onset of response/toxicity;
\item temps_obseff, temps_obstox: end of observation period for the patient.
}
}
\description{
Perform \code{n_sim} trials with interim analyses after recruitment of \code{ana_inter} patients.
There can be 1 or more treatment arms with or without a control group. Recruitment of patients is
deemed uniform with a mean interpatient of \code{interpatient}. Efficacy and toxicity are generated
with multinomial distribution and their time of event with the specified distribution. But for now,
TOP is only single-arm.
}
\examples{
# Generate 20 trials of 40 patients with a mean accrual of 2 patients per month
# and an observation window of 120 days for efficacy and 90 days for toxicity with
# an uniform distribution
gen_patients_multinomTOP(n_sim = 20,
                         loi_gen_patients = "unif",
                         ana_inter = rep(10, 4),
                         max_teff = 120, max_ttox = 90, interpatient = 15,
                         multinom_ttt = c(.15, .3, .15, .4))
}
