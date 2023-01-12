#' @title
#' Generate a time of entry or of event
#'
#' @description
#' The function generates times for patients with the specified law. There are actually
#' 3 laws supported : uniform, Weibull and log-logistic.
#' It is assumed that the time is between 0 and \code{tmax}. With uniform distribution, 50% of the patients will have
#' a time between 0 and \code{tmax}/2; and with Weibull/log-logistic distributions, \code{prop_demi}% of
#' the patients will have a time between 0 and \code{tmax}/2.
#'
#' @return A vector of \code{n} times between 0 and \code{tmax} simulated according to the wanted distribution.
#'
#' @param loi The distribution law (uniform, weibull or log-logistic).
#' @param n An positive integer giving the number of patients.
#' @param tmax A positive real giving the maximum time of the observation window.
#' @param ptmax For Weibull and log-logistic distributions, the proportion of patients with endpoint
#' at \code{tmax}. (But all n times will be between 0 and \code{tmax}. It is used to determine the 2
#' parameters of the distribution).
#' @param prop_demi For Weibull and log-logistic distributions, the proportion of \code{ptmax} at \code{tmax}/2.
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' # Generate the time of arrival of 10 successive patients with a uniform
#' # recruitment of mean 30 days
#' time_pat <- gen_temps_patient(loi = "unif",
#'                               n = 10L,
#'                               tmax = 2 * 30)
#' cumsum(time_pat)
#'
#' # Generate the time of efficacy of 20 patients after their recruitment
#' # with a Weibull distribution between 0 and 100 days with 75% of the responses
#' # in the first half of the window.
#' gen_temps_patient(loi = "weib",
#'                   n = 20L,
#'                   tmax = 100, ptmax = .8,
#'                   prop_demi = .75)
gen_temps_patient <- function(loi = c("uniform", "weibull", "log-logistic"),
                              n = 1L,
                              tmax, ptmax = 1,
                              prop_demi = .5) {

  loi <- match.arg(loi, c("uniform", "weibull", "log-logistic"))
  if (!is.numeric(tmax) || tmax <= 0)
    stop("\"tmax\" should be a positive real number.", call. = FALSE)
  if (!is.numeric(n) || n <= 0 || n %% 1 != 0)
    stop("\"n\" should be a strictly positive integer.", call. = FALSE)

  if (loi == "uniform") {
    # Uniform law so 50% of patients at Tmax/2

    obs_patient <- runif(n, min = 0, max = tmax)

  } else {

    if (!is.numeric(prop_demi) || prop_demi < 0 || prop_demi > 1)
      stop("\"prop_demi\" should be a positive real number between 0 and 1.", call. = FALSE)
    if (!is.numeric(ptmax) || ptmax < 0 || ptmax > 1)
      stop("\"ptmax\" should be a positive real number between 0 and 1.", call. = FALSE)

    if (loi == "weibull") {
      # Weibull distibution:
      # Cumulative distribution function: F(t) = 1 - exp(-(lambda * t) ^ alpha)
      # We have F(Tmax) = ptmax and F(Tmax/2) = ptdemi = prop_demi * ptmax
      # ==> prop_demi % in first half of the following

      ptdemi      <- prop_demi * ptmax
      alpha       <- log(log(1 - ptmax) / log(1 - ptdemi)) / log(2)
      lambda      <- (-log(1 - ptmax)) ^ (1 / alpha) / tmax
      stmax       <- exp(-(lambda * tmax) ^ alpha)
      obs_patient <- (-log(runif(n, stmax, 1))) ^ (1 / alpha) / lambda

    } else if (loi == "log-logistic") {
      # Log-logistic distribution:
      # Cumulative ditribution function: F(t) = 1 / (1 + (lambda * t) ^ (-alpha))
      # We have F(Tmax) = ptmax and F(Tmax/2) = ptdemi = prop_demi * ptmax
      # ==> prop_demi % in first half of the following

      ptdemi      <- prop_demi * ptmax
      alpha       <- log((1 / ptdemi - 1) / (1 / ptmax - 1)) / log(2)
      lambda      <- (ptmax / (1 - ptmax)) ^ (1 / alpha) / tmax
      stmax       <- 1 / (1 + (lambda * tmax) ^ alpha)
      obs_patient <- ((1 / runif(n, stmax, 1) - 1)) ^ (1 / alpha) / lambda

    }

  }

  return(obs_patient)

}


#' @title
#' Generate a list of trials with time to event for each patient
#'
#' @description
#' Perform \code{n_sim} trials with interim analyses after recruitment of \code{ana_inter} patients.
#' There can be 1 or more treatment arms with or without a control group. Recruitment of patients is
#' deemed uniform with a mean interpatient of \code{interpatient}. Efficacy and toxicity are generated
#' with multinomial distribution and their time of event with the specified distribution. But for now,
#' TOP is only single-arm.
#'
#' @return
#' A list of \code{n_sim} elements. Each element is a data.frame of 14 columns :
#' \itemize{
#'   \item ttt: the treatment arm;
#'   \item V1, V2, V3 and V4: the Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox endpoints;
#'   \item nb_pat: the index of the patient in his arm;
#'   \item analyse: the index of the interim analysis;
#'   \item eff, tox: 1 if the patient experience response/toxicity, 0 otherwise;
#'   \item temps_recrutement: time of the patient's start in the trial;
#'   \item temps_eff, temps_tox: time of onset of response/toxicity;
#'   \item temps_obseff, temps_obstox: end of observation period for the patient.
#' }
#'
#'
#' @param n_sim The number of simulated trials.
#' @param ana_inter A vector of the supplementray number of patients at each interim analysis.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio A vector of the same length of \code{multinom_ttt} or length 1 giving the ratio of patients between control and treatment groups.
#' @param multinom_cont Vector of length 4 giving the law of control group (Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox). NULL if no control arm.
#' @param multinom_ttt List of vector of length 4 giving the law of treatment groups (Eff/Tox, Eff/NoTox, NoEff/Tox and NoEff/NoTox). 1 vector if 1 group is accepted.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints.
#' @param interpatient Mean time between 2 patients.
#' @param loi_gen_patients The law of generation of time to event (uniform, weibull, log-logistic).
#' @param max_teff The length of the observation window for efficacy.
#' @param prop_demieff The proportion of responses at \code{max_teff}/2.
#' @param max_ttox The length of the observation window for toxicity.
#' @param prop_demitox The proportion of toxicities at \code{max_ttox}/2.
#' @param seed The seed used to generate the responses/toxicities with the multinomial distribution and the times.
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # Generate 20 trials of 40 patients with a mean accrual of 2 patients per month
#' # and an observation window of 120 days for efficacy and 90 days for toxicity with
#' # an uniform distribution
#' gen_patients_multinomTOP(n_sim = 20,
#'                          loi_gen_patients = "unif",
#'                          ana_inter = rep(10, 4),
#'                          max_teff = 120, max_ttox = 90, interpatient = 15,
#'                          multinom_ttt = c(.15, .3, .15, .4))
gen_patients_multinomTOP <- function(n_sim,
                                     ana_inter, ana_inter_tox = NULL,
                                     interpatient, loi_gen_patients = "unif",
                                     max_teff, prop_demieff = .5, max_ttox, prop_demitox = .5,
                                     rand_ratio = NULL,
                                     multinom_cont = NULL,
                                     multinom_ttt = list(),
                                     mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                     seed = 1024) {

  # Checks arguments
  if (length(n_sim) != 1 || !is.numeric(n_sim) || n_sim <= 0 || n_sim %% 1 != 0)
    stop("\"n_sim\" shoud be a strictly positive integer.", call. = FALSE)
  if (any(!is.numeric(ana_inter)) || any(ana_inter < 0) || any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.matrix(mat_beta_xi))
    stop("\"mat_beta_xi\" must be a matrix.", call. = FALSE)
  if (any(dim(mat_beta_xi) != c(2, 4)))
    stop("The matrix \"mat_beta_xi\" should be a 2x4 matrix.", call. = FALSE)
  if (length(interpatient) != 1 || !is.numeric(interpatient) || interpatient <= 0)
    stop("\"interpatient\" should be an positive real.", call. = FALSE)
  if (length(max_teff) != 1 || !is.numeric(max_teff) || max_teff <= 0)
    stop("\"max_teff\" should be an positive real.", call. = FALSE)
  if (length(max_ttox) != 1 || !is.numeric(max_ttox) || max_ttox <= 0)
    stop("\"max_ttox\" should be an positive real.", call. = FALSE)

  if (length(multinom_ttt) == 0)
    stop("No treatment probability specified., call. = FALSE")
  if (!is.list(multinom_ttt)) {
    multinom_ttt <- list(ttt1 = multinom_ttt)
    warning("\"multinom_ttt\" should be a list. Converted to a list with the command list(multinom_ttt) assuming only 1 treatment.", call. = FALSE, immediate. = TRUE)
    message(multinom_ttt)
  }
  if (any(lapply(multinom_ttt, length) != 4))
    stop("Each law of probability in \"multinom_ttt\" must be of length 4.", call. = FALSE)
  if (any(unlist(lapply(
    multinom_ttt,
    FUN = function(x) !dplyr::near(sum(x), 1)
  ))))
    stop("Each law of probability in \"multinom_ttt\" must sum to 1.", call. = FALSE)
  if (is.null(multinom_cont)) {
    message("No control group.")
    if (is.null(names(multinom_ttt))) {
      noms_ttt <- c(paste0("ttt", seq_len(length(multinom_ttt))))
    } else {
      noms_ttt <- c(names(multinom_ttt))
    }
    proba <- multinom_ttt
    names(proba) <- noms_ttt
  } else {
    if (length(multinom_cont) != 4)
      stop("The probability distribution of controls should be a vector of length 4.", call. = FALSE)
    if (!dplyr::near(sum(multinom_cont), 1))
      stop("The probability distribution of controls should sum to 1.", call. = FALSE)
    if (is.null(names(multinom_ttt))) {
      noms_ttt <- c("ttt0", paste0("ttt", seq_len(length(multinom_ttt))))
    } else {
      noms_ttt <- c("ttt0", names(multinom_ttt))
    }
    proba <- append(list(multinom_cont), multinom_ttt)
    names(proba) <- noms_ttt
  }
  # ttt0 is always for controls, and ttt1, 2, ... for evaluated treatments.
  # If no control treatment, start at ttt1.

  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.")
    anas_inters_cum <- sort(union(cumsum(ana_inter), cumsum(ana_inter_tox)))
    ana_inter     <- c(anas_inters_cum[1], diff(anas_inters_cum))
  }

  if (is.null(rand_ratio)) {
    message("Without specifying \"rand_ratio\", all groups are assumed to be of same size.")
    rand_ratio <- rep(1, length(multinom_ttt) + !is.null(multinom_cont))
  } else if (length(rand_ratio) == 1 & length(multinom_ttt) != 1) {
    message("Each arm will have \"rand_ratio\" x \"ana_inter\" patients at each interim analysis.")
    if (!is.numeric(rand_ratio) || rand_ratio < 0)
      stop("\"rand_ratio\" should be a positive real number indicating the ratio between the number of patients in treatment arms and in control arm.", call. = FALSE)
    if (is.null(multinom_cont)) {
      rand_ratio <- rep(rand_ratio, length(multinom_ttt))
    } else {
      rand_ratio <- c(1, rep(rand_ratio, length(multinom_ttt)))
    }
  } else {
    if (!is.numeric(rand_ratio) || any(rand_ratio < 0))
      stop("\"rand_ratio\" should be a vector of positive real numbers indicating the ratio between the number of patients in treatment arms and in control arm.", call. = FALSE)
    if (length(rand_ratio) != length(multinom_ttt))
      stop("\"rand_ratio\" and \"multinom_ttt\"' must have the same length.", call. = FALSE)
    if (!is.null(multinom_cont)) {
      rand_ratio <- c(1, rand_ratio)
    }
  }

  anas_inters <- lapply(rand_ratio, FUN = function(x, y) {x * y}, y = ana_inter)
  if (lapply(anas_inters, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
    warning("\"rand_ratio\" x \"ana_inter\" must only return integers.
         The numbers were converted to integers using the ceiling function.",
             call. = FALSE, immediate. = TRUE)
    anas_inters <- lapply(anas_inters,  ceiling)
  }
  names(anas_inters) <- noms_ttt
  nmax <- lapply(anas_inters, sum)
  ana_inter_cum <- lapply(anas_inters, cumsum)

  set.seed(seed)
  on.exit(set.seed(NULL), add = TRUE) # Reset the seed at end of function to avoid change of environment outside the execution of the function

  liste_essais <- purrr::map(
    .x = seq_len(n_sim),
    proba = proba,
    nmax = nmax,
    .f = function(x, proba, nmax) {
      purrr::map_dfr(
        .x = names(proba),
        .f = function(name) {
          tableau <- rmultinom(n = nmax[[name]], size = 1, prob = proba[[name]]) %>%
            t() %>%
            as.data.frame()
          tableau <- cbind(ttt = name, tableau, nb_pat = seq_len(nmax[[name]]))
          tableau$analyse <- NA
          for (i in rev(ana_inter_cum[[name]])) {tableau$analyse <- ifelse(tableau$nb_pat <= i, match(i, ana_inter_cum[[name]]), tableau$analyse)}
          tableau$eff <- tableau$V1 + tableau$V2
          tableau$tox <- tableau$V1 + tableau$V3
          tableau$temps_recrutement <- cumsum(gen_temps_patient(loi = "unif", n = nmax[[name]], tmax = interpatient * 2))
          tableau$temps_eff <- gen_temps_patient(loi = loi_gen_patients, n = nmax[[name]], tmax = max_teff, ptmax = sum(proba[[name]] * mat_beta_xi[1, ]), prop_demi = prop_demieff)
          tableau$temps_eff <- ifelse(tableau$eff == 1, tableau$temps_recrutement + tableau$temps_eff, 1e+12)
          tableau$temps_obseff <- tableau$temps_recrutement + max_teff
          tableau$temps_tox <- gen_temps_patient(loi = loi_gen_patients, n = nmax[[name]], tmax = max_ttox, ptmax = 1 - sum(proba[[name]] * mat_beta_xi[2, ]), prop_demi = prop_demitox)
          tableau$temps_tox <- ifelse(tableau$tox == 1, tableau$temps_recrutement + tableau$temps_tox, 1e+12)
          tableau$temps_obstox <- tableau$temps_recrutement + max_ttox
          return(tableau)
        }
      )
    }
  )

  return(liste_essais)

}
