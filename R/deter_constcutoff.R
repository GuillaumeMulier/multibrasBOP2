#' Return the max number of responses/non toxicities before being greater than cutoff.
#'
#' @param cutoff Thresholds for efficacy and toxicity.
#' @param prior Lax of prior (by default, under H0).
#' @param ntr Number of enrolled patients in 1 arm.
#' @param phi Theoretical thresholds to which we compare the efficacy and toxicity (by default, values under H0).
#' @param mat_beta_xi Link matrix between endpoints and outcomes.
max_resp_const <- function(cutoff,
                           prior,
                           ntr,
                           ntr_cont = NULL,
                           phi,
                           delta = c(0, 0),
                           mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE)) {

  # Prior vector
  vec_priors <- c(sum(prior[mat_beta_xi[1, ] == 1]), sum(prior[mat_beta_xi[2, ] == 1]))

  if (is.null(ntr_cont)) { # Comparison with a reference value
    ncut <- numeric(length = length(phi))

    for (i in seq_along(phi)) {
      ncut[i] <- sum(1 - pbeta(phi[i],
                               0:ntr + vec_priors[i],
                               ntr:0 + sum(prior) - vec_priors[i]) < cutoff[i]) - 1
    }

    return(ncut)

  } else { # Comparison with a control group
    tableau <- expand.grid(y_cont = (seq_len(ntr_cont + 1) - 1),
                           y_ttt = (seq_len(ntr + 1) - 1))
    tableau$p_eff <- NA
    tableau$p_notox <- NA

    for (i in (seq_len(ntr_cont + 1) - 1)) { # Loop over responses/non toxicities in control group

      for (j in (seq_len(ntr + 1) - 1)) { # Loop over responses/non toxicities in treatment group

        # P(pi_k > pi_0 + delta | Dn) from Hobbs, 2018
        tableau[tableau$y_cont == i & tableau$y_ttt == j, "p_eff"] <- integrate(
                    f = function(x) {
                      (1 - pbeta(x + delta[1], vec_priors[1] + j, 1 - vec_priors[1] + ntr - j)) *
                        dbeta(x, vec_priors[1] + i, 1 - vec_priors[1] + ntr_cont - i)
                    },
                    lower = 0, upper = 1 - delta[1])$value

        tableau[tableau$y_cont == i & tableau$y_ttt == j, "p_notox"] <- integrate(
                    f = function(x) {
                      (1 - pbeta(x + delta[2], vec_priors[2] + j, 1 - vec_priors[2] + ntr - j)) *
                        dbeta(x, vec_priors[2] + i, 1 - vec_priors[2] + ntr_cont - i)
                    },
                    lower = 0, upper = 1 - delta[2])$value
      }

    }
    tableau$deci_eff <- tableau$p_eff < cutoff[1] # Compared to Cn, 2 thresholds this time
    tableau$deci_notox <- tableau$p_notox < cutoff[2]
    tableau <- dplyr::arrange(tableau, y_cont, y_ttt)

    return(tableau)

  }

}


#' Stopping boundaries for constant threshold
#'
#' Get the table of stopping boundaries for a given couple of thresholds.
#'
#' Generation of a table for stopping boundaries. It can be:
#' \itemize{
#'   \item for a design vs a reference value: a data.frame with the number of patients at each interim analysis with efficacy's thresholds (stop for futility if
#'   less or equal responses than the threshold), non toxicity's thresholds (stop for toxicity if less or equal non toxicities than the threshold) and toxicity's thresholds
#'   (stop for toxicity if more or equal toxicities than this threshold). If all groups aren't of the same size for all treatments, data in long format with a column
#'   ttt indicating what treatment group the tresholds are for.
#'   \item for a design vs a control group: a data.frame with the efficacy's thresholds (stop for futility if less or equal responses than the threshold when there are
#'   y_cont responses in the control group) and non toxicity's threshold (stop for toxicity if less or equal non toxicities than the threshold when there are
#'   y_cont non toxicities in the control group). If all groups aren't of the same size for all treatments, data in long format with a column
#'   ttt indicating what treatment group the tresholds are for.
#' }
#'
#' @encoding UTF-8
#'
#' @param ana_inter Vector giving the number of additional patients at each interim analysis.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value). 1 by default for all groups.
#' @param seuils The 2 constant probability thresholds for efficacy and toxicity.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. If no control group, leave NULL. Should be between -1 and 1.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # To get the table of stopping boundaries vs a given reference value
#' get_stopbound_const(ana_inter = rep(10, 4),
#'                     seuils = c(.67, .39),
#'                     p_n = c(.15, .3, .15, .4))
#'
#' @export
get_stopbound_const <- function(ana_inter,
                                rand_ratio = NULL,
                                seuils,
                                p_n,
                                prior = NULL,
                                mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                phi = NULL,
                                delta = NULL) {
  if (!dplyr::near(sum(p_n), 1))
    stop("Vector \"p_n\" should sum to 1.", call. = FALSE)
  if (length(p_n) != 4)
    stop("\"p_n\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (is.null(prior)) {
    prior <- p_n
  } else {
    if (length(prior) != 4)
      stop("\"prior\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
    if (!dplyr::near(sum(prior), 1))
      warning("\"prior\" should sum to 1 as it is a law of probability.")
  }
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (is.null(phi))
    phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  if (any(seuils > 1) | any(seuils < 0))
    stop("Parameters \"seuils\" should be between 0 and 1.", call. = FALSE)
  if (any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" represents the number of supplementary patients at each interim analysis and should thus be composed of integers.", call. = FALSE)
  if (is.null(rand_ratio)) {
    message("Without specified \"rand_ratio\" all groups will have the same size (same as control group).")
    rand_ratio <- 1
  } else if (!is.numeric(rand_ratio)) {
    stop("\"rand_ratio\" should be a number or a vector of numbers.", call. = FALSE)
  } else if (any(rand_ratio <= 0)) {
    stop("\"rand_ratio\" must be positive.", call. = FALSE)
  }
  if (length(unique(rand_ratio)) == 1 & length(rand_ratio) > 1) {
    rand_ratio <- rand_ratio[1]
  }

  if (length(rand_ratio) == 1) { # All groups have the same number of patients

    # Reassessment of the number of patients
    ana_inter_ttt <- ana_inter * rand_ratio
    if (any(!dplyr::near(ana_inter_ttt %% 1, 0))) {
      warning("\"rand_ratio\" x \"ana_inter\" isn't always an integer.
              Converted to integers with function ceiling.",
              call. = FALSE, immediate. = TRUE)
      ana_inter_ttt <- ceiling(ana_inter_ttt)
    }

    N <- sum(ana_inter_ttt)
    ana_inter_cum <- cumsum(ana_inter_ttt)

    # No control group
    if (is.null(delta)) {
      stopbound <- lapply(
        X = ana_inter_cum,
        FUN = function(nb) {
          max_resp_const(
            cutoff = seuils,
            prior = prior,
            mat_beta_xi = mat_beta_xi,
            ntr = nb,
            phi = phi
          )
        }
      )
      stopbound <- do.call("cbind", stopbound)
      stopbound <- stopbound %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(
          nb_ana = seq_along(ana_inter),
          tot_pat = ana_inter_cum,
          .before = 1
        ) %>%
        dplyr::mutate(seuil_tox = tot_pat - V2)
      names(stopbound)[3:4] <- c("seuil_eff", "seuil_notox")

    } else { # There is a control group

      # Sample size with control group
      ana_inter_cum_cont <- cumsum(ana_inter)
      N0 <- sum(ana_inter)

      stopbound <- purrr::pmap(
        .l = list(n0 = ana_inter_cum_cont, nt = ana_inter_cum),
        .f = function(n0, nt) {
          max_resp_const(
            cutoff = seuils,
            prior = prior,
            mat_beta_xi = mat_beta_xi,
            ntr = nt,
            ntr_cont = n0,
            phi = phi,
            delta = delta
          )
        }
      )

      stopbound <- lapply(stopbound, function(x) {
        x %>%
          dplyr::group_by(y_cont) %>%
          dplyr::summarise(seuil_eff = sum(deci_eff) - 1,
                           seuil_notox = sum(deci_notox) - 1)
      })
      stopbound <- lapply(seq_along(stopbound), function(x) cbind(stopbound[[x]], nb_ana = x))
      stopbound <- do.call("rbind", stopbound)

    }

  } else { # All groups don't have the same number of patients

    ana_inter_ttt <- lapply(rand_ratio, function(x) ana_inter * x)
    if (lapply(ana_inter_ttt, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
      warning("\"rand_ratio\" x \"ana_inter\" aren't always integers.
              Converted to integers with function ceiling.", call. = FALSE, immediate. = TRUE)
      ana_inter_ttt <- lapply(ana_inter_ttt,  ceiling)
    }
    if (is.null(names(ana_inter_ttt))) {
      names(ana_inter_ttt) <- paste0("ttt", seq_along(ana_inter_ttt))
    }

    N <- lapply(ana_inter_ttt, sum)
    ana_inter_cum <- lapply(ana_inter_ttt, cumsum)

    if (is.null(delta)) { # Vs reference value

      stopbound <- purrr::pmap(
        .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt),
        .f = function(.x, .y, .z) {
          threshold <- lapply(
            X = .x,
            FUN = function(nb) {
              max_resp_const(
                cutoff = seuils,
                prior = prior,
                mat_beta_xi = mat_beta_xi,
                ntr = nb,
                phi = phi
              )
            }
          )
          threshold <- do.call("cbind", threshold)
          threshold <- threshold %>%
            t() %>%
            as.data.frame() %>%
            dplyr::mutate(nb_ana = seq_along(.z),
                          tot_pat = .x,
                          .before = 1) %>%
            dplyr::mutate(seuil_tox = tot_pat - V2)
          names(threshold)[3:4] <-
            c("seuil_eff", "seuil_notox")
          return(threshold)
        }
      )

      stopbound <- map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))

    } else { # Vs control group

      ana_inter_cum_cont <- cumsum(ana_inter)
      N0 <- sum(ana_inter)

      stopbound <- purrr::pmap(
        .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt),
        .f = function(.x, .y, .z) {
          threshold <-
            purrr::pmap(
              .l = list(n0 = ana_inter_cum_cont, nt = .x),
              .f = function(n0, nt) {
                max_resp_const(
                  cutoff = seuils,
                  prior = prior,
                  mat_beta_xi = mat_beta_xi,
                  ntr = nt,
                  ntr_cont = n0,
                  phi = phi,
                  delta = delta
                )
              }
            )
          threshold <- lapply(threshold,
                              function(x) {
                                x %>%
                                  dplyr::group_by(y_cont) %>%
                                  dplyr::summarise(
                                    seuil_eff = sum(deci_eff) - 1,
                                    seuil_notox = sum(deci_notox) - 1
                                  )
                              })
          threshold <- lapply(seq_along(threshold),
                              function(x) {
                                cbind(threshold[[x]], nb_ana = x)
                              })
          threshold <- do.call("rbind", threshold)
        }
      )

      stopbound <- map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))

    }

  }

  return(stopbound)
}


#' Get the operating characteristics for a given couple gamma_eff/gamma_tox et a given list of trials/
#'
#' @param n_sim Number of simulated trials.
#' @param ana_inter Vector giving the number of additional patients at each interim analysis.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value). 1 by default for all groups.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param seuils The 2 constant probability thresholds for efficacy and toxicity.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. If no control group, leave NULL. Should be between -1 and 1.
#' @param liste_patients Data.frame of generated trials with \code{gen_patients_multinom}.
#'
#' @importFrom magrittr %>%
getoc_tox_const <- function(ana_inter,
                            rand_ratio = NULL,
                            prior = NULL,
                            p_n,
                            seuils,
                            mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                            phi = NULL,
                            delta = c(0, 0),
                            liste_patients) {

  if (is.null(prior)) prior <- p_n

  if (is.null(phi)) phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))

  if ("ttt0" %nin% unique(liste_patients$ttt)) { # Vs a reference value

    # Stopbound table for each analysis and each treatment
    stopbound <- suppressWarnings(
      get_stopbound_const(
        ana_inter = ana_inter,
        rand_ratio = rand_ratio,
        seuils = seuils,
        p_n = p_n,
        prior = prior,
        mat_beta_xi = mat_beta_xi,
        phi = phi
      )
    )

    # Add stopping criteria to the patients' counts
    if ("ttt" %in% names(stopbound)) {
      liste_patients <- dplyr::left_join(liste_patients, stopbound, by = c("nb_ana", "ttt", "tot_pat"))
    } else {
      liste_patients <- dplyr::left_join(liste_patients, stopbound, by = c("nb_ana", "tot_pat"))
    }

    # Decision for each trial
    liste_patients <- liste_patients %>%
      dplyr::mutate(decision = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Continue",
                                                (nb_ana == max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Accept the treatment",
                                                (nb_ana != max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Early stopping",
                                                (nb_ana == max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Stopping")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

  } else { # Vs control group

    # Stopbound table for each variable and interim analysis
    stopbound <- suppressWarnings(
      get_stopbound_const(
        ana_inter = ana_inter,
        rand_ratio = rand_ratio,
        seuils = seuils,
        p_n = p_n,
        prior = prior,
        mat_beta_xi = mat_beta_xi,
        delta = delta,
        phi = phi
      )
    )

    # Convert to wide format to oppose control group and each treatment group
    liste_patients <- liste_patients %>%
      tidyr::pivot_wider(names_from = ttt, values_from = efftox:tot_notox) %>%
      tidyr::pivot_longer(cols = tidyselect::matches("_ttt[^0]$"), names_pattern = "^(.*)_(ttt\\d+)$", names_to = c("colonne", "ttt")) %>%
      tidyr::pivot_wider(names_from = colonne, values_from = value)

    # Adding of stopping criteria
    if ("ttt" %in% names(stopbound)) {
      liste_patients <- liste_patients %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_notox), by = c("nb_ana", "ttt", "tot_eff_ttt0" = "y_cont")) %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_eff), by = c("nb_ana", "ttt", "tot_notox_ttt0" = "y_cont"))
    } else {
      liste_patients <- liste_patients %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_notox), by = c("nb_ana", "tot_eff_ttt0" = "y_cont")) %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_eff), by = c("nb_ana", "tot_notox_ttt0" = "y_cont"))
    }

    # Decision for each trial
    liste_patients <- liste_patients %>%
      dplyr::mutate(decision = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Continue",
                                    (nb_ana == max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Accept the treatment",
                                    (nb_ana != max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Early stopping",
                                    (nb_ana == max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Stopping")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  # Output of the function
  if (length(unique(liste_patients$ttt[liste_patients$ttt != "ttt0"])) == 1) {
    return(
      list(
        nb_accept = mean(liste_patients$decision == "Accept the treatment"),
        nb_patients_moy = mean(liste_patients$tot_pat),
        arret_premature_moy = mean(liste_patients$decision == "Early stopping"),
        nb_rejet_moy = mean(stringr::str_detect(liste_patients$decision, "Stopping"))
      )
    )
  } else {
    rejet_glob <- liste_patients %>%
      dplyr::mutate(decision_h0 = decision == "Accept the treatment") %>%
      dplyr::group_by(n_simu) %>%
      dplyr::summarise(decision_h0 = sum(decision_h0) > 0) %>%
      dplyr::pull(decision_h0) %>%
      mean()

    rejet_ttt1 <- liste_patients %>%
      dplyr::filter(ttt == "ttt1") %>%
      dplyr::mutate(decision = decision == "Accept the treatment") %>%
      dplyr::pull(decision) %>%
      mean()

    return(list(rejet_glob = rejet_glob,
                rejet_ttt1 = rejet_ttt1))

  }

}


#' Thresholds optimization
#'
#' Optimize gamma_eff and gamma_tox (constant thoughout the trial). The stopping rule can be written as stop if Pr(mat_beta_xi * p_ttt <= phi | Dn) > 1 - gamma_eff or 1 - gamma_tox.
#'
#' Optimize gamma_eff and gamma_tox:
#' \itemize{
#'   \item for 1 arm: screening of every couple gamma_eff/gamma_tox by simulating monoarm trials under H0 to assess the control of type I error rate and under
#'   H1 to maximize the power.
#'   \item for multi-arm trials: screening of every couple gamma_eff/gamma_tox by simulating trials with all arms under H0 to assess the control of FWER and under
#'   LFC to maximize the power (defined as the acceptation of the treatment at least in the effective and not toxic arm).
#' }
#'
#' @return
#' Display of the selected couple in console.
#' When stocked in an object, a list containing:
#' \itemize{
#'   \item a 1st element with the chosen couple and the type I error rate and power;
#'   \item a 2nd element with the same informations if we take the same value for efficacy and toxicity;
#'   \item a 3rd element with the matrix used to chose the optimal couple.
#' }
#'
#' @encoding UTF-8
#'
#' @param alpha Desired type I error rate.
#' @param n_bras Number of treatment arms.
#' @param bonf TRUE to apply a Bonferroni correction.
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1. Set to NULL (default value) for an uncontrolled setting.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param p_a Law of probability under the efficacy/non toxicity hypothesis.
#' @param nsim_oc Number of simulated trials.
#' @param seed Seed for \code{rmultinom}.
#' @param methode Method used to generate the data:
#' \itemize{
#'   \item 1: trial by trial;
#'   \item 2: trial by trial, patient by patient;
#'   \item 3: whole in 1.
#' }
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value). Not currently implemented.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param seq_eff,seq_tox Vectors of values used to choose the couple while optimizing the couple gamma_eff/gamma_tox.
#' @param affich_mat Yes to display in console the matrix used for the choice, No for no display and default option NULL will ask you.
#'
#' @examples
#' # In the following, not running the code because it can take some time
#' if (FALSE) {
#' # For the monoarm version of BOP2 with constant thresholds
#' carac_opt <- deter_constcutoff(alpha = 0.1,
#'                                ana_inter = c(10, rep(5, 6)),
#'                                p_n = c(0.15, 0.3, 0.15, 0.4),
#'                                p_a = c(0.18, 0.42, 0.02, 0.38),
#'                                methode = 3L,
#'                                affich_mat = "No")
#' carac_opt[[1]]
#'
#' # For a trial with 3 arms
#' carac_opt <- deter_constcutoff(alpha = 0.1,
#'                                ana_inter = rep(15, 4),
#'                                p_n = c(0.1, 0.2, 0.2, 0.5),
#'                                p_a = c(0.15, 0.5, 0.05, 0.3),
#'                                n_bras = 3,
#'                                affich_mat = "No",
#'                                bonf = FALSE)
#' carac_opt[[1]]
#' }
#'
#' @export
deter_constcutoff <- function(alpha,
                              ana_inter,
                              n_bras = 1,
                              rand_ratio = NULL,
                              bonf = FALSE,
                              mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                              phi = NULL,
                              delta = NULL,
                              prior = NULL,
                              p_n,
                              p_a,
                              nsim_oc = 10000,
                              seq_eff = seq(.01, .99, by = .02),
                              seq_tox = seq(.01, .99, by = .02),
                              seed = 1024,
                              methode = 2L,
                              affich_mat = NULL,
                              Parallele = FALSE) {

  # Check arguments
  if (is.null(prior)) {
    prior <- p_n
    message("\"p_n\" taken as prior.")
  } else {
    if (length(prior) != 4)
      stop("\"prior\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
    if (!dplyr::near(sum(prior), 1))
      warning("Vector \"prior\" should sum to 1.")
  }
  if (alpha > 1 | alpha < 0)
    stop("Alpha-risk should be a real between 0 and 1.", call. = FALSE)
  if (alpha < 0.01 | alpha > 0.3)
    warning(paste0("The specified alpha-risk is a little extreme: ", alpha, "."), call. = FALSE, immediate. = TRUE)
  if (n_bras <= 0 | n_bras %% 1 != 0)
    stop("Number of arms should be an integer >= 1.", call. = FALSE)
  if (!is.null(affich_mat))
    affich_mat <- match.arg(affich_mat, c("Yes", "yes", "No", "no"))
  if (methode %nin% c(1L, 2L, 3L))
    stop(
      "Choice between 3 methods with an integer:
              * 1L = trial by trial;
              * 2L = trial by trial, patient by patient;
              * 3L = whole in 1.
      Speed is following: 3L > 1L > 2L.",
      call. = FALSE
    )
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (length(p_n) != 4)
    stop("\"p_n\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (!dplyr::near(sum(p_n), 1))
    stop("Vector \"p_n\" should sum to 1.", call. = FALSE)
  if (length(p_a) != 4)
    stop("\"p_a\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (!dplyr::near(sum(p_a), 1))
    stop("Vector \"p_a\" should sum to 1.", call. = FALSE)
  if (length(nsim_oc) != 1 || !is.numeric(nsim_oc) || nsim_oc < 0 || nsim_oc %% 1 != 0)
    stop("\"nsim_oc\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (any(delta > 1) | any(delta < -1))
    stop("\"delta\" should be a real between -1 and 1.")

  if (is.null(rand_ratio)) {
    rand_ratio <- 1
  } else if (!is.numeric(rand_ratio)) {
    stop("\"rand_ratio\" should be (a) number(s).", call. = FALSE)
  } else if (any(rand_ratio <= 0)) {
    stop("\"rand_ratio\" should be positive.", call. = FALSE)
  } else if (length(rand_ratio) == 1 & n_bras != 1) {
    message(paste0("All treatment groups will have ", rand_ratio, " x \"ana_inter\"/\"ana_inter_tox\" patients at each interim analysis."))
  } else if (length(rand_ratio) != n_bras) {
    stop("\"rand_ratio\" should have ", n_bras, " specified values or only 1.", call. = FALSE)
  }

  prompt_affich <- function(affich_mat = NULL) {
    while (is.null(affich_mat)) {
      affich_mat <-
        readline(prompt = "Do you want to display the C_/gamma matrix? (Yes/No) > ")
      if (!stringr::str_detect(affich_mat, "Yes|yes|No|no")) {
        affich_mat <- NULL
        message("Answer with Yes or No.")
      }

    }

    return(affich_mat)

  }

  if (is.null(affich_mat))
    affich_mat <- prompt_affich()

  # Bonferroni correction for type I error rate
  if (bonf & n_bras > 1) {
    alpha <- alpha / n_bras
    n_bras <- 1
  }

  if (n_bras > 1) {
    # Reconstruction of LFC
    multinom_ttt_lfc <-
      matrix(c(p_a, rep(p_n, n_bras - 1)), nrow = 4) %>%
      as.data.frame() %>%
      as.list()
    names(multinom_ttt_lfc) <- paste0("ttt", seq_len(n_bras))

    # Reconstruction of the global null hypothesis
    multinom_ttt_H0 <- matrix(rep(p_n, n_bras), nrow = 4) %>%
      as.data.frame() %>%
      as.list()
    names(multinom_ttt_H0) <- paste0("ttt", seq_len(n_bras))

  }

  # Matrix of possible values for gamma_eff and gamma_tox
  matrice_carac <- expand.grid(seq_eff, seq_tox) %>%
    as.matrix()

  # Monoarm BOP2
  if (n_bras == 1) {
    # Generation of dataset with the chosen method
    suppressMessages(
      tab_h0 <-
        gen_patients_multinom(
          n_sim = nsim_oc,
          ana_inter = ana_inter,
          rand_ratio = rand_ratio,
          multinom_ttt = list(p_n),
          multinom_cont = if (is.null(delta)) {NULL} else {p_n},
          seed = seed,
          methode = methode
        )
    )
    suppressMessages(suppressWarnings(
      tab_h1 <-
        gen_patients_multinom(
          n_sim = nsim_oc,
          ana_inter = ana_inter,
          rand_ratio = rand_ratio,
          multinom_ttt = list(p_a),
          multinom_cont = if (is.null(delta)) {NULL} else {p_n},
          seed = seed,
          methode = methode
        )
    ))

  } else {
    # Generation of dataset with the chosen method
    suppressMessages(
      tab_h0 <-
        gen_patients_multinom(
          n_sim = nsim_oc,
          ana_inter = ana_inter,
          rand_ratio = rand_ratio,
          multinom_ttt = multinom_ttt_H0,
          multinom_cont = if (is.null(delta)) {NULL} else {p_n},
          seed = seed,
          methode = methode
        )
    )
    suppressMessages(suppressWarnings(
      tab_h1 <-
        gen_patients_multinom(
          n_sim = nsim_oc,
          ana_inter = ana_inter,
          rand_ratio = rand_ratio,
          multinom_ttt = multinom_ttt_lfc,
          multinom_cont = if (is.null(delta)) {NULL} else {p_n},
          seed = seed,
          methode = methode
        )
    ))

  }

  # Get the type I error rate and power for each couple
  # if (future::availableCores() == 1 | !Parallele) {
  #   plan_or <- future::plan(future::sequential)
  # } else {
  #   plan_or <- future::plan(future::multisession, workers = future::availableCores() - 1)
  # }
  liste_mat <- purrr::map(
    .x = seq_len(nrow(matrice_carac)),
    .f = function(x) {
      if (x %% 100 == 0) cat(x, "/", nrow(matrice_carac), "\n")
      oc_tox_n <- getoc_tox_const(
        ana_inter = ana_inter,
        rand_ratio = rand_ratio,
        seuils = c(matrice_carac[x, 1], matrice_carac[x, 2]),
        liste_patients = tab_h0,
        prior = prior,
        p_n = p_n,
        phi =  phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      oc_tox_a <- getoc_tox_const(
        ana_inter = ana_inter,
        rand_ratio = rand_ratio,
        seuils = c(matrice_carac[x, 1], matrice_carac[x, 2]),
        liste_patients = tab_h1,
        prior = prior,
        p_n = p_n,
        phi =  phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      if (n_bras == 1) {
        return(c(
          matrice_carac[x, 1],
          matrice_carac[x, 2],
          oc_tox_n$nb_accept,
          oc_tox_a$nb_accept
        ))
      } else {
        return(c(
          matrice_carac[x, 1],
          matrice_carac[x, 2],
          oc_tox_n$rejet_glob,
          oc_tox_a$rejet_ttt1
        ))
      }
    }
  )
  liste_mat <- do.call("rbind", liste_mat)
  # future::plan(plan_or)

  # Column 4 = P(Reject of H0 | H1) and column 3 = P(Reject H0 | H0)
  # Determine the optimal couple
  carac_optimale <- liste_mat[liste_mat[, 3] <= alpha,]
  carac_optimale <- carac_optimale[carac_optimale[, 4] == max(carac_optimale[, 4]),]
  if (is.matrix(carac_optimale)) carac_optimale <- carac_optimale[1,]

  carac_optimale <- as.double(carac_optimale)
  gamma_eff <- carac_optimale[1]
  gamma_tox <- carac_optimale[2]
  alpha_calc <- carac_optimale[3]
  puissance_calc <- carac_optimale[4]

  # Couple if gamma_eff = gamma_tox
  carac_optimale_sensvie <- liste_mat[liste_mat[, 1] == liste_mat[, 2],]
  carac_optimale_sensvie <- carac_optimale_sensvie[carac_optimale_sensvie[, 3] <= alpha,]
  carac_optimale_sensvie <- carac_optimale_sensvie[carac_optimale_sensvie[, 4] == max(carac_optimale_sensvie[, 4]),]
  if (is.matrix(carac_optimale_sensvie)) carac_optimale_sensvie <- carac_optimale_sensvie[1,]
  carac_optimale_sensvie <- as.double(carac_optimale_sensvie)
  gamma_eff_sensvie <- carac_optimale_sensvie[1]
  gamma_tox_sensvie <- carac_optimale_sensvie[2]
  alpha_calc_sensvie <- carac_optimale_sensvie[3]
  puissance_calc_sensvie <- carac_optimale_sensvie[4]

  cat("The optimal couple is: gamma_eff =", gamma_eff, "and gamma_tox =", gamma_tox, "\n")
  cat("Estimated type I error rate: ", alpha_calc, " and estimated power: ", puissance_calc, "\n")
  if (affich_mat == "Oui") print(matrice_carac)

  invisible(list(
    seuils = c(
      gamma_eff = gamma_eff,
      gamma_tox = gamma_tox,
      alpha_calc = alpha_calc,
      puissance_calc = puissance_calc
    ),
    seuils_sensvie = c(
      gamma_eff_sensvie = gamma_eff_sensvie,
      gamma_tox_sensvie = gamma_tox_sensvie,
      alpha_calc_sensvie = alpha_calc_sensvie,
      puissance_calc_sensvie = puissance_calc_sensvie
    ),
    mat = liste_mat
  ))

}
