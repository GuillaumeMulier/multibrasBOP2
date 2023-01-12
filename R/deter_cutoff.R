#' Get the boundaries for 1 number of patients. Give the maximum number of patients (or TRUE/FALSE in case of controlled design)
#' that satisfy Pr(p > phi | Dn) < 1 - Cn.
#'
#' @param cutoff Cn (in reality we give 1 - Cn).
#' @param prior The prior's law (by default the law under the inefficacy/toxicity hypothesis).
#' @param ntr Number of patients in treatment arm.
#' @param ntr_cont Number of patients in control arm.
#' @param delta Minimum desired difference for effect of treatment (by default 0 for efficacy and toxicity).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis).
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints.
max_resp <- function(cutoff,
                     prior,
                     ntr,
                     ntr_cont = NULL,
                     phi,
                     delta = c(0, 0),
                     mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE)) {

  # Priors vector
  vec_priors <- c(sum(prior[mat_beta_xi[1, ] == 1]), sum(prior[mat_beta_xi[2, ] == 1]))

  if (is.null(ntr_cont)) { # Comparison to a reference value
    ncut <- numeric(length = length(phi))
    for (i in seq_along(phi)) {
      ncut[i] <- sum(1 - pbeta(phi[i],
                               0:ntr + vec_priors[i],
                               ntr:0 + sum(prior) - vec_priors[i]) < cutoff) - 1
    }

    return(ncut)

  } else { # Comparison to a control group
    tableau <- expand.grid(y_cont = (seq_len(ntr_cont + 1) - 1),
                           y_ttt = (seq_len(ntr + 1) - 1))
    tableau$p_eff <- NA
    tableau$p_notox <- NA

    for (i in (seq_len(ntr_cont + 1) - 1)) { # Loop over all possible values of responses/non toxicities in control group

      for (j in (seq_len(ntr + 1) - 1)) { # Loop over all possible values of responses/non toxicities in treatment group

        tableau[tableau$y_cont == i & tableau$y_ttt == j, "p_eff"] <- integrate(
                    f = function(x) {
                      (1 - pbeta(x + delta[1], vec_priors[1] + j, 1 - vec_priors[1] + ntr - j)) *
                        dbeta(x, vec_priors[1] + i, 1 - vec_priors[1] + ntr_cont - i)
                    },
                    lower = 0, upper = 1 - delta[1]
                  )$value

        tableau[tableau$y_cont == i & tableau$y_ttt == j, "p_notox"] <- integrate(
                    f = function(x) {
                      (1 - pbeta(x + delta[2], vec_priors[2] + j, 1 - vec_priors[2] + ntr - j)) *
                        dbeta(x, vec_priors[2] + i, 1 - vec_priors[2] + ntr_cont - i)
                    },
                    lower = 0, upper = 1 - delta[2]
                  )$value

      }

    }
    tableau$deci_eff <- tableau$p_eff < cutoff
    tableau$deci_notox <- tableau$p_notox < cutoff
    tableau <- dplyr::arrange(tableau, y_cont, y_ttt)

    return(tableau)

  }

}


#' Table of the stopping boundaries
#'
#' Given a couple &lambda;/&gamma;, give the table of stopping boundaries at each interim analysis.
#'
#' Generate the stopping boundaries table different given the situation:
#' \itemize{
#'   \item against a fixed reference value: a data.frame with thresholds at each interim analysis. nb_ana give the number of patients at each
#'   interim analysis, ttt the treatment group and tot_pat the number of recruited patients. seuil_eff is the threshold for futility (stop for futility
#'   if the number of responses is less or equal to this threshold), seuil_notox and seuil_tox for toxicity (stop for toxicity if the number of non toxicities
#'   is less or equal to seuil_notox or if the number of toxicities is greater or equal to seuil_tox).
#'   \item against a control treatment: a data.frame with thresholds at each interim analysis with the values for each possible number of responses/non toxicities
#'   in the control group. seuil_eff and seuil_notox works in the same way as above.
#' }
#'
#' @encoding UTF-8
#'
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value).
#' @param C_,gamm The 2 tuning parameters of criterion Cn.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. If no control group, leave NULL. Should be between -1 and 1.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # The table from the article of Zhou et al. de 2017
#' get_stopbound(ana_inter = c(10, rep(5, 6)),
#'               C_ = 0.625, gamm = 0.98,
#'               p_n = c(0.15, 0.3, 0.15, 0.4))
#'
#' # The table for a larger trial with 3 arms
#' get_stopbound(ana_inter = rep(15, 4),
#'               C_ = 0.78, gamm = 0.9,
#'               p_n = c(0.15, 0.3, 0.15, 0.4))
#'
#' @export
get_stopbound <- function(ana_inter,
                          ana_inter_tox = NULL,
                          rand_ratio = NULL,
                          C_,
                          gamm,
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
      stop("\"prior\" should sum to 1 as it is a law of probability.", call. = FALSE)
  }
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (is.null(phi))
    phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  if (length(C_) > 1 || !is.numeric(C_) || C_ > 1 || C_ < 0)
    stop("Parameter \"C_\" should be between 0 and 1.", call. = FALSE)
  if (length(gamm) > 1 || !is.numeric(gamm) || gamm < 0)
    stop("Parameter \"gamm\" should be a positive real.", call. = FALSE)
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

  if (length(unique(rand_ratio)) == 1 & length(rand_ratio) > 1) rand_ratio <- rand_ratio[1]

  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.", call. = FALSE)
    ana_eff         <- cumsum(ana_inter)
    ana_tox         <- cumsum(ana_inter_tox)
    anas_inters_cum <- sort(union(ana_eff, ana_tox))
    anas_inters     <- c(anas_inters_cum[1], diff(anas_inters_cum))
  } else {
    ana_eff         <- cumsum(ana_inter)
    ana_tox         <- cumsum(ana_inter)
    anas_inters     <- ana_inter
  }

  if (length(rand_ratio) == 1) { # All groups have the same number of patients

    # Reassessment of the number of patients
    ana_inter_ttt <- anas_inters * rand_ratio
    if (any(!dplyr::near(ana_inter_ttt %% 1, 0))) {
      warning("\"rand_ratio\" x \"ana_inter\" isn't always an integer.
              Converted to integers with function ceiling.",
        call. = FALSE, immediate. = TRUE)
      ana_inter_ttt <- ceiling(ana_inter_ttt)
    }

    N <- sum(ana_inter_ttt)
    ana_inter_cum <- cumsum(ana_inter_ttt)

    if (is.null(delta)) { # Vs reference value
      stopbound <- lapply(
        X = ana_inter_cum,
        FUN = function(n) {
          max_resp(
            cutoff = C_ * (n / N) ^ gamm,
            prior = prior,
            mat_beta_xi = mat_beta_xi,
            ntr = n,
            phi = phi
          )
        }
      )
      stopbound <- do.call("cbind", stopbound)
      stopbound <- stopbound %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(
          nb_ana = seq_along(ana_inter_ttt),
          tot_pat = ana_inter_cum,
          .before = 1
        ) %>%
        dplyr::mutate(seuil_tox = tot_pat - V2)
      names(stopbound)[3:4] <- c("seuil_eff", "seuil_notox")
      stopbound$seuil_eff[stopbound$tot_pat %nin% ana_eff] <- -1
      stopbound$seuil_notox[stopbound$tot_pat %nin% ana_tox] <- -1
      stopbound$seuil_tox[stopbound$tot_pat %nin% ana_tox] <- stopbound$tot_pat[stopbound$tot_pat %nin% ana_tox] + 1

    } else { # Vs control group

      ana_inter_cum_cont <- cumsum(anas_inters)
      N0 <- sum(anas_inters)
      stopbound <- purrr::pmap(
        .l = list(n0 = ana_inter_cum_cont, nt = ana_inter_cum),
        .f = function(n0, nt) {
          max_resp(
            cutoff = C_ * ((n0 + nt) / (N0 + N)) ^ gamm,
            prior = prior,
            mat_beta_xi = mat_beta_xi,
            ntr = nt,
            ntr_cont = n0,
            phi = phi,
            delta = delta
          )
        }
      )
      stopbound <- lapply(stopbound,
                          function(x) {
                            x %>%
                              dplyr::group_by(y_cont) %>%
                              dplyr::summarise(seuil_eff = sum(deci_eff) - 1,
                                               seuil_notox = sum(deci_notox) - 1)
                          })
      stopbound <- lapply(seq_along(stopbound), function(x) cbind(stopbound[[x]], nb_ana = x))
      stopbound <- do.call("rbind", stopbound) %>%
        dplyr::mutate(tot_pat = ana_inter_cum_cont[match(nb_ana, seq_along(ana_inter_cum_cont))],
                      tot_pat_ttt = ana_inter_cum[match(nb_ana, seq_along(ana_inter_cum))])
      stopbound$seuil_eff[stopbound$tot_pat %nin% ana_eff] <- -1
      stopbound$seuil_notox[stopbound$tot_pat %nin% ana_tox] <- -1
    }

  } else { # Groups without the same number of patients

    ana_inter_ttt <- lapply(rand_ratio, function(x) anas_inters * x)
    ana_eff_ttt <- lapply(rand_ratio, function(x) ana_eff * x)
    names(ana_eff_ttt) <- paste0("ttt", seq_along(ana_eff_ttt))
    ana_tox_ttt <- lapply(rand_ratio, function(x) ana_tox * x)
    names(ana_tox_ttt) <- paste0("ttt", seq_along(ana_tox_ttt))
    if (lapply(ana_inter_ttt, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
      warning("\"rand_ratio\" x \"ana_inter\" isn't always an integer.
              Converted to integers with function ceiling.",
              call. = FALSE, immediate. = TRUE)
      ana_inter_ttt <- lapply(ana_inter_ttt,  ceiling)
      ana_eff_ttt <- lapply(ana_eff_ttt,  ceiling)
      ana_tox_ttt <- lapply(ana_tox_ttt,  ceiling)
    }
    if (is.null(names(ana_inter_ttt))) names(ana_inter_ttt) <- paste0("ttt", seq_along(ana_inter_ttt))

    N <- lapply(ana_inter_ttt, sum)
    ana_inter_cum <- lapply(ana_inter_ttt, cumsum)

    if (is.null(delta)) { # Vs a reference value

      stopbound <- purrr::pmap(
        .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt, ana_eff_ttt = ana_eff_ttt, ana_tox_ttt = ana_tox_ttt),
        .f = function(.x, .y, .z, ana_eff_ttt, ana_tox_ttt) {
          threshold <- lapply(
            X = .x,
            FUN = function(n) {
              max_resp(
                cutoff = C_ * (n / .y) ^ gamm,
                prior = prior,
                mat_beta_xi = mat_beta_xi,
                ntr = n,
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
          names(threshold)[3:4] <- c("seuil_eff", "seuil_notox")
          threshold$seuil_eff[threshold$tot_pat %nin% ana_eff_ttt] <- -1
          threshold$seuil_notox[threshold$tot_pat %nin% ana_tox_ttt] <- -1
          threshold$seuil_tox[threshold$tot_pat %nin% ana_tox_ttt] <- threshold$tot_pat[threshold$tot_pat %nin% ana_tox_ttt] + 1
          return(threshold)
        }
      )

      stopbound <- purrr::map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))

    } else { # Vs control group

      ana_inter_cum_cont <- cumsum(anas_inters)
      N0 <- sum(anas_inters)

      stopbound <- purrr::pmap(
          .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt, ana_eff_ttt = ana_eff_ttt, ana_tox_ttt = ana_tox_ttt),
          .f = function(.x, .y, .z, ana_eff_ttt, ana_tox_ttt) {
            threshold <- purrr::pmap(
                .l = list(n0 = ana_inter_cum_cont, nt = .x),
                .f = function(n0, nt) {
                  max_resp(
                    cutoff = C_ * ((n0 + nt) / (N0 + .y)) ^ gamm,
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
            threshold <- lapply(seq_along(threshold), function(x) cbind(threshold[[x]], nb_ana = x))
            threshold <- do.call("rbind", threshold) %>%
              dplyr::mutate(tot_pat = ana_inter_cum_cont[match(nb_ana, seq_along(ana_inter_cum_cont))],
                            tot_pat_ttt = .x[match(nb_ana, seq_along(.x))])
            threshold$seuil_eff[threshold$tot_pat_ttt %nin% ana_eff_ttt] <- -1
            threshold$seuil_notox[threshold$tot_pat_ttt %nin% ana_tox_ttt] <- -1
            return(threshold)
          }
        )

      stopbound <- purrr::map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))

    }

  }

  return(stopbound)

}


#' Operating characteristics for a couple &lambda;/&gamma; given a trials' list.
#'
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value).
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param C_,gamm The 2 hyperparameters of criterion Cn.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param liste_patients Generated trials via function \code{gen_patients_multinom}.
#'
#' @importFrom magrittr %>%
#'
#' @encoding UTF-8
getoc_tox <- function(ana_inter,
                      ana_inter_tox = NULL,
                      rand_ratio = NULL,
                      prior = NULL,
                      p_n,
                      C_,
                      gamm,
                      mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                      phi = NULL,
                      delta = c(0, 0),
                      liste_patients) {

  if (is.null(prior)) prior <- p_n

  if (is.null(phi)) phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))

  if ("ttt0" %nin% unique(liste_patients$ttt)) { # Vs a reference value

    # Generation of stopbound table at each interim analysis and for each group
    stopbound <- suppressWarnings(
      get_stopbound(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        C_ = C_,
        gamm = gamm,
        prior = prior,
        p_n = p_n,
        mat_beta_xi = mat_beta_xi,
        phi = phi
      )
    )

    # Add the stopping boundaries to the dataset
    if ("ttt" %in% names(stopbound)) {
      liste_patients <- dplyr::left_join(liste_patients, stopbound, by = c("nb_ana", "ttt", "tot_pat"))
    } else {
      liste_patients <- dplyr::left_join(liste_patients, stopbound, by = c("nb_ana", "tot_pat"))
    }

    # Decision rules
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

    # Generation of stopbound table at each interim analysis and for each group
    stopbound <- suppressWarnings(
      get_stopbound(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        C_ = C_,
        gamm = gamm,
        p_n = p_n,
        prior = prior,
        mat_beta_xi = mat_beta_xi,
        delta = delta
      )
    )

    # Wide format to compare control with treatment groups
    liste_patients <- liste_patients %>%
      tidyr::pivot_wider(names_from = ttt, values_from = efftox:tot_notox) %>%
      tidyr::pivot_longer(cols = tidyselect::matches("_ttt[^0]$"), names_pattern = "^(.*)_(ttt\\d+)$", names_to = c("colonne", "ttt")) %>%
      tidyr::pivot_wider(names_from = colonne, values_from = value)

    # Adding the stopping boundaries
    if ("ttt" %in% names(stopbound)) {
      liste_patients <- liste_patients %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_notox), by = c("nb_ana", "ttt", "tot_eff_ttt0" = "y_cont")) %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_eff), by = c("nb_ana", "ttt", "tot_notox_ttt0" = "y_cont"))
    } else {
      liste_patients <- liste_patients %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_notox), by = c("nb_ana", "tot_eff_ttt0" = "y_cont")) %>%
        dplyr::left_join(stopbound %>% dplyr::select(-seuil_eff), by = c("nb_ana", "tot_notox_ttt0" = "y_cont"))
    }

    # Decision rules for the n_sim trials
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


#' Threshold's optimization
#'
#' Optimization of lambda and gamma in threshold Cn = 1 - &lambda; * (n / N) ^ &gamma;.
#'
#' 2 situations:
#' \itemize{
#'   \item for 1 arm: screening of the couples &lambda;/&gamma; by simulating monoarm trials under inefficacy/toxicity hypothesis to assess the control of &alpha; risk; and under
#'   the efficacy/non toxicity hypothesis to assess the power;
#'   \item for multiple arms with Bonferroni correction: same as for 1 arm with a risk of &alpha;/n_bras;
#'   \item for multiple arms without Bonferroni correction: screening of the couples &lambda;/&gamma; by simulating multiarm trials under Global Null Hypothesis to assess the control of FWER; and under
#'   the Least Favourable Configuration to assess the power.
#' }
#'
#' @return
#' Displaying in console of the result.
#' If stocked in an object, a list:
#' \itemize{
#'   \item the optimal couple &lambda;/&gamma; and &alpha;-risk and power a priori;
#'   \item the matrix of the choices.
#' }
#'
#' @param alpha Maximal type I error rate.
#' @param n_bras Number of treatment arms.
#' @param bonf TRUE to apply a Bonferroni correction.
#' @param nsim_oc Number of simulated trials.
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value).
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param p_a Law of probability under the efficacy/non toxicity hypothesis.
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param cut_seq,power_seq Vectors of values used to choose the couple while optimizing the Cn.
#' @param seed Seed for \code{rmultinom}.
#' @param methode Method used to generate the data:
#' \itemize{
#'   \item 1: trial by trial;
#'   \item 2: trial by trial, patient by patient;
#'   \item 3: whole in 1.
#' }
#' @param affich_mat Yes to display in console the matrix used for the choice, No for no display and default option NULL will ask you.
#'
#' @encoding UTF-8
#'
#' @examples
#' # The trial mimicked in the article of Zhou et al., 2017
#' carac_opt <- deter_cutoff(alpha = 0.1,
#'                           ana_inter = c(10, rep(5, 6)),
#'                           p_n = c(0.15, 0.3, 0.15, 0.4),
#'                           p_a = c(0.18, 0.42, 0.02, 0.38),
#'                           methode = 3L,
#'                           affich_mat = "No")
#' carac_opt[[1]]
#'
#' # A trial with 3 arms with the same hypotheses
#' carac_opt <- deter_cutoff(alpha = 0.1,
#'                           ana_inter = rep(15, 4),
#'                           p_n = c(0.1, 0.2, 0.2, 0.5),
#'                           p_a = c(0.15, 0.5, 0.05, 0.3),
#'                           n_bras = 3,
#'                           affich_mat = "No",
#'                           bonf = FALSE)
#' carac_opt[[1]]
#'
#' @export
deter_cutoff <- function(alpha = .1,
                         n_bras = 1,
                         bonf = FALSE,
                         nsim_oc = 10000,
                         ana_inter,
                         ana_inter_tox = NULL,
                         rand_ratio = NULL,
                         p_n,
                         prior = NULL,
                         p_a,
                         phi = NULL,
                         delta = NULL,
                         mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                         cut_seq = seq(.5, 1, by = .005),
                         power_seq = seq(0, 1, by = .01),
                         seed = 1024,
                         methode = 2L,
                         affich_mat = NULL) {

  # Check arguments
  if (alpha > 1 | alpha < 0)
    stop("Alpha-risk should be a real between 0 and 1.", call. = FALSE)
  if (alpha < 0.01 | alpha > 0.3)
    warning(paste0("The specified alpha-risk is a little extreme: ", alpha, "."), call. = FALSE, immediate. = TRUE)
  if (n_bras <= 0 | n_bras %% 1 != 0)
    stop("Number of arms should be an integer >= 1.", call. = FALSE)
  if (!is.null(affich_mat)) affich_mat <- match.arg(affich_mat, c("Yes", "yes", "No", "no"))
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
  if (!dplyr::near(sum(p_n), 1))
    stop("Vector \"p_n\" should sum to 1.", call. = FALSE)
  if (length(p_n) != 4)
    stop("\"p_n\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (is.null(prior)) {
    prior <- p_n
    message("\"p_n\" taken as prior.")
  } else {
    if (length(prior) != 4)
      stop("Vector \"prior\" should be of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
    if (!dplyr::near(sum(prior), 1))
      stop("Vector \"prior\" should sum to 1.", call. = FALSE)
  }
  if (!dplyr::near(sum(p_a), 1))
    stop("Vector \"p_a\" should sum to 1.", call. = FALSE)
  if (length(p_a) != 4)
    stop("\"p_a\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (length(nsim_oc) != 1 || !is.numeric(nsim_oc) || nsim_oc < 0 || nsim_oc %% 1 != 0)
    stop("\"nsim_oc\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.", call. = FALSE)
  }

  if (is.null(delta)) {
    message("Optimization vs a reference value.")
    # Determine phi. If no specified value, value under H0.
    # Will be the thresholds for efficacy and toxicity. Careful in the following, sometimes we take phi of non toxicity to be close to 0 and have more precision.
    if (is.null(phi)) phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  } else {
    message("Optimization vs control arm.")
    if (any(delta > 1) | any(delta < -1))
      stop("\"delta\" should be a real between -1 and 1.", call. = FALSE)
  }

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
      affich_mat <- readline(prompt = "Do you want to display the C_/gamma matrix? (Yes/No) > ")
      if (!stringr::str_detect(affich_mat, "Yes|yes|No|no")) {
        affich_mat <- NULL
        message("Answer with Yes or No.")
      }
    }
    return(affich_mat)
  }

  if (is.null(affich_mat)) affich_mat <- prompt_affich()

  # Bonferroni corection for multiarm settings
  if (bonf & n_bras > 1) {
    alpha <- alpha / n_bras
    n_bras <- 1
  }

  if (n_bras > 1) {
    # Reconstitution of LFC
    multinom_ttt_lfc <- matrix(c(p_a, rep(p_n, n_bras - 1)), nrow = 4) %>%
      as.data.frame() %>%
      as.list()
    names(multinom_ttt_lfc) <- paste0("ttt", seq_len(n_bras))

    # Reconstitution of Global Null Hypothesis
    multinom_ttt_H0 <- matrix(rep(p_n, n_bras), nrow = 4) %>%
      as.data.frame() %>%
      as.list()
    names(multinom_ttt_H0) <- paste0("ttt", seq_len(n_bras))
  }

  debut_cut <- 1 # Starting point for gamma
  # Because the Cn function is monotonic with lambda and gamma, no need to evaluate all couples
  matrice_carac <- matrix(numeric(0), ncol = 4, nrow = length(power_seq) + 2)

  # Monoarm BOP2
  if (n_bras == 1) {
    # Generation of datasets under LFC and GNH
    suppressMessages(
      tab_h0 <- gen_patients_multinom(
        n_sim = nsim_oc,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        multinom_ttt = list(p_n),
        multinom_cont = if (is.null(delta)) {NULL} else {p_n},
        seed = seed,
        methode = methode)
    )
    suppressMessages(suppressWarnings(
      tab_h1 <- gen_patients_multinom(
        n_sim = nsim_oc,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        multinom_ttt = list(p_a),
        multinom_cont = if (is.null(delta)) {NULL} else {p_n},
        seed = seed,
        methode = methode)
    ))

  } else { # Multiarm BOP2

    # Generation of datasets under LFC and GNH
    suppressMessages(
      tab_h0 <- gen_patients_multinom(
        n_sim = nsim_oc,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        multinom_ttt = multinom_ttt_H0,
        multinom_cont = if (is.null(delta)) {NULL} else {p_n},
        seed = seed,
        methode = methode)
    )
    suppressMessages(suppressWarnings(
      tab_h1 <- gen_patients_multinom(
        n_sim = nsim_oc,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        multinom_ttt = multinom_ttt_lfc,
        multinom_cont = if (is.null(delta)) {NULL} else {p_n},
        seed = seed,
        methode = methode)
    ))

  }

  # Loop over lambda and gamma
  for (j in seq_along(cut_seq)) {
    for (k in debut_cut:length(power_seq)) {
      oc_tox_n <- getoc_tox(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        C_ = cut_seq[j],
        gamm = power_seq[k],
        liste_patients = tab_h0,
        prior = prior,
        p_n = p_n,
        phi =  phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      if (n_bras == 1) { # Type I error rate is crossed so stop
        if (oc_tox_n$nb_accept > alpha + 0.00) break
      } else {
        if (oc_tox_n$rejet_glob > alpha + 0.00) break
      }

      oc_tox_a <- getoc_tox(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        C_ = cut_seq[j],
        gamm = power_seq[k],
        liste_patients = tab_h1,
        prior = prior,
        p_n = p_n,
        phi = phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      if (n_bras == 1) {
        matrice_carac[k,] <- c(cut_seq[j], power_seq[k], oc_tox_n$nb_accept, oc_tox_a$nb_accept)
      } else {
        matrice_carac[k,] <- c(cut_seq[j], power_seq[k], oc_tox_n$rejet_glob, oc_tox_a$rejet_ttt1)
      }

    }
    debut_cut <- k # Restart of the loop to last value of gamma evaluated
    if (debut_cut == length(power_seq)) {
      matrice_carac <- matrice_carac[!is.na(matrice_carac[, 1]),]
      break
    }
  }

  # Column 4 = P(Reject of H0 | H1) and column 3 = P(Reject H0 | H0)
  # Determine the optimal couple
  carac_optimale <- matrice_carac[matrice_carac[, 4] == max(matrice_carac[, 4]), ]
  if (is.matrix(carac_optimale)) carac_optimale <- carac_optimale[1, ]

  C_ <- carac_optimale[1]
  gamma <- carac_optimale[2]
  alpha_calc <- carac_optimale[3]
  puissance_calc <- carac_optimale[4]

  cat(paste0("Optimal couple is: lambda = ", C_, " and gamma = ", gamma, ".\n"))
  cat(paste0("Calculated alpha-risk is ", alpha_calc, " and power is ", puissance_calc, ".\n"))
  if (affich_mat %in% c("Yes", "yes")) print(matrice_carac)

  invisible(list(
    c(C_ = C_,
      gamma = gamma,
      alpha_calc = alpha_calc,
      puissance_calc = puissance_calc),
    mat = matrice_carac
  ))

}
