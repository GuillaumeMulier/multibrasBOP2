
#' Table of the stopping boundaries for BOP with thresholds dependant on the number of active arms
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
#' of comparison to a reference value). Not yet implemented.
#' @param C_holm,gamm_holm The 2 tuning parameters of criterion Cn optimized with the correction dependant on the number of active arms.
#' @param C_mono,gamm_mono The 2 tuning parameters of criterion Cn optimized for a single arm design.
#' @param m The maximum number of active arms.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. If no control group, leave NULL. Should be between -1 and 1.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # Example for uncontrolled design
#' get_stopbound_holm2(ana_inter = rep(15, 4),
#'                     C_ = 0.535, gamm = 0.98, m = 3,
#'                     p_n = c(0.15, 0.3, 0.15, 0.4))
#'
#' # Example for controlled design
#' get_stopbound_holm2(ana_inter = rep(15, 4),
#'                     C_ = 0.65, gamm = 0.8, m = 3,
#'                     p_n = c(0.15, 0.3, 0.15, 0.4),
#'                     delta = c(0, 0))
#'
#' @export
get_stopbound_holm2 <- function(ana_inter,
                                ana_inter_tox = NULL,
                                rand_ratio = NULL,
                                C_mono, gamm_mono,
                                C_holm, gamm_holm,
                                m,
                                p_n,
                                prior = NULL,
                                mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                phi = NULL,
                                delta = NULL) {

  if (is.null(phi))
    phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  if (is.null(rand_ratio)) {
    message("Without specified \"rand_ratio\" all groups will have the same size (same as control group).")
    rand_ratio <- 1
  } else if (!is.numeric(rand_ratio)) {
    stop("\"rand_ratio\" should be a number or a vector of numbers.", call. = FALSE)
  } else if (any(rand_ratio <= 0)) {
    stop("\"rand_ratio\" must be positive.", call. = FALSE)
  }
  if (is.null(prior)) {
    prior <- p_n
  } else {
    if (length(prior) != 4)
      stop("\"prior\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
    if (!dplyr::near(sum(prior), 1))
      stop("\"prior\" should sum to 1 as it is a law of probability.", call. = FALSE)
  }
  if (length(unique(rand_ratio)) == 1 & length(rand_ratio) > 1) rand_ratio <- rand_ratio[1]

  if (!is.null(ana_inter_tox)) {
    if (sum(ana_inter) != sum(ana_inter_tox))
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
          tab <- list()
          for (i in seq_len(m + 1) - 1) {
            nu <- min(1 + m - i, m)
            seuil <- (1 - C_holm / nu) * (n / N) ^ gamm_holm # Cutoff from Ruitao to take into account the number of active arms
            # seuil <- 1 - (1 - C_holm / nu) * (n / N) ^ gamm_holm # Cutoff from Ruitao to take into account the number of active arms
            # On aurait un seuil CN qui vaut 1 - (nu - C_holm) / nu * (n / N) ^ gamma_holm
            if (n == N) seuil <- max(seuil, C_mono) # At the end of the study, to ensure the control of type I error rate take a more conservative cutoff including BOP2 threshold for single arm
            if (i != 0) {
              tab[[match(i, seq_len(m + 1) - 1)]] <- c(i, max_resp(
                cutoff = seuil,
                prior = prior,
                mat_beta_xi = mat_beta_xi,
                ntr = n,
                phi = phi
              ))
            } else {
              tab[[match(i, seq_len(m + 1) - 1)]] <- c(i, -1, -1)
            }
          }
          tab <- do.call("rbind", tab)
          return(tab)
        }
      )
      stopbound <- do.call("rbind", stopbound)
      stopbound <- stopbound %>%
        as.data.frame() %>%
        dplyr::mutate(
          nb_ana = rep(seq_along(ana_inter_ttt), each = m + 1),
          tot_pat = rep(ana_inter_cum, each = m + 1),
          .before = 1
        ) %>%
        dplyr::mutate(seuil_tox = tot_pat - V3)
      names(stopbound)[3:5] <- c("nb_act", "seuil_eff", "seuil_notox")
      stopbound$seuil_eff[stopbound$tot_pat %nin% ana_eff] <- -1
      stopbound$seuil_notox[stopbound$tot_pat %nin% ana_tox] <- -1
      stopbound$seuil_tox[stopbound$tot_pat %nin% ana_tox] <- stopbound$tot_pat[stopbound$tot_pat %nin% ana_tox] + 1

    } else { # Vs control group

      ana_inter_cum_cont <- cumsum(anas_inters)
      N0 <- sum(anas_inters)
      stopbound <- purrr::pmap(
        .l = list(n0 = ana_inter_cum_cont, nt = ana_inter_cum),
        .f = function(n0, nt) {
          tab <- list()
          for (i in seq_len(m)) {
            nu <- min(1 + m - i, m)
            seuil <- ((nu - C_holm) / nu) * ((n0 + nt) / (N0 + N)) ^ gamm_holm # Cutoff from Ruitao to take into account the number of active arms
            # seuil <- 1 - (1 - C_holm / nu) * ((n0 + nt) / (N0 + N)) ^ gamm_holm # Cutoff from Ruitao to take into account the number of active arms
            if (nt == N) seuil <- max(seuil, C_mono) # At the end of the study, to ensure the control of type I error rate take a more conservative cutoff including BOP2 threshold for single arm
            tab[[i]] <- cbind(
              nb_act = i,
              max_resp(
                cutoff = seuil,
                prior = prior,
                mat_beta_xi = mat_beta_xi,
                ntr = nt,
                ntr_cont = n0,
                phi = phi,
                delta = delta
              )
            )
          }
          tab <- do.call("rbind", tab)
          return(tab)
        }
      )
      stopbound <- lapply(stopbound,
                          function(x) {
                            x %>%
                              dplyr::group_by(nb_act, y_cont) %>%
                              dplyr::summarise(seuil_eff = sum(deci_eff) - 1,
                                               seuil_notox = sum(deci_notox) - 1,
                                               .groups = "drop") %>%
                              dplyr::bind_rows(tibble::tibble(nb_act = 0,
                                                              y_cont = seq(0, max(.$y_cont)),
                                                              seuil_eff = -1,
                                                              seuil_notox = -1))
                          })
      stopbound <- lapply(seq_along(stopbound), function(x) cbind(stopbound[[x]], nb_ana = x))
      stopbound <- do.call("rbind", stopbound) %>%
        dplyr::mutate(tot_pat = ana_inter_cum_cont[match(nb_ana, seq_along(ana_inter_cum_cont))],
                      tot_pat_ttt = ana_inter_cum[match(nb_ana, seq_along(ana_inter_cum))])
      stopbound$seuil_eff[stopbound$tot_pat %nin% ana_eff] <- -1
      stopbound$seuil_notox[stopbound$tot_pat %nin% ana_tox] <- -1
    }

  } else { # Groups without the same number of patients
    # Not implemented yet

    # ana_inter_ttt <- lapply(rand_ratio, function(x) anas_inters * x)
    # ana_eff_ttt <- lapply(rand_ratio, function(x) ana_eff * x)
    # names(ana_eff_ttt) <- paste0("ttt", seq_along(ana_eff_ttt))
    # ana_tox_ttt <- lapply(rand_ratio, function(x) ana_tox * x)
    # names(ana_tox_ttt) <- paste0("ttt", seq_along(ana_tox_ttt))
    # if (lapply(ana_inter_ttt, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
    #   warning("\"rand_ratio\" x \"ana_inter\" isn't always an integer.
    #           Converted to integers with function ceiling.",
    #           call. = FALSE, immediate. = TRUE)
    #   ana_inter_ttt <- lapply(ana_inter_ttt,  ceiling)
    #   ana_eff_ttt <- lapply(ana_eff_ttt,  ceiling)
    #   ana_tox_ttt <- lapply(ana_tox_ttt,  ceiling)
    # }
    # if (is.null(names(ana_inter_ttt))) names(ana_inter_ttt) <- paste0("ttt", seq_along(ana_inter_ttt))
    #
    # N <- lapply(ana_inter_ttt, sum)
    # ana_inter_cum <- lapply(ana_inter_ttt, cumsum)
    #
    # if (is.null(delta)) { # Vs a reference value
    #
    #   stopbound <- purrr::pmap(
    #     .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt, ana_eff_ttt = ana_eff_ttt, ana_tox_ttt = ana_tox_ttt),
    #     .f = function(.x, .y, .z, ana_eff_ttt, ana_tox_ttt) {
    #       threshold <- lapply(
    #         X = .x,
    #         FUN = function(n) {
    #           max_resp(
    #             cutoff = C_ * (n / .y) ^ gamm,
    #             prior = prior,
    #             mat_beta_xi = mat_beta_xi,
    #             ntr = n,
    #             phi = phi
    #           )
    #         }
    #       )
    #       threshold <- do.call("cbind", threshold)
    #       threshold <- threshold %>%
    #         t() %>%
    #         as.data.frame() %>%
    #         dplyr::mutate(nb_ana = seq_along(.z),
    #                       tot_pat = .x,
    #                       .before = 1) %>%
    #         dplyr::mutate(seuil_tox = tot_pat - V2)
    #       names(threshold)[3:4] <- c("seuil_eff", "seuil_notox")
    #       threshold$seuil_eff[threshold$tot_pat %nin% ana_eff_ttt] <- -1
    #       threshold$seuil_notox[threshold$tot_pat %nin% ana_tox_ttt] <- -1
    #       threshold$seuil_tox[threshold$tot_pat %nin% ana_tox_ttt] <- threshold$tot_pat[threshold$tot_pat %nin% ana_tox_ttt] + 1
    #       return(threshold)
    #     }
    #   )
    #
    #   stopbound <- purrr::map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))
    #
    # } else { # Vs control group
    #
    #   ana_inter_cum_cont <- cumsum(anas_inters)
    #   N0 <- sum(anas_inters)
    #
    #   stopbound <- purrr::pmap(
    #     .l = list(.x = ana_inter_cum, .y = N, .z = ana_inter_ttt, ana_eff_ttt = ana_eff_ttt, ana_tox_ttt = ana_tox_ttt),
    #     .f = function(.x, .y, .z, ana_eff_ttt, ana_tox_ttt) {
    #       threshold <- purrr::pmap(
    #         .l = list(n0 = ana_inter_cum_cont, nt = .x),
    #         .f = function(n0, nt) {
    #           max_resp(
    #             cutoff = C_ * ((n0 + nt) / (N0 + .y)) ^ gamm,
    #             prior = prior,
    #             mat_beta_xi = mat_beta_xi,
    #             ntr = nt,
    #             ntr_cont = n0,
    #             phi = phi,
    #             delta = delta
    #           )
    #         }
    #       )
    #       threshold <- lapply(threshold,
    #                           function(x) {
    #                             x %>%
    #                               dplyr::group_by(y_cont) %>%
    #                               dplyr::summarise(
    #                                 seuil_eff = sum(deci_eff) - 1,
    #                                 seuil_notox = sum(deci_notox) - 1
    #                               )
    #                           })
    #       threshold <- lapply(seq_along(threshold), function(x) cbind(threshold[[x]], nb_ana = x))
    #       threshold <- do.call("rbind", threshold) %>%
    #         dplyr::mutate(tot_pat = ana_inter_cum_cont[match(nb_ana, seq_along(ana_inter_cum_cont))],
    #                       tot_pat_ttt = .x[match(nb_ana, seq_along(.x))])
    #       threshold$seuil_eff[threshold$tot_pat_ttt %nin% ana_eff_ttt] <- -1
    #       threshold$seuil_notox[threshold$tot_pat_ttt %nin% ana_tox_ttt] <- -1
    #       return(threshold)
    #     }
    #   )
    #
    #   stopbound <- purrr::map_dfr(seq_along(stopbound), ~ cbind(ttt = names(stopbound)[.x], stopbound[[.x]]))
    #
    # }
    #
  }

  return(stopbound)

}



#' Operating characteristics for a couple &lambda_holm;/&gamma_holm; given a trials' list and the same hyperparameters for a monoarm trial.
#'
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value). Not yet implemented.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param C_holm,gamm_holm The 2 tuning parameters of criterion Cn optimized with the correction dependant on the number of active arms.
#' @param C_mono,gamm_mono The 2 tuning parameters of criterion Cn optimized for a single arm design.
#' @param m The maximum number of active arms.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints (efficacy and non toxicity /!\).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param liste_patients Generated trials via function \code{gen_patients_multinom}.
#'
#' @importFrom magrittr %>%
#'
#' @encoding UTF-8
getoc_tox_holm2 <- function(ana_inter,
                            ana_inter_tox = NULL,
                            rand_ratio = NULL,
                            prior = NULL,
                            p_n,
                            m,
                            C_holm, gamm_holm,
                            C_mono, gamm_mono,
                            mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                            phi = NULL,
                            delta = c(0, 0),
                            liste_patients) {

  if (is.null(prior)) prior <- p_n

  if (is.null(phi)) phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))

  ana_eff         <- cumsum(ana_inter)
  ana_tox         <- cumsum(ana_inter_tox)
  anas_inters_cum <- sort(union(ana_eff, ana_tox))
  anas_inters     <- c(anas_inters_cum[1], diff(anas_inters_cum))

  if ("ttt0" %nin% unique(liste_patients$ttt)) { # Vs a reference value

    # Generation of stopbound table at each interim analysis and for each group
    stopbound <- suppressWarnings(
      get_stopbound_holm2(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        m = m,
        C_mono = C_mono,
        gamm_mono = gamm_mono,
        C_holm = C_holm,
        gamm_holm = gamm_holm,
        prior = prior,
        p_n = p_n,
        mat_beta_xi = mat_beta_xi,
        phi = phi
      )
    )

    # Splitting the data to perform separately each trial
    liste_patients <- split(liste_patients, liste_patients$n_simu)
    liste_patients <- purrr::map(
      .x = liste_patients,
      .f = function(data) {
        data$nb_act <- NA_integer_
        data$continuer <- NA
        data$seuil_eff <- data$seuil_notox <- NA_real_
        for (i in seq_along(anas_inters_cum)) {
          if (i == 1) data$nb_act[data$nb_ana == i] <- bras_act <- m else data$nb_act[data$nb_ana == i] <- bras_act <- sum(data$continuer[data$nb_ana == (i - 1)])
          data$seuil_eff[data$nb_ana == i] <- stopbound$seuil_eff[stopbound$nb_ana == i & stopbound$nb_act == bras_act]
          data$seuil_notox[data$nb_ana == i] <- stopbound$seuil_notox[stopbound$nb_ana == i & stopbound$nb_act == bras_act]
          if (i == 1) {
            data$continuer[data$nb_ana == i] <- data$tot_eff[data$nb_ana == i] > data$seuil_eff[data$nb_ana == i] & data$tot_notox[data$nb_ana == i] > data$seuil_notox[data$nb_ana == i]
          } else {
            data$continuer[data$nb_ana == i] <- data$continuer[data$nb_ana == (i - 1)]
            data$continuer[data$continuer & data$nb_ana == i] <- data$tot_eff[data$continuer & data$nb_ana == i] > data$seuil_eff[data$continuer & data$nb_ana == i] &
              data$tot_notox[data$continuer & data$nb_ana == i] > data$seuil_notox[data$continuer & data$nb_ana == i]
          }
        }
        return(data)
      }
    ) %>%
      do.call(what = "rbind", args = .)

    # Decision rules
    liste_patients <- liste_patients %>%
      dplyr::mutate(decision = dplyr::case_when((nb_ana != max(nb_ana)) & continuer ~ "Continue",
                                                (nb_ana == max(nb_ana)) & continuer ~ "Acceptation traitement",
                                                (nb_ana != max(nb_ana)) & !continuer ~ "Arrêt précoce",
                                                (nb_ana == max(nb_ana)) & !continuer ~ "Arrêt")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

  } else { # Vs control group

    # Generation of stopbound table at each interim analysis and for each group
    stopbound <- suppressWarnings(
      get_stopbound_holm2(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        m = m,
        C_holm = C_holm,
        gamm_holm = gamm_holm,
        C_mono = C_mono,
        gamm_mono = gamm_mono,
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

    # Adding the stopping boundaries and performing the trials
    liste_patients <- split(liste_patients, liste_patients$n_simu)
    liste_patients <- purrr::map(
      .x = liste_patients,
      .f = function(data) {
        data$nb_act <- NA_integer_
        data$continuer <- NA
        data$seuil_eff <- data$seuil_notox <- NA_real_
        for (i in seq_along(anas_inters_cum)) {
          if (i == 1) data$nb_act[data$nb_ana == i] <- bras_act <- m else data$nb_act[data$nb_ana == i] <- bras_act <- sum(data$continuer[data$nb_ana == (i - 1)])
          data$seuil_eff[data$nb_ana == i] <- stopbound$seuil_eff[stopbound$nb_ana == i & stopbound$nb_act == bras_act & stopbound$y_cont == data$tot_eff_ttt0[data$nb_ana == i][1]]
          data$seuil_notox[data$nb_ana == i] <- stopbound$seuil_notox[stopbound$nb_ana == i & stopbound$nb_act == bras_act & stopbound$y_cont == data$tot_notox_ttt0[data$nb_ana == i][1]]
          if (i == 1) {
            data$continuer[data$nb_ana == i] <- data$tot_eff[data$nb_ana == i] > data$seuil_eff[data$nb_ana == i] & data$tot_notox[data$nb_ana == i] > data$seuil_notox[data$nb_ana == i]
          } else {
            data$continuer[data$nb_ana == i] <- data$continuer[data$nb_ana == (i - 1)]
            data$continuer[data$continuer & data$nb_ana == i] <- data$tot_eff[data$continuer & data$nb_ana == i] > data$seuil_eff[data$continuer & data$nb_ana == i] &
              data$tot_notox[data$continuer & data$nb_ana == i] > data$seuil_notox[data$continuer & data$nb_ana == i]
          }
        }
        return(data)
      }
    ) %>%
      do.call(what = "rbind", args = .)

    # Decision rules for the n_sim trials
    liste_patients <- liste_patients %>%
      dplyr::mutate(decision = dplyr::case_when((nb_ana != max(nb_ana)) & continuer ~ "Continue",
                                                (nb_ana == max(nb_ana)) & continuer ~ "Acceptation traitement",
                                                (nb_ana != max(nb_ana)) & !continuer ~ "Arrêt précoce",
                                                (nb_ana == max(nb_ana)) & !continuer ~ "Arrêt")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  # Output of the function
  rejet_glob <- liste_patients %>%
    dplyr::mutate(decision_h0 = decision == "Acceptation traitement") %>%
    dplyr::group_by(n_simu) %>%
    dplyr::summarise(decision_h0 = sum(decision_h0) > 0) %>%
    dplyr::pull(decision_h0) %>%
    mean()

  rejet_ttt1 <- liste_patients %>%
    dplyr::filter(ttt == "ttt1") %>%
    dplyr::mutate(decision = decision == "Acceptation traitement") %>%
    dplyr::pull(decision) %>%
    mean()

  return(list(rejet_glob = rejet_glob,
              rejet_ttt1 = rejet_ttt1))

}


#' Threshold's optimization
#'
#' Optimization of lambda and gamma in threshold Cn = 1 - (&eta; - &lambda;) / &eta; * (n / N) ^ &gamma;, with &eta; = K + 1 - k, K the number of arms (not counting control arm) and k the number of still active arms.
#' To ensure to control the type I error rate in the individual arms, at the final analysis, the threshold taken is min(CN, 1 - &lamba;_monoarm) ensuring a not too loose cutoff.
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
#' of comparison to a reference value). Not yet implemented.
#' @param p_n Law of probability under the inefficacy/toxicity hypothesis.
#' @param prior Prior for the law of probability (by default, skeptical: the inefficacy/toxicity hypothesis).
#' @param p_a Law of probability under the efficacy/non toxicity hypothesis.
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param cut_seq_holm,power_seq_holm Vectors of values used to choose the couple while optimizing the Cn.
#' @param cut_seq_mono,power_seq_mono Vectors of values used to choose the couple while optimizing the Cn for the monoarm trial. If only 1 number is supplied, won't optimise
#' but take the given value as hyperparameter.
#' @param m The maximum number of active arms.
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
#' if (FALSE) { # Not run because of time it can take
#' # A trial with 3 arms with the same hypotheses as article of Zhou et al. 2017
#' carac_opt <- deter_cutoff_holm2(alpha = 0.1,
#'                                 ana_inter = rep(15, 4),
#'                                 p_n = c(0.1, 0.2, 0.2, 0.5),
#'                                 p_a = c(0.15, 0.5, 0.05, 0.3),
#'                                 n_bras = 3,
#'                                 affich_mat = "No")
#' carac_opt[[1]]
#' }
#'
#' @export
deter_cutoff_holm2 <- function(alpha = .1,
                               n_bras = 3,
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
                               cut_seq_mono = seq(.5, 1, by = .005),
                               power_seq_mono = seq(0, 1, by = .01),
                               cut_seq_holm = seq(.5, 1, by = .005),
                               power_seq_holm = seq(0, 1, by = .01),
                               seed = 1024,
                               methode = 2L,
                               affich_mat = NULL) {

  if (is.null(prior)) {
    prior <- p_n
    message("\"p_n\" taken as prior.")
  }

  if (is.null(delta)) {
    message("Optimization vs a reference value.")
    # Determine phi. If no specified value, value under H0.
    # Will be the thresholds for efficacy and toxicity. Careful in the following, sometimes we take phi of non toxicity to be close to 0 and have more precision.
    if (is.null(phi)) phi <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  } else {
    message("Optimization vs control arm.")
    # if (any(delta > 1) | any(delta < -1))
  }

  if (is.null(rand_ratio)) {
    rand_ratio <- 1
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
  if (bonf) {
    alpha <- alpha / n_bras
    n_bras <- 1
  }

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

  # Compute threshold for monoarm BOP2
  if (length(cut_seq_mono) == 1 & length(power_seq_mono) == 1) {
    C_mono <- cut_seq_mono
    gamm_mono <- power_seq_mono
  } else {
    carac_mono <- deter_cutoff(
      alpha = alpha,
      n_bras = 1,
      bonf = FALSE,
      nsim_oc = nsim_oc,
      ana_inter = ana_inter,
      ana_inter_tox = ana_inter_tox,
      rand_ratio = rand_ratio,
      p_n = p_n,
      prior = prior,
      p_a = p_a,
      phi = phi,
      delta = delta,
      mat_beta_xi = mat_beta_xi,
      cut_seq = cut_seq_mono,
      power_seq = power_seq_mono,
      seed = seed,
      methode = methode,
      affich_mat = affich_mat
    )
    C_mono <- carac_mono[[1]][["C_"]]
    gamm_mono <- carac_mono[[1]][["gamma"]]
  }

  # Now determine the threshold's parameter for the new cutoff
  debut_cut <- 1 # Starting point for gamma
  # Because the Cn function is monotonic with lambda and gamma, no need to evaluate all couples
  matrice_carac <- matrix(numeric(0), ncol = 4, nrow = length(power_seq_holm) + 2)
  # matrice_carac <- matrix(numeric(0), ncol = 4, nrow = length(power_seq_holm) * length(cut_seq_holm))

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

  cut_seq_holm <- sort(cut_seq_holm, decreasing = TRUE)
  # Due to other form of Cn, Cn is now increasing with lambda, so in order to decrease alpha risk at next iteration we need
  # to lower lambda. Thus working in reverse order of lambda in the loop
  # Loop over lambda and gamma
  for (j in seq_along(cut_seq_holm)) {
    for (k in debut_cut:length(power_seq_holm)) {
      # for (k in 1:length(power_seq_holm)) {
      # print(paste0("C=", j, "/G=", k))
      oc_tox_n <- getoc_tox_holm2(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        m = n_bras,
        C_holm = cut_seq_holm[j],
        gamm_holm = power_seq_holm[k],
        C_mono = C_mono,
        gamm_mono = gamm_mono,
        liste_patients = tab_h0,
        prior = prior,
        p_n = p_n,
        phi =  phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      if (oc_tox_n$rejet_glob > alpha + 0.00) break

      oc_tox_a <- getoc_tox_holm2(
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        m = n_bras,
        C_holm = cut_seq_holm[j],
        gamm_holm = power_seq_holm[k],
        C_mono = C_mono,
        gamm_mono = gamm_mono,
        liste_patients = tab_h1,
        prior = prior,
        p_n = p_n,
        phi = phi,
        delta = delta,
        mat_beta_xi = mat_beta_xi
      )

      matrice_carac[k,] <- c(cut_seq_holm[j], power_seq_holm[k], oc_tox_n$rejet_glob, oc_tox_a$rejet_ttt1)
      # matrice_carac[ (j - 1) * length(power_seq_holm) + k,] <- c(cut_seq_holm[j], power_seq_holm[k], oc_tox_n$rejet_glob, oc_tox_a$rejet_ttt1)

    }
    debut_cut <- k # Restart of the loop to last value of gamma evaluated
    if (debut_cut == length(power_seq_holm)) break
  }

  matrice_carac <- matrice_carac[!is.na(matrice_carac[, 1]),]

  # Column 4 = P(Reject of H0 | H1) and column 3 = P(Reject H0 | H0)
  # Determine the optimal couple
  carac_optimale <- matrice_carac[matrice_carac[, 4] == max(matrice_carac[, 4]), ]
  if (is.matrix(carac_optimale)) carac_optimale <- carac_optimale[1, ]

  C_holm <- carac_optimale[1]
  gamma_holm <- carac_optimale[2]
  alpha_calc <- carac_optimale[3]
  puissance_calc <- carac_optimale[4]

  cat(paste0("Optimal couple is: lambda = ", C_holm, " and gamma = ", gamma_holm, ".\n"))
  cat(paste0("Optimal couple for mono-arm BOP2 is: lambda = ", C_mono, " and gamma = ", gamm_mono, ".\n"))
  cat(paste0("Calculated alpha-risk is ", alpha_calc, " and power is ", puissance_calc, ".\n"))
  if (affich_mat %in% c("Yes", "yes")) print(matrice_carac)

  invisible(list(
    c(C_holm = C_holm,
      gamma_holm = gamma_holm,
      C_mono = C_mono,
      gamma_mono = gamm_mono,
      alpha_calc = alpha_calc,
      puissance_calc = puissance_calc),
    mat = matrice_carac
  ))

}
