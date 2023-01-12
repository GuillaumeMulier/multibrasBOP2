#' Operating characteristics
#'
#' Obtain the operating characteristics with a constant threshold vs a reference value (or control group).
#'
#' Optimize the couple gamma_eff/gamma_tox, then simulate \code{nsim_essais} trials under \code{p_reel} in order to determine the operating characteristics.
#' It is possible to optimise the couple before with the function \code{deter_constcutoff} and supply the obtained values in arguments
#' \code{seq_eff} and \code{seq_tox}.
#' The stopping rule can be written as stop if Pr(mat_beta_xi * p_ttt <= phi | Dn) > 1 - gamma_eff or 1 - gamma_tox.
#'
#' Works for mono- and multiarm trials.
#'
#' For multiarm settings, \code{bonf = TRUE} to perform Bonferroni correction (less time, but more conservative or even couple).
#'
#' @return
#' \itemize{
#'   \item for monoarm, a list:
#'   \itemize{
#'     \item caracteristique: data.frame with the global informations (type I error rate, inefficacy/toxicity and efficacy/non toxicity hypotheses,
#'     simulation probabilities, method used to generate the data,
#'     Cn, proportion of reject of inefficacy/toxicity hypothesis and early stoppings, early stoppings for futility/toxicity, mean number of stoppings (for
#'     futility and toxicity), mean number of enrolled patients, responses and toxicities);
#'     \item essais: data.frame of the simulated trials with 1 line for each trial (the final analysis).
#'   }
#'   \item for multiarm, a list:
#'   \itemize{
#'     \item carac_globales: data.frame with the global informations (FWER, inefficacy/toxicity and efficacy/non toxicity hypotheses,
#'     simulation probabilities, method used to generate the data, Cn, proportion of global reject of inefficacy/toxicity hypothesis);
#'     \item carac_bras: data.frame with the informations for each arm (proportion of reject of the inefficacy/toxicity hypothesis,
#'     mean number of patients/responses/toxicities in each arm, proportion of early stoppings (and for futility/toxicity) and stoppings
#'     (for futility and toxicity too));
#'     \item essais: data.frame of the simulated trials.
#'   }
#' }
#'
#' @encoding UTF-8
#'
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param alpha Maximal type I error rate (FWER in case of multiarm).
#' @param bonf TRUE to apply a Bonferroni correction in case of multiarm settings.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param p_reel,p_n,p_a Vectors of length 4 giving the law of probability under the 'real' scenario, the inefficacy/toxicity hypothesis and the efficacy/non toxicity hypothesis.
#' @param nsim_oc,nsim_essai Number of simulated trials to optimize CN and get the operating characteristics.
#' @param tableau_essais Data.frame of simulated trials with \code{gen_patients_multinom}.
#' @param seq_eff,seq_tox Vectors of values used to choose the couple while optimizing the couple gamma_eff/gamma_tox. If a unique value is specified for each parameter, the grid search won't be performed.
#' @param seed_cut,seed_pat Seeds for \code{rmultinom} (for Cn and operating characteristics).
#' @param methode Method used to generate the data:
#' \itemize{
#'   \item 1: trial by trial;
#'   \item 2: trial by trial, patient by patient;
#'   \item 3: whole in 1.
#' }
#'
#' @examples
#' # For BOP2 simulated trials in Zhou, 2017
#' opcharac_efftox_const(ana_inter = c(10, rep(5, 6)),
#'                       alpha = 0.1,
#'                       p_reel = c(0.18, 0.42, 0.02, 0.38),
#'                       p_n = c(0.15, 0.3, 0.15, 0.4), p_a = c(0.18, 0.42, 0.02, 0.38),
#'                       nsim_oc = 10000, nsim_essai = 10000,
#'                       seed_cut = 1024, seed_pat = 1993,
#'                       methode = 3L)
#'
#' # For a 3 arm-trial
#' opcharac_efftox(ana_inter = rep(15, 4),
#'                 alpha = 0.1,
#'                 p_reel = list(ttt1 = c(0.15, 0.5, 0.05, 0.3), ttt2 = c(0.15, 0.3, 0.25, 0.3), ttt3 = c(0.05, 0.15, 0.15, 0.65)),
#'                 p_n = c(0.1, 0.2, 0.2, 0.5), p_a = c(0.15, 0.5, 0.05, 0.3),
#'                 nsim_oc = 10000, nsim_essai = 10000)
#'
#' @export
opcharac_efftox_const <- function(ana_inter,
                                  alpha = 0.1,
                                  mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                  phi = NULL,
                                  delta = NULL,
                                  p_reel,
                                  bonf = FALSE,
                                  p_n,
                                  p_a,
                                  nsim_oc = 10000,
                                  nsim_essai = 10000,
                                  tableau_essais = NULL,
                                  seq_eff = seq(.01, .99, by = .02),
                                  seq_tox = seq(.01, .99, by = .02),
                                  seed_cut = 1024,
                                  seed_pat = 1993,
                                  methode = 2L) {
  # Check arguments
  if (alpha > 1 | alpha < 0)
    stop("Alpha-risk should be a real between 0 and 1.", call. = FALSE)
  if (alpha < 0.01 | alpha > 0.3)
    warning(paste0("The specified alpha-risk is a little extreme: ", alpha, "."), call. = FALSE, immediate. = TRUE)
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (!dplyr::near(sum(p_n), 1))
    stop("Vector \"p_n\" should sum to 1.", call. = FALSE)
  if (length(p_n) != 4)
    stop("\"p_n\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (!is.list(p_reel))
    p_reel <- list(ttt1 = p_reel)
  if (any(lapply(p_reel, length) != 4))
    stop("Vectors \"p_reel\" shoud be of length 4 for each treatment:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (any(unlist(lapply(lapply(p_reel, sum),
                        FUN = function(x) !dplyr::near(x, 1)))))
    stop("Vectors \"p_reel\" should sum to 1 for each treatment.", call. = FALSE)
  if (!dplyr::near(sum(p_a), 1))
    stop("Vector \"p_a\" should sum to 1.", call. = FALSE)
  if (length(p_a) != 4)
    stop("\"p_a\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (length(nsim_oc) != 1 || !is.numeric(nsim_oc) || nsim_oc < 0 || nsim_oc %% 1 != 0)
    stop("\"nsim_oc\" should be a positive integer.", call. = FALSE)
  if (length(nsim_essai) != 1 || !is.numeric(nsim_essai) || nsim_essai < 0 || nsim_essai %% 1 != 0)
    stop("\"nsim_essai\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)

  # Reconstitution of phi. If no specified value, the values under H0 are taken.
  if (is.null(phi)) {
    phitox <- c(sum(p_n * mat_beta_xi[1, ]), 1 - sum(p_n * mat_beta_xi[2, ]))
    phi <- phitox
    phinotox <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  } else {
    phitox <- phi
    phinotox <- c(phi[1], 1 - phi[2])
  }

  # Number of arms
  n_bras <- length(p_reel)

  # Optimization of couple gamma_eff/gamma_tox
  if (length(seq_eff) == 1 & length(seq_tox) == 1) {
    gamma_eff <- seq_eff
    gamma_tox <- seq_tox
  } else {
    carac_optimales <- deter_constcutoff(
      alpha = alpha,
      n_bras = n_bras,
      bonf = bonf,
      ana_inter = ana_inter,
      mat_beta_xi = mat_beta_xi,
      phi = phinotox,
      delta = delta,
      p_n = p_n,
      p_a = p_a,
      nsim_oc = nsim_oc,
      seq_eff = seq_eff,
      seq_tox = seq_tox,
      seed = seed_cut,
      methode = methode,
      affich_mat = "Non"
    )
    gamma_eff <- carac_optimales[[1]]["gamma_eff"] # Store the 2 thresholds
    gamma_tox <- carac_optimales[[1]]["gamma_tox"]
  }

  nmax <- sum(ana_inter)

  # resp_list give the number of responses and toxicities at each interim analysis for each trial
  if (is.null(tableau_essais)) {
    resp_list <- gen_patients_multinom(
      n_sim = nsim_essai,
      ana_inter = ana_inter,
      multinom_ttt = p_reel,
      multinom_cont = if (is.null(delta)) {NULL} else {p_n},
      seed = seed_pat,
      methode = methode
    )
  } else {
    resp_list <- tableau_essais
  }

  phi_eff <- phitox[1] # Phi for efficacy
  phi_tox <- phitox[2] # Phi for toxicity
  prior_eff <- sum(p_n * mat_beta_xi[1,])
  prior_tox <- sum(p_n * (1 - mat_beta_xi[2,]))

  # Stopping boundaries
  stopbound <- get_stopbound_const(
    ana_inter = ana_inter,
    seuils = c(gamma_eff, gamma_tox),
    p_n = p_n,
    phi = phinotox,
    delta = delta
  )

  if ("ttt0" %nin% unique(resp_list$ttt)) { # Vs a reference value

    # Adding the stopping boundaries to resp_list
    resp_list <- resp_list %>%
      dplyr::mutate(nb_ana = as.integer(nb_ana)) %>%
      dplyr::left_join(stopbound, by = c("nb_ana", "tot_pat"))

    # Decision rules
    resp_list <- resp_list %>%
      dplyr::mutate(
        prior_eff = prior_eff,
        prior_tox = prior_tox,
        decision = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Continue",
                                    (nb_ana == max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Accept the treatment",
                                    (nb_ana != max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Early stopping",
                                    (nb_ana == max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Stopping"),
        decision_eff = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff ~ "Continue",
                                        (nb_ana == max(nb_ana)) & tot_eff > seuil_eff ~ "Accept the treatment",
                                        (nb_ana != max(nb_ana)) & tot_eff <= seuil_eff ~ "Early stopping (futility)",
                                        (nb_ana == max(nb_ana)) & tot_eff <= seuil_eff ~ "Stopping (futility)"),
        decision_tox = dplyr::case_when((nb_ana != max(nb_ana)) & tot_notox > seuil_notox ~ "Continue",
                                        (nb_ana == max(nb_ana)) & tot_notox > seuil_notox ~ "Accept the treatment",
                                        (nb_ana != max(nb_ana)) & tot_notox <= seuil_notox ~ "Early stopping (toxicity)",
                                        (nb_ana == max(nb_ana)) & tot_notox <= seuil_notox ~ "Stopping (toxicity)")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

  } else { # Vs a control group

    # Wide format to compare control and treatment groups
    resp_list <- resp_list %>%
      tidyr::pivot_wider(names_from = ttt, values_from = efftox:tot_notox) %>%
      tidyr::pivot_longer(cols = tidyselect::matches("_ttt[^0]$"), names_pattern = "^(.*)_(ttt\\d+)$", names_to = c("colonne", "ttt")) %>%
      tidyr::pivot_wider(names_from = colonne, values_from = value)

    # Adding stopping boundaries
    resp_list <- resp_list %>%
      dplyr::left_join(stopbound %>% dplyr::select(-seuil_notox), by = c("nb_ana", "tot_eff_ttt0" = "y_cont")) %>%
      dplyr::left_join(stopbound %>% dplyr::select(-seuil_eff), by = c("nb_ana", "tot_notox_ttt0" = "y_cont"))

    # Decision rules
    resp_list <- resp_list %>%
      dplyr::mutate(decision = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Continue",
                                                (nb_ana == max(nb_ana)) & tot_eff > seuil_eff & tot_notox > seuil_notox ~ "Accept the treatment",
                                                (nb_ana != max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Early stopping",
                                                (nb_ana == max(nb_ana)) & (tot_eff <= seuil_eff | tot_notox <= seuil_notox) ~ "Stopping"),
                    decision_eff = dplyr::case_when((nb_ana != max(nb_ana)) & tot_eff > seuil_eff ~ "Continue",
                                                    (nb_ana == max(nb_ana)) & tot_eff > seuil_eff ~ "Accept the treatment",
                                                    (nb_ana != max(nb_ana)) & tot_eff <= seuil_eff ~ "Early stopping (futility)",
                                                    (nb_ana == max(nb_ana)) & tot_eff <= seuil_eff ~ "Stopping (futility)"),
                    decision_tox = dplyr::case_when((nb_ana != max(nb_ana)) & tot_notox > seuil_notox ~ "Continue",
                                                    (nb_ana == max(nb_ana)) & tot_notox > seuil_notox ~ "Accept the treatment",
                                                    (nb_ana != max(nb_ana)) & tot_notox <= seuil_notox ~ "Early stopping (toxicity)",
                                                    (nb_ana == max(nb_ana)) & tot_notox <= seuil_notox ~ "Stopping (toxicity)")) %>%
      dplyr::filter(decision != "Continue") %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

  }

  # Output of the function
  if (n_bras == 1) {
    return(list(
      caracteristique = data.frame(
        alpha = alpha,
        p_n = paste0("(", paste(p_n, collapse = "/"), ")"),
        p_a = paste0("(", paste(p_a, collapse = "/"), ")"),
        p_reel = paste0("(", paste(unlist(p_reel), collapse = "/"), ")"),
        methode = methode,
        rejet_h0 = mean(resp_list$decision == "Accept the treatment"),
        arret_precoce = mean(resp_list$decision == "Early stopping"),
        arret_precoce_fut = mean(resp_list$decision_eff == "Early stopping (futility)"),
        arret_precoce_tox = mean(resp_list$decision_tox == "Early stopping (toxicity)"),
        arret_fut = mean(stringr::str_detect(resp_list$decision_eff, "futility")),
        arret_tox = mean(stringr::str_detect(resp_list$decision_tox, "toxicity")),
        nb_pts = mean(resp_list$tot_pat),
        nb_eff = mean(resp_list$tot_eff),
        nb_tox = mean(resp_list$tot_pat - resp_list$tot_notox),
        gamma_eff = gamma_eff,
        gamma_tox = gamma_tox
      ),
      essais = resp_list
    ))
  } else {
    return(list(
      carac_globales = data.frame(
        alpha = alpha,
        p_n = paste0("(", paste(p_n, collapse = "/"), ")"),
        p_a = paste0("(", paste(p_a, collapse = "/"), ")"),
        p_reel = paste(lapply(seq_along(p_reel), FUN = function(x) paste0(x, ":(", paste(p_reel[[x]], collapse = "/"), ")")), collapse = ", "),
        methode = methode,
        gamma_eff = gamma_eff,
        gamma_tox = gamma_tox,
        rejet_glob = resp_list %>%
          dplyr::group_by(n_simu) %>%
          dplyr::summarise(rejet_glob = sum(decision == "Accept the treatment") > 0, .groups = "drop") %>%
          dplyr::summarise(rejet_glob = mean(rejet_glob)) %>%
          dplyr::pull(rejet_glob)
      ),
      carac_bras = data.frame(
        summarise_decision(resp_list, ttt, decision, "Accept the treatment", rejet_h0) %>%
          dplyr::left_join(summarise_decision(resp_list, ttt, decision, "Early stopping", arret_precoce), by = "ttt") %>%
          dplyr::left_join(summarise_detect(resp_list, ttt, decision, "Stopping|stopping", arret), by = "ttt") %>%
          dplyr::left_join(summarise_decision(resp_list, ttt, decision_eff, "Early stopping (futility)", arret_precoce_fut), by = "ttt") %>%
          dplyr::left_join(summarise_decision(resp_list, ttt, decision_tox, "Early stopping (toxicity)", arret_precoce_tox), by = "ttt") %>%
          dplyr::left_join(summarise_detect(resp_list, ttt, decision_eff, "futility", arret_fut), by = "ttt") %>%
          dplyr::left_join(summarise_detect(resp_list, ttt, decision_tox, "toxicity", arret_tox), by = "ttt") %>%
          dplyr::left_join(summarise_ttt(resp_list, ttt, tot_pat), by = "ttt") %>%
          dplyr::left_join(summarise_ttt(resp_list, ttt, tot_eff), by = "ttt") %>%
          dplyr::left_join(summarise_ttt(resp_list, ttt, tot_pat - tot_notox, tot_tox), by = "ttt")
      ),
      essais = resp_list
    ))
  }

}
