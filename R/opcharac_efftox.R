#' Operating characteristics
#'
#' Obtain the operating characteristics for BOP2 vs a reference value or a control group.
#'
#' Optimize Cn and then simulate \code{nsim_essais} with law \code{p_reel} to determine the operating characteristics. It is possible to optimize
#' Cn before with \code{deter_cutoff} and then supply these values as argument in \code{power_seq} and \code{cut_seq}. It allows reduction of the execution
#' if we want to simulate several scenarios with the same Cn.
#'
#' Works for mono- and multiarm trials.
#'
#' For multiarm settings, \code{bonf = TRUE} to perform Bonferroni correction (less time, but more conservative or even Cn).
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
#' @param alpha Maximal type I error rate (FWER in case of multiarm).
#' @param bonf TRUE to apply a Bonferroni correction in case of multiarm settings.
#' @param nsim_oc,nsim_essai Number of simulated trials to optimize CN and get the operating characteristics.
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value). Not yet implemented.
#' @param p_reel,p_n,prior,p_a Vectors of length 4 giving the law of probability under the 'real' scenario, the inefficacy/toxicity hypothesis, the prior
#' (by default \code{p_n}) and the efficacy/non toxicity hypothesis.
#' @param phi Thresholds for efficacy and toxicity (by default the values under inefficacy/toxicity hypothesis). If specified, the values for efficacy and non toxicity.
#' @param delta Minimal desirable difference of effect if compared with a control arm. Should be between -1 and 1.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param cut_seq,power_seq Vectors of values used to choose the couple while optimizing the Cn. If a unique value is specified for each parameter, the grid search won't be performed.
#' @param seed_cut,seed_pat Seeds for \code{rmultinom} (for Cn and operating characteristics).
#' @param methode Method used to generate the data:
#' \itemize{
#'   \item 1: trial by trial;
#'   \item 2: trial by trial, patient by patient;
#'   \item 3: whole in 1.
#' }
#' @param tableau_essais Data.frame of simulated trials with \code{gen_patients_multinom}. If NULL, the function will simulated them using \code{p_reel}.
#'
#' @examples
#' # For BOP2 simulated trials in Zhou, 2017
#' opcharac_efftox(ana_inter = c(10, rep(5, 6)),
#'                 alpha = 0.1,
#'                 p_reel = c(0.18, 0.42, 0.02, 0.38),
#'                 p_n = c(0.15, 0.3, 0.15, 0.4), p_a = c(0.18, 0.42, 0.02, 0.38),
#'                 nsim_oc = 10000, nsim_essai = 10000,
#'                 seed_cut = 1024, seed_pat = 1024,
#'                 methode = 3L)
#'
#' # For a 3 arm-trial
#' opcharac_efftox(ana_inter = rep(15, 4),
#'                 alpha = 0.1,
#'                 p_reel = list(ttt1 = c(0.15, 0.5, 0.05, 0.3), ttt2 = c(0.15, 0.3, 0.25, 0.3), ttt3 = c(0.05, 0.15, 0.15, 0.65)),
#'                 p_n = c(0.1, 0.2, 0.2, 0.5), p_a = c(0.15, 0.5, 0.05, 0.3),
#'                 nsim_oc = 10000, nsim_essai = 10000)
#'
#' @export
opcharac_efftox <- function(alpha = .1,
                            bonf = FALSE,
                            nsim_oc = 10000,
                            nsim_essai = 10000,
                            ana_inter,
                            ana_inter_tox = NULL,
                            rand_ratio = 1,
                            p_n,
                            prior = NULL,
                            p_a,
                            p_reel,
                            phi = NULL,
                            delta = NULL,
                            mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                            cut_seq = seq(0.5, 0.95, by = 0.005),
                            power_seq = seq(0, 1, by = 0.01),
                            seed_cut = 1024,
                            seed_pat = 1993,
                            methode = 2L,
                            tableau_essais = NULL) {
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
  if (is.null(prior)) {
    prior <- p_n
    message("\"p_n\" taken as prior.")
  } else {
    if (length(prior) != 4)
      stop("Vector \"prior\" should be of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
    if (!dplyr::near(sum(prior), 1))
      warning("Vector \"prior\" should sum to 1.")
  }
  if (!dplyr::near(sum(p_a), 1))
    stop("Vector \"p_a\" should sum to 1.", call. = FALSE)
  if (length(p_a) != 4)
    stop("\"p_a\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (!is.list(p_reel))
    p_reel <- list(ttt1 = p_reel)
  if (any(lapply(p_reel, length) != 4))
    stop("Vectors \"p_reel\" shoud be of length 4 for each treatment:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (any(unlist(lapply(lapply(p_reel, sum),
                        FUN = function(x) !dplyr::near(x, 1)))))
    stop("Vectors \"p_reel\" should sum to 1 for each treatment.", call. = FALSE)
  if (length(nsim_oc) != 1 || !is.numeric(nsim_oc) || nsim_oc < 0 || nsim_oc %% 1 != 0)
    stop("\"nsim_oc\" should be a positive integer.", call. = FALSE)
  if (length(nsim_essai) != 1 || !is.numeric(nsim_essai) || nsim_essai < 0 || nsim_essai %% 1 != 0)
    stop("\"nsim_essai\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.", call. = FALSE)
  }

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
  if (!is.null(delta)) {
    if (any(delta > 1) | any(delta < -1))
      stop("\"delta\" should be a real between -1 and 1.", call. = FALSE)
  }
  if (lapply(lapply(rand_ratio,
                    function(x) ana_inter * x),
             function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
    warning("\"rand_ratio\" x \"ana_inter\" must only return integers.
         The numbers were converted to integers using the ceiling function.",
            call. = FALSE, immediate. = TRUE)
    anas_inters <- lapply(anas_inters,  ceiling)
  }

  # Optimization of Cn if needed
  if (length(cut_seq) == 1 & length(power_seq) == 1) {
    C_ <- cut_seq
    gamm <- power_seq
  } else {
    carac_optimales <- suppressMessages(suppressWarnings(
      deter_cutoff(
        alpha = alpha,
        n_bras = n_bras,
        bonf = bonf,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        mat_beta_xi = mat_beta_xi,
        phi = phinotox,
        delta = delta,
        p_n = p_n,
        prior = prior,
        p_a = p_a,
        nsim_oc = nsim_oc,
        cut_seq = cut_seq,
        power_seq = power_seq,
        seed = seed_cut,
        methode = methode,
        affich_mat = "No"
      )
    ))

    C_ <- carac_optimales[[1]]["C_"] # Store the 2 hyperparameters of Cn
    gamm <- carac_optimales[[1]]["gamma"]
  }

  # resp_list give the number of responses and toxicities at each step of each trial
  if (is.null(tableau_essais)) {
    resp_list <- suppressMessages(suppressWarnings(
      gen_patients_multinom(
        n_sim = nsim_essai,
        ana_inter = ana_inter,
        ana_inter_tox = ana_inter_tox,
        rand_ratio = rand_ratio,
        multinom_ttt = p_reel,
        multinom_cont = if (is.null(delta)) {NULL} else {p_n},
        seed = seed_pat,
        methode = methode
      )
    ))
  } else {
    resp_list <- tableau_essais
  }

  phi_eff <- phitox[1] # Phi for efficacy
  phi_tox <- phitox[2] # Phi for toxicity
  prior_eff <- sum(prior * mat_beta_xi[1,])
  prior_tox <- sum(prior * (1 - mat_beta_xi[2,]))

  # Stopping boundaries
  stopbound <- suppressMessages(suppressWarnings(
    get_stopbound(
      ana_inter = ana_inter,
      ana_inter_tox = ana_inter_tox,
      rand_ratio = rand_ratio,
      C_ = C_,
      gamm = gamm,
      p_n = p_n,
      prior = prior,
      phi = phinotox,
      delta = delta
    )
  ))

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

  } else { # Vs control group

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
      dplyr::mutate(
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
        C = C_,
        gamma = gamm
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
        C_ = C_,
        gamma = gamm,
        rejet_glob = resp_list %>%
          dplyr::group_by(n_simu) %>%
          dplyr::summarise(
            rejet_glob = sum(decision == "Accept the treatment") > 0,
            .groups = "drop"
          ) %>%
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
