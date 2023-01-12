#' Generate data
#'
#' Generate the data of \code{n_sim} trials.
#'
#' Uses the multinomial law to simulate \code{n_sim} trials with \code{rmultinom}.
#'
#' With or without control group, and with 1 final analysis or possible multiple interim analyses.
#'
#' Data are in long format.
#'
#' By default we chose method 2 because it seemed to be less biased (but all method are very close).
#'
#' @return A data.frame with 1 row by interim analysis and treatment group, and 10 variables:
#' \itemize{
#'   \item n_sim: the index of the trial;
#'   \item nb_ana: the index of the interim analysis;
#'   \item ttt: the treatment arm;
#'   \item efftox, effnotox, noefftox, noeffnotox: cumulative number of each outcome at each interim analysis for each trial;
#'   \item tot_pat, tot_eff, tot_notox: the total of patients, responses and non toxicities at each interim analysis for each trial.
#' }
#'
#' @encoding UTF-8
#'
#' @param n_sim The desired number of simulated trials.
#' @param ana_inter Vector giving the number of additional patients at each interim analysis. If ana_inter_tox is supplied, represents the vector for efficacy interim analyses.
#' @param ana_inter_tox Vector of the number of additional patients at each interim analysis for toxicity. If analyses for effficacy and toxicity occur at the same
#' number of patients, set it to NULL.
#' @param rand_ratio Number or numeric vector representing the ratio between the number of patients in treatment arms and control arm (or some reference sample size in case
#' of comparison to a reference value).
#' @param multinom_cont Vector of length 4 describing the multinomial law followed by the control group ("Eff/Tox", "Eff/NoTox", "NoEff/Tox", "NoEff/NoTox").
#' If no control group, let it set to default NULL.
#' @param multinom_ttt A list of vectors of length 4 describing the multinomial law followed by each treatment group ("Eff/Tox", "Eff/NoTox", "NoEff/Tox", "NoEff/NoTox").
#' If only one treatment group, one can just specify a vector of length 4.
#' @param mat_beta_xi The link matrix between outcomes of the multinomial distribution and endpoints. (Second row is for non toxicity instead of toxicity).
#' @param seed The seed for \code{rmultinom}.
#' @param methode Method used to generate the data:
#' \itemize{
#'   \item 1: trial by trial;
#'   \item 2: trial by trial, patient by patient;
#'   \item 3: whole in 1.
#' }
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' # 100 trials of 40 patients with an analysis after each enrollment of 10 patients. Only 1 treatment group.
#' gen_patients_multinom(n_sim = 100,
#'                       ana_inter = rep(10, 4),
#'                       multinom_ttt = rep(0.25, 4))
#'
#' # 100 trials of 40 patients with an analysis after each enrollment of 10 patients. 3 treatment groups.
#' gen_patients_multinom(n_sim = 100,
#'                       ana_inter = rep(10, 4),
#'                       multinom_ttt = list(rep(0.25, 4), c(0.2, 0.4, 0.2, 0.2), c(0.1, 0.4, 0.05, 0.45)))
#'
#' # 100 trials of 40 patients with an analysis after each enrollment of 10 patients. 1 treamtent group with 1 control group.
#' gen_patients_multinom(n_sim = 100,
#'                       ana_inter = rep(10, 4),
#'                       multinom_ttt = rep(0.25, 4),
#'                       multinom_cont = c(0.2, 0.3, 0.1, 0.4))
#'
#' # 100 trials of 40 patients with an analysis after each enrollment of 10 patients. 3 treatment groups with a control group.
#' gen_patients_multinom(n_sim = 100,
#'                       ana_inter = rep(10, 4),
#'                       multinom_ttt = list(rep(0.25, 4), c(0.2, 0.4, 0.2, 0.2), c(0.1, 0.4, 0.05, 0.45)),
#'                       multinom_cont = c(0.2, 0.3, 0.1, 0.4))
#'
#' # 100 trials of 40 patients with an analysis after each enrollment of 10 patients. 3 treatment groups with a
#' # 40 patients control group. Each group has a different number of patients
#' # (for the 0.75 ratio, 0.75 * 10 = 7.5 rounded to 8 additional patients at each interim analysis).
#' gen_patients_multinom(n_sim = 100,
#'                       rand_ratio = c(1, .75, 1.5),
#'                       ana_inter = rep(10, 4),
#'                       multinom_ttt = list(rep(0.25, 4), c(0.2, 0.4, 0.2, 0.2), c(0.1, 0.4, 0.05, 0.45)),
#'                       multinom_cont = c(0.2, 0.3, 0.1, 0.4))
#' @export
gen_patients_multinom <- function(n_sim,
                                  ana_inter,
                                  ana_inter_tox = NULL,
                                  rand_ratio = NULL,
                                  multinom_cont = NULL,
                                  multinom_ttt = list(),
                                  mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                  seed = 1024,
                                  methode = 2L) {
  # Checks arguments
  if (length(n_sim) != 1 || !is.numeric(n_sim) || n_sim < 0 || n_sim %% 1 != 0)
    stop("\"n_sim\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (methode %nin% c(1L, 2L, 3L))
    stop(
      "Choice between 3 methods with an integer:
              * 1L = trial by trial;
              * 2L = trial by trial, patient by patient;
              * 3L = whole in 1.
      Speed is following: 3L > 1L > 2L.",
      call. = FALSE
    )

  # Checking the specified probabilities
  if (length(multinom_ttt) == 0)
    stop("No specified treatment.", call. = FALSE)
  if (!is.list(multinom_ttt)) {
    multinom_ttt <- list(ttt1 = multinom_ttt)
    warning("\"multinom_ttt\" should be a list. Converted with list(multinom_ttt) assuming that there is only 1 treatment group.", call. = FALSE, immediate. = TRUE)
  }
  if (any(lapply(multinom_ttt, length) != 4))
    stop("Each law in \"multinom_ttt\" should be a vector of length 4.", call. = FALSE)
  if (any(unlist(lapply(multinom_ttt, FUN = function(x) !dplyr::near(sum(x), 1)))))
    stop("For each treatment the law shoudl sum to 1.", call. = FALSE)
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
      stop("A vector of length 4 must be provided if you want to specify the control group.", call. = FALSE)
    if (!dplyr::near(sum(multinom_cont), 1))
      stop("Control law should sum to 1.", call. = FALSE)
    if (is.null(names(multinom_ttt))) {
      noms_ttt <- c("ttt0", paste0("ttt", seq_len(length(multinom_ttt))))
    } else {
      noms_ttt <- c("ttt0", names(multinom_ttt))
    }
    proba <- append(list(multinom_cont), multinom_ttt)
    names(proba) <- noms_ttt
  }
  # ttt0 is always for control group, and ttt1, 2, ... will be evaluated treatments.
  # If no control group, start at ttt1.

  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.")
    anas_inters_cum <- sort(union(cumsum(ana_inter), cumsum(ana_inter_tox)))
    ana_inter     <- c(anas_inters_cum[1], diff(anas_inters_cum))
  }

  # Get the sample sizes at each interim analysis
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

  anas_inters <- lapply(
    rand_ratio,
    FUN = function(x, y) {x * y},
    y = ana_inter
  ) # Number of additional patients at each interim analysis for each group
  if (lapply(anas_inters, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
    warning("\"rand_ratio\" x \"ana_inter\" must only return integers.
         The numbers were converted to integers using the ceiling function.",
         call. = FALSE, immediate. = TRUE)
    anas_inters <- lapply(anas_inters,  ceiling)
  }
  names(anas_inters) <- noms_ttt
  nmax <- lapply(anas_inters, sum)

  # Number of subjects at each interim analysis for each arm (including the control arm)
  sujet_ana <- lapply(anas_inters, function(x) {
    data.frame(nb_ana = as.integer(seq_along(x)),
               nb_patients = x)
  })
  set.seed(seed)
  on.exit(set.seed(NULL), add = TRUE) # Reset the seed to not impair the environment with the use of the function

  if (methode == 3L) {

    liste_patients <- list()
    for (i in noms_ttt) {
      resp_list <- lapply(
        seq_along(anas_inters[[i]]),
        FUN = function(n) {
          n_ana <- anas_inters[[i]][n]
          mat <- rmultinom(n_sim, size = n_ana, prob = proba[[i]])
          mat <- apply(
            X = mat,
            MARGIN = 2,
            FUN = function(p) {
              efftox <- p[1]
              effnotox <- p[2]
              noefftox <- p[3]
              noeffnotox <- p[4]
              eff <- sum(mat_beta_xi[1,] * p)
              notox <- sum(mat_beta_xi[2,] * p)
              return(c(efftox, effnotox, noefftox, noeffnotox, eff, notox))
            }
          ) %>%
            t() %>%
            as.data.frame()
          colnames(mat) <- paste0(c("efftox", "effnotox", "noefftox", "noeffnotox", "eff", "notox"), "_ana", n)
          mat$n_simu <- seq_len(n_sim)
          return(mat)
        }
      )
      resp_list <- purrr::reduce(.x = resp_list,
                                 .f = dplyr::left_join,
                                 by = "n_simu")

      # Adding total counts at each interim analysis
      resp_list <- resp_list %>%
        tidyr::pivot_longer(
          cols = -n_simu,
          names_pattern = "^(efftox|effnotox|noefftox|noeffnotox|eff|notox)_ana(\\d+)$",
          names_to = c("critere", "nb_ana")
        ) %>%
        tidyr::pivot_wider(names_from = critere, values_from = value) %>%
        dplyr::mutate(nb_ana = as.integer(nb_ana)) %>%
        dplyr::left_join(sujet_ana[[i]], by = "nb_ana") %>%
        dplyr::group_by(n_simu) %>%
        dplyr::mutate(
          efftox = cumsum(efftox),
          effnotox = cumsum(effnotox),
          noefftox = cumsum(noefftox),
          noeffnotox = cumsum(noeffnotox),
          tot_eff = cumsum(eff),
          tot_notox = cumsum(notox),
          tot_pat = cumsum(nb_patients)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-eff,-notox,-nb_patients)
      liste_patients[[i]] <- resp_list %>%
        dplyr::mutate(ttt = i, .before = 3) %>%
        dplyr::relocate(tot_pat, .before = 8)
    }

    liste_patients <- do.call("rbind", liste_patients) %>%
      dplyr::arrange(n_simu, ttt, nb_ana)

  } else if (methode == 1) {

    tab_patients <- lapply(names(sujet_ana), function(x) cbind(sujet_ana[[x]], ttt = x)) %>%
      do.call(what = "rbind", arg = .)

    # Sample size for each treatment at each interim analysis
    liste_patients <- expand.grid(
      n_simu = seq_len(n_sim),
      nb_ana = seq_len(length(ana_inter)),
      ttt = noms_ttt
    ) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::left_join(tab_patients, by = c("nb_ana", "ttt")) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(tot_pat = cumsum(nb_patients)) %>%
      dplyr::ungroup()

    # Efficacy and toxicity probabilities
    liste_patients <- liste_patients %>%
      dplyr::mutate(
        p_mult = purrr::map(ttt, ~ proba[[match(.x, names(proba))]]),
        n_mult = purrr::map2(p_mult, nb_patients, function(x, y) {
          mat <- as.double(rmultinom(n = 1, size = y, prob = x))
          names(mat) <- c("efftox", "effnotox", "noefftox", "noeffnotox")
          return(mat)
        })
      ) %>%
      tidyr::unnest_wider(n_mult) %>%
      dplyr::mutate(tot_eff = efftox + effnotox,
                    tot_notox = effnotox + noeffnotox) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(dplyr::across(efftox:tot_notox, ~ cumsum(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-nb_patients, -p_mult) %>%
      dplyr::relocate(tot_pat, .before = 9)

  } else if (methode == 2) {

    anas_inters_cum <- lapply(anas_inters, cumsum)
    liste_patients <- list()

    for (p in noms_ttt) {

      liste_patients[[p]] <- do.call("rbind",
                                     lapply(
                                       seq_len(n_sim),
                                       function(x) {
                                         tableau <- rmultinom(nmax[[p]], size = 1, prob = proba[[p]]) %>%
                                           t() %>%
                                           as.data.frame()
                                         tableau <- cbind(n_simu = x, ttt = p, tableau, tot_pat = seq_len(nmax[[p]]))
                                         tableau$nb_ana <- NA
                                         for (i in rev(anas_inters_cum[[p]])) {tableau$nb_ana <- ifelse(tableau$tot_pat <= i, match(i, anas_inters_cum[[p]]), tableau$nb_ana)}
                                         names(tableau)[3:6] <- c("efftox", "effnotox", "noefftox", "noeffnotox")
                                         return(tableau)
                                       }
                                     ))
    }

    liste_patients <- do.call("rbind", liste_patients) %>%
      dplyr::group_by(n_simu, ttt, nb_ana) %>%
      dplyr::summarise(dplyr::across(c(efftox:noeffnotox), ~ sum(.x)),
                       tot_pat = dplyr::last(tot_pat),
                       .groups = "drop") %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(dplyr::across(c(efftox:noeffnotox), ~ cumsum(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(tot_eff = efftox + effnotox,
                    tot_notox = effnotox + noeffnotox)

  }

  return(liste_patients)

}





#' Resample the interim analyses
#'
#' Resample the interim analyses based on the table with interim analyses after each patients.
#'
#' @param tab_patients Dataset generated with \code{gen_patients_multinom} with ana_inter = rep(1, X).
#' @param nouv_ana_inter The new vector of additional patients at each interim analysis. If several treatments,
#' you can either give 1 vector that will be applied to each group or give a nammed list with names corresponding to treatments.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' tab1 <- gen_patients_multinom(n_sim = 100,
#'                               ana_inter = rep(1, 40),
#'                               multinom_ttt = list(rep(0.25, 4), c(0.2, 0.4, 0.2, 0.4), c(0.1, 0.4, 0.05, 0.45))
#' tab2 <- modif_tab_patients(tab_patients = tab1,
#'                            nouv_ana_inter = rep(10, 4))
#' @export
modif_tab_patients <- function(tab_patients,
                               nouv_ana_inter) {

  if (any(!names(tab_patients) %in% c("n_simu", "nb_ana", "ttt", "efftox", "effnotox", "noefftox", "noeffnotox", "tot_pat", "tot_eff", "tot_notox")))
    stop("\"tab_patients\" should be created with function \"gen_patients_multinom\".", call. = FALSE)
  if (any(unique(tab_patients$tot_pat) != (1:max(tab_patients$tot_pat))))
    stop("\"tab_patients\" should be a table containing patient's analyses 1 by 1.", call. = FALSE)
  if (!is.list(nouv_ana_inter)) nouv_ana_inter <- list(nouv_ana_inter)
  if (length(nouv_ana_inter) == 1 & length(unique(tab_patients$ttt)) != 1) {
    message("\"nouv_ana_inter\" will be applied for all groups.")
  } else if (length(nouv_ana_inter) != length(unique(tab_patients$ttt))) {
    stop("If it's not the same for all arms, you should specify one set of interim analyses by arm.", call. = FALSE)
  }

  # Get the cumulative number of patients
  nouv_ana_inter <- lapply(nouv_ana_inter, cumsum)

  if (length(nouv_ana_inter) == 1) {
    if (max(nouv_ana_inter[[1]]) > max(tab_patients$tot_pat))
      warning("Some interim analyses won't be performed because \"tab_patients\" doesn't contain enough patients.", call. = FALSE, immediate. = TRUE)
    tab_patients <- tab_patients %>%
      dplyr::left_join(data.frame(n_ana = seq_along(nouv_ana_inter[[1]]), tot_pat = nouv_ana_inter[[1]]), by = "tot_pat") %>%
      dplyr::filter(!is.na(n_ana)) %>%
      dplyr::mutate(nb_ana = n_ana) %>%
      dplyr::select(-n_ana)
  } else {
    tab_ana <- lapply(seq_along(nouv_ana_inter),
           FUN = function(x) {
             data.frame(ttt = names(nouv_ana_inter)[x],
                        n_ana = seq_along(nouv_ana_inter[[x]]),
                        tot_pat = nouv_ana_inter[[x]])
           }) %>%
      do.call(what = "rbind", args = .)

    test <- tab_patients %>%
      dplyr::group_by(ttt) %>%
      dplyr::summarise(maximum = max(tot_pat)) %>%
      dplyr::right_join(tab_ana, by = "ttt") %>%
      dplyr::group_by(ttt) %>%
      dplyr::filter(tot_pat == max(tot_pat)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(test = maximum <= tot_pat) %>%
      dplyr::pull(test)
    if (sum(is.na(test)) | !all(test, na.rm = TRUE))
      warning("Some interim analyses won't be performed because \"tab_patients\" doesn't contain enough patients.", call. = FALSE, immediate. = TRUE)

    tab_patients <- tab_patients %>%
      dplyr::left_join(tab_ana, by = c("tot_pat", "ttt")) %>%
      dplyr::filter(!is.na(n_ana)) %>%
      dplyr::mutate(nb_ana = n_ana) %>%
      dplyr::select(-n_ana)
  }

  return(tab_patients)

}
