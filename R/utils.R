"%nin%" <- function(x, y) match(x, y, nomatch = 0L) == 0L

#' @importFrom rlang {{
#' @importFrom magrittr %>%
summarise_ttt <- function(data, groupe, colonne, intitule = {{colonne}}) {
  data %>%
    dplyr::group_by({{groupe}}) %>%
    dplyr::summarise("{{intitule}}" := mean({{colonne}}))
}

#' @importFrom rlang {{
#' @importFrom magrittr %>%
summarise_decision <- function(data, groupe, colonne, char_decision, intitule) {
  data %>%
    dplyr::group_by({{groupe}}) %>%
    dplyr::summarise("{{intitule}}" := mean({{colonne}} == char_decision))
}

#' @importFrom rlang {{
#' @importFrom magrittr %>%
summarise_detect <- function(data, groupe, col_decision, char_decision, intitule) {
  data %>%
    dplyr::group_by({{groupe}}) %>%
    dplyr::summarise("{{intitule}}" := mean(stringr::str_detect({{col_decision}}, char_decision)))
}
