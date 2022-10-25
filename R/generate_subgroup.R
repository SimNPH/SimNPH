# TODO: refactor/reparametrize to subgroup, prevalence and hr in subgroup and total pop

#' Generate Dataset with different treatment effect in subgroup
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condidtion has to contain the following columns:
#'
#'   * n_trt number of paitents in treatment arm
#'   * n_ctrl number of patients in control arm
#'   * hazard_ctrl hazard in the control arm
#'   * hazard_trt hazard in the treatment arm for not cured patients
#'   * hazard_subgroup hazard in for cured patients
#'   * prevalence proportion of cured patients
#'
#' @return
#' For generate_subgroup: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_subgroup simulates a dataset with a mixture of cured
#'   patients
#'
#' @examples
#' one_simulation <- merge(
#'     assumptions_subgroup(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1) |>
#'   generate_cure_mixture()
#' head(one_simulation)
#' tail(one_simulation)
generate_subgroup <- function(condition, fixed_objects=NULL){
  # if t_max is not given in fixed_objects
  if(is.null(fixed_objects) || (!hasName(fixed_objects, "t_max"))){
    # set t_max to 1-1/10000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(10000) / condition$hazard_ctrl,
      log(10000) / condition$hazard_trt
    )
  } else {
    t_max <- fixed_objects$t_max
  }

  if (condition$prevalence < 0 || condition$prevalence > 1) {
    stop(gettext("Subgroup prevalence has to be between 0 and 1"))
  }

  data_trt <- data.frame(
    t = nph::rSurv_fun(
      condition$n_trt,
      nph::pop_pchaz(
        Tint = c(0, t_max),
        lambdaMat1 = matrix(c(
          condition$hazard_trt,
          condition$hazard_subgroup
        ), ncol=1),
        lambdaMat2 = matrix(c(0, 0), ncol=1),
        lambdaProgMat = matrix(c(0, 0), ncol=1),
        p=c(1-condition$prevalence, condition$prevalence)
      )
    ),
    trt = 1,
    evt = TRUE
  )

  data_ctrl <- data.frame(
    t = nph::rSurv_fun(
      condition$n_ctrl,
      nph::pchaz(
        c(0, t_max),
        c(condition$hazard_ctrl)
      )
    ),
    trt = 0,
    evt = TRUE
  )

  rbind(data_trt, data_ctrl)
}

#' Create an empty assumtions data.frame for generate_cure_mixture
#'
#' @return For assumptions_subgroup: a design tibble with default values invisibly
#'
#' @details assumptions_subgroup prints the code to generate a default
#'   design tibble for use with generate_cure_mixture and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_subgroup generate default assumptions tibble
#'
#' @examples
#' Design <- assumptions_subgroup()
#' Design
assumptions_subgroup <- function(){
  skel <- "expand.grid(
  hazard_ctrl       = 0.003795467,       # hazard under control (med. survi. 6m)
  hazard_trt        = 0.001265156,       # hazard under treatment (med. surv. 18m)
  hazard_subgroup   = 9.488668e-05,      # hazard cured (med. surv. 20y)
  prevalence        = seq(0, 1, by=0.2), # proportion of patients belonging to subgroup
  random_withdrawal = 0.01               # rate of random withdrawal
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}
