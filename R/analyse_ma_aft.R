#' Analyse Dataset with model averaging techniques
#'
#' @param level confidence level for CI computation
#' @param dist Which models should be fit to the data?
#' @param alternative alternative hypothesis for the tests "two.sided"
#'   or "one.sided"
#'   @param n_boot Number of bootstrapped datasets
#'
#' @return an analysis function that returns a data.frame
#'
#' @export
#'
#' @details the columns in the return are the two-sided p-value for the test of
#'   equal medians. The estimated medians in the treatment and control group and
#'   the estimated difference in median survival with confidence intervals.
#'
#'   The estimates and tests are comstructed by fitting seperate Weibull
#'   regression models in the treatment and control groups and then estimating
#'   the medians and respective variances with the delta-method.
#'
#' @author Andrew Hooker
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(3) |>
#'   tail(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_ma()(condition, dat)

analyse_ma_aft <- function(level=0.95, alternative = "two.sided",
                       dist=c("gengamma","genf","weibull",
                              "gamma","llogis",
                              "lnorm"),
                       n_boot=1000) {

  # Checks -----
  stopifnot(alternative %in% c("two.sided", "one.sided"))
  stopifnot(all(dist %in% c("gengamma","genf","weibull",
                            "gamma","llogis",
                            "lnorm")))

  capture_eval <- function(expr, env = parent.frame()){
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(eval(expr, envir = env), error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(value = res, warning = warn, error = err)
  }

  # wrap in SimDesign function
  function(condition, dat, fixed_objects = NULL){


    dat <- dat |> dplyr::mutate(trt=factor(trt))


    # bootstrap the dataset ----
    ### control group
    TRT0 <- dat |> dplyr::filter(trt==0)
    bootstraped_survival_datasets_TRT0<-list()
    for (i in 1:n_boot){
      TRTn<- TRT0[sample(x = nrow(TRT0), size = nrow(TRT0), replace = T), ]
      bootstraped_survival_datasets_TRT0[[i]]<-TRTn
    }

    ###active group
    TRT1 <- dat |> dplyr::filter(trt==1)
    bootstraped_survival_datasets_TRT1<-list()
    for (i in 1:n_boot){
      TRTn<- TRT1[sample(x = nrow(TRT1), size = nrow(TRT1), replace = T), ]
      bootstraped_survival_datasets_TRT1[[i]]<-TRTn
    }

    ### merge the lists containing bootstraped datasets from arm TRT0 and TRT1
    FUN_merge_lists <- function(X,Y) {
      merged_df<-rbind(X, Y)
      return(merged_df)
    }

    TRT_merged <- mapply(bootstraped_survival_datasets_TRT0,
                         bootstraped_survival_datasets_TRT1,
                         FUN = FUN_merge_lists, SIMPLIFY = FALSE)

    FUN_fit_models <- function(data){
      results <- c()
      for(dist_i in dist){
        # dist_i <- "weibull"

        fit_tmp <- capture_eval(
          flexsurv::flexsurvreg(Surv(t, evt)~trt,
                                data=data,dist=dist_i)
          )$value
        if(is.null(fit_tmp)) next


        # fit_tmp <-flexsurv::flexsurvreg(Surv(t, evt)~trt,
        #                                 data=data,dist=dist_i)
        aic_tmp <- AIC(fit_tmp)
        cov_par_tmp <- coefficients(fit_tmp)[["trt1"]]

        # the time acceleration factor.
        # if TAF > 1  indicate that higher covariate values are associated
        # with a higher risk of the event, or shorter times to the event.
        # 1/TAF is the ratio of expected times to the event between covariate
        # values of 1 and 0.
        if(dist_i %in% c("weibull","llogis","lnorm",
                         "gengamma","genf"))
          taf_tmp <- exp(-cov_par_tmp)
        if(dist_i %in% c("gamma","exp","gompertz"))
          taf_tmp <- exp(cov_par_tmp)

        res <- tibble::tibble(fit=list(fit_tmp),
                              aic=aic_tmp,
                              taf=taf_tmp,
                              dist=dist_i)
        results <- dplyr::bind_rows(results,res)
      }

      results_best <- results |> dplyr::slice(which.min(aic))
      return(results_best)
    }

    # FUN_fit_models(dat)

    fit_res <-
      lapply(TRT_merged,
             FUN = FUN_fit_models)

    fit_res_combined <- purrr::map_dfr(fit_res, dplyr::bind_rows)

    # Look at one result
    # fit_tmp <- results[3,"fit"][[1]][[1]]
    # fit_tmp$res
    # coef(fit_tmp)
    # plot(fit_tmp,col = c("red","blue"),
    #      col.obs = c("red","blue"),lty.obs = "dotted")


    # compute median and percentiles
    res_med_ci <- fit_res_combined |>
      dplyr::reframe(qs=quantile(
        taf,
        c((1-level)/2, 0.5, level+(1-level)/2)),
        prob = c((1-level)/2, 0.5, level+(1-level)/2))
    # res_med_ci

    # 2-sided Statistcal significance of the  TAF
    # res_sig_diff <- res_med_ci |>
    #   filter(prob!=0.5) |>
    #   dplyr::mutate(across(everything(),~.x-1)) |>
    #   dplyr::mutate(across(everything(),sign)) |>
    #   dplyr::reframe(sig_diff=sum(qs)) |>
    #   dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))
    # res_bestfit_sig_diff

    p_val_greater <- fit_res_combined |>
      dplyr::mutate(val_gt_null=dplyr::if_else(taf<=1,1,0)) |>
      dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

    p_val_less <- fit_res_combined |>
      dplyr::mutate(val_lt_null=dplyr::if_else(taf>=1,1,0)) |>
      dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

    p_val_2sid <- min(p_val_greater,p_val_less)*2

    p_value <-
      switch(alternative,
             two.sided = p_val_2sid,
             one.sided = p_val_less)

    # Final considerations ----
    return(list(
      p = p_value[[1]],
      alternative = alternative,
      taf_est = res_med_ci[[2,1]],
      taf_lower = res_med_ci[[1,1]],
      taf_upper = res_med_ci[[3,1]],
      CI_level = level,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    ))
  }
}
