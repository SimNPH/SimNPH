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

analyse_ma <- function(level=0.95, alternative = "two.sided",
                       dist=c("weibull","exponential",
                              "lognormal","gompertz",
                              "gamma"),
                       n_boot=1000,
                       values=c("median","6mo","12mo","24mo","aic_diff")) {

  # Checks -----
  stopifnot(alternative %in% c("two.sided", "one.sided"))
  stopifnot(all(dist %in% c("weibull","exponential",
                        "lognormal","gompertz",
                        "gamma")))
  stopifnot(all(values %in% c("median","6mo","12mo","24mo","aic_diff")))


  # wrap in SimDesign function
  function(condition, dat, fixed_objects = NULL){

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

    # Fitting the datasets to Weibul model ----
    if("weibull" %in% dist){

      FUN_fit_weibull <- function(X){
        survreg(Surv(t, evt)~1,X, dist="weibull")
      }

      ##control group
      fit_results_TRT0 <-
        lapply(bootstraped_survival_datasets_TRT0,
               FUN_fit_weibull)

      ##active group
      fit_results_TRT1 <-
        lapply(bootstraped_survival_datasets_TRT1,
               FUN_fit_weibull)

      AIC_values_weib_TRT0<-lapply(fit_results_TRT0, AIC)
      AIC_values_weib_TRT1<-lapply(fit_results_TRT1, AIC)
    }

    # Fitting the datasets to exponential model ----
    if("exponential" %in% dist){

      FUN_fit_exponential <- function(X){
        survreg(Surv(t, evt)~1,X,dist="exponential")
      }

      ## control group
      fit_results_exp_TRT0 <-
        lapply(bootstraped_survival_datasets_TRT0,
               FUN_fit_exponential)

      ##active group
      fit_results_exp_TRT1 <-
        lapply(bootstraped_survival_datasets_TRT1,
               FUN_fit_exponential)

      AIC_values_exp_TRT0<-lapply(fit_results_exp_TRT0, AIC)
      AIC_values_exp_TRT1<-lapply(fit_results_exp_TRT1, AIC)
    }

    # Fitting the datasets to lognormal model ----
    if("lognormal" %in% dist){

      FUN_fit_lognorm <- function(X){
        survreg(Surv(t, evt)~1,X,dist="lognormal")
      }

      ## control group
      fit_results_logn_TRT0 <-
        lapply(bootstraped_survival_datasets_TRT0,
               FUN_fit_lognorm)

      ##active group
      fit_results_logn_TRT1 <-
        lapply(bootstraped_survival_datasets_TRT1,
               FUN_fit_lognorm)

      AIC_values_logn_TRT0<-lapply(fit_results_logn_TRT0, AIC)
      AIC_values_logn_TRT1<-lapply(fit_results_logn_TRT1, AIC)
    }

    # Fitting the datasets to Gompertz model ----
    if("gompertz" %in% dist){

      FUN_fit_gom <- function(X){
        result <-flexsurv::flexsurvreg(Surv(t, evt)~1,data=X,dist="gom")
        return(result)
      }

      ## control group
      fit_results_gom_TRT0 <-
        lapply(bootstraped_survival_datasets_TRT0,
               FUN_fit_gom)

      ##active group
      fit_results_gom_TRT1 <-
        lapply(bootstraped_survival_datasets_TRT1,
               FUN_fit_gom)

      AIC_values_gom_TRT1 <- lapply(fit_results_gom_TRT1, AIC)
      AIC_values_gom_TRT0 <- lapply(fit_results_gom_TRT0, AIC)
    }

    # Fitting the datasets to Gamma model ----
    if("gamma" %in% dist){

      FUN_fit_gamma <- function(X){
        result <-flexsurv::flexsurvreg(Surv(t, evt)~1,data=X,dist="gamma")
        return(result)
      }

      ## control group
      fit_results_gamma_TRT0 <-
        lapply(bootstraped_survival_datasets_TRT0,
               FUN_fit_gamma)

      ##active group
      fit_results_gamma_TRT1 <-
        lapply(bootstraped_survival_datasets_TRT1,
               FUN_fit_gamma)

      AIC_values_gamma_TRT1 <- lapply(fit_results_gamma_TRT1, AIC)
      AIC_values_gamma_TRT0 <- lapply(fit_results_gamma_TRT0, AIC)
    }

    # Gather AIC values ----
    ### Combining the AIC values in a single dataset 'AIC_values_combined'
    AIC_values_combinedv2 <-
      tibble::tibble(AIC_exp_TRT0=unlist(AIC_values_exp_TRT0),
                     AIC_exp_TRT1=unlist(AIC_values_exp_TRT1),
                     AIC_logn_TRT0=unlist(AIC_values_logn_TRT0),
                     AIC_logn_TRT1=unlist(AIC_values_logn_TRT1),
                     AIC_weib_TRT0=unlist(AIC_values_weib_TRT0),
                     AIC_weib_TRT1=unlist(AIC_values_weib_TRT1),
                     AIC_gom_TRT0=unlist(AIC_values_gom_TRT0),
                     AIC_gom_TRT1=unlist(AIC_values_gom_TRT1),
                     AIC_gamma_TRT0=unlist(AIC_values_gamma_TRT0),
                     AIC_gamma_TRT1=unlist(AIC_values_gamma_TRT1))

    #### Adding exp, logn, gamma, gompertz, and weibull model AIC values together
    AIC_values_sumv2 <- AIC_values_combinedv2 |>
      dplyr::reframe(
        AIC_exp = (AIC_exp_TRT0 + AIC_exp_TRT1),
        AIC_logn =(AIC_logn_TRT0 + AIC_logn_TRT1),
        AIC_weib =(AIC_weib_TRT0 + AIC_weib_TRT1),
        AIC_gom = (AIC_gom_TRT0  + AIC_gom_TRT1),
        AIC_gamma = (AIC_gamma_TRT0 + AIC_gamma_TRT1)
      )

    ###obtaining the lowest AIC value for each bootstraped dataset
    AIC_values_sumv2 <- AIC_values_sumv2 |>
      dplyr::mutate(AIC_min=apply(AIC_values_sumv2, 1, FUN = min))


    ### Identifying the best fit model by substracting the min AIC value. The best fit model will have a 0 as their AIC value
    df_zeros<- AIC_values_sumv2 |>
      dplyr::mutate(AIC_exp = AIC_exp - AIC_min,
                    AIC_logn = AIC_logn - AIC_min,
                    AIC_weib = AIC_weib - AIC_min,
                    AIC_gom  = AIC_gom  - AIC_min,
                    AIC_gamma = AIC_gamma - AIC_min)


    ## Identifying the locations of zeros (corresponds to best fit models)
    zero_locations<-which(df_zeros == 0, arr.ind = TRUE)
    zero_locations_df <-data.frame(zero_locations)
    zero_locations_df<-zero_locations_df[order(zero_locations_df$row),]

    ### in best_fit_model_index 1-exponential, 2-lognormal, 3-weibull, 4-gompertz, 5 - gamma
    best_fit_model_index <- data.frame(zero_locations_df$col)

    # median time calculations ----
    if("median" %in% values){

      ### Weibul model
      if("weibull" %in% dist){

        #### Weibull distribution
        #### Median
        FUN_median <- function(fit_results){

          #   survreg's scale  =    1/(rweibull shape)
          #   survreg's intercept = log(rweibull scale)
          sh <- 1/fit_results$scale
          sc <- exp(fit_results$coefficients[[1]])

          med_tte<- qweibull(p = 0.5, shape = sh, scale = sc)
          return(med_tte)
        }

        ### active  group
        median_survival_TRT1<-lapply(fit_results_TRT1, FUN_median)

        ### control group
        median_survival_TRT0 <- lapply(fit_results_TRT0, FUN_median)

      }

      ### exponential model
      if("exponential" %in% dist){

        ####Exponential distribution
        ####Median
        FUN_median_exp <- function(X){
          ra <- (1/exp(X$coefficients[[1]]))
          med_tte <-qexp(p = 0.5, rate = ra)
          return(med_tte)
        }

        ### active group
        median_survival_exp_TRT1<-lapply(fit_results_exp_TRT1, FUN_median_exp)

        ### control group
        median_survival_exp_TRT0 <- lapply(fit_results_exp_TRT0, FUN_median_exp)

      }

      ### lognormal model
      if("lognormal" %in% dist){

        ####Lognormal distribution
        ####Median

        FUN_median_logn <- function(fit_results){
          mu <- fit_results$coefficients[[1]]
          sigma <- fit_results$scale
          med_tte <-qlnorm(p = 0.5, meanlog = mu, sdlog = sigma,lower.tail = TRUE, log.p = FALSE)
          return(med_tte)
        }

        ###control and active group
        median_survival_logn_TRT1<-lapply(fit_results_logn_TRT1, FUN_median_logn)
        median_survival_logn_TRT0 <- lapply(fit_results_logn_TRT0, FUN_median_logn)

      }

      ### Gompertz model
      if("gompertz" %in% dist){
        FUN_median_gom <- function(X){
          sh <- X$coefficients[[1]]
          ra <- exp(X$coefficients[[2]])
          med_tte<-flexsurv::qgompertz(p=0.5, shape = sh, rate = ra)
          return(med_tte)
        }

        ###control and active group
        median_survival_gom_TRT1<-lapply(fit_results_gom_TRT1, FUN_median_gom)
        median_survival_gom_TRT0 <- lapply(fit_results_gom_TRT0, FUN_median_gom)
      }

      ### Gamma model
      if("gamma" %in% dist){

        FUN_median_gamma <- function(X){
          ## coefficients are returned when X$coefficients are in log scale
          shape <-exp(X$coefficients[[1]])
          rate <-exp(X$coefficients[[2]])


          med_tte<- qgamma(p = 0.5, shape = shape, rate = rate)
          return(med_tte)
        }

        ###control and active group
        median_survival_gamma_TRT1<-lapply(fit_results_gamma_TRT1, FUN_median_gamma)
        median_survival_gamma_TRT0 <- lapply(fit_results_gamma_TRT0, FUN_median_gamma)
      }

      ### Gather all medians in one dataframe, and add them together
      All_medians <- tibble::tibble(
        logn_trt1_med = unlist(median_survival_logn_TRT1),
        exp_trt1_med  = unlist(median_survival_exp_TRT1),
        weib_trt1_med = unlist(median_survival_TRT1),
        gom_trt1_med = unlist(median_survival_gom_TRT1),
        gamma_trt1_med = unlist(median_survival_gamma_TRT1),
        exp_trt0_med  = unlist(median_survival_exp_TRT0),
        logn_trt0_med = unlist(median_survival_logn_TRT0),
        weib_trt0_med = unlist(median_survival_TRT0),
        gom_trt0_med = unlist(median_survival_gom_TRT0),
        gamma_trt0_med = unlist(median_survival_gamma_TRT0))


      Median_difference <- data.frame(Data_set_no = 1:n_boot, median_diff = NA)

      for (i in 1:n_boot){
        Median_difference$median_diff[i] <-
          ifelse(best_fit_model_index[i,] == 1,
                 All_medians$exp_trt1_med[i] - All_medians$exp_trt0_med[i],
                 ifelse(best_fit_model_index[i,] == 2,
                        All_medians$logn_trt1_med[i] - All_medians$logn_trt0_med[i],
                        ifelse(best_fit_model_index[i,] == 3,
                               All_medians$weib_trt1_med[i] - All_medians$weib_trt0_med[i],
                               ifelse(best_fit_model_index[i,] == 4,
                                      All_medians$gom_trt1_med[i] - All_medians$gom_trt0_med[i],
                                      All_medians$gamma_trt1_med[i] - All_medians$gamma_trt0_med[i] ))))
      }

      # compute confidence interval best fit model medians difference
      res_bestfit_ci <- Median_difference |>
        dplyr::reframe(qs=quantile(median_diff,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))
      # res_bestfit_ci

      # Statistcal significance of the median difference
      res_bestfit_sig_diff <- res_bestfit_ci |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))
      # res_bestfit_sig_diff

      p_val_greater <- Median_difference |>
        dplyr::mutate(val_gt_null=dplyr::if_else(median_diff<=0,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_less <- Median_difference |>
        dplyr::mutate(val_lt_null=dplyr::if_else(median_diff>=0,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_2sid <- min(p_val_greater,p_val_less)*2

      p_value <-
        switch(alternative,
               two.sided = p_val_2sid,
               one.sided = p_val_greater)

      # best fit model plots
      # p5 <- ggplot(Median_difference, aes(median_diff)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 0,col="red") +
      #   geom_vline(xintercept = res_bestfit_ci$qs)
      # p5

    }

    # 6mo survival calculations ----
    if("6mo" %in% values){

      ### Weibul model
      if("weibull" %in% dist){

        #### Weibull distribution
        #### survival probability
        FUN_surv_prob <- function(fit_results){

          #   survreg's scale  =    1/(rweibull shape)
          #   survreg's intercept = log(rweibull scale)
          sh <- 1/fit_results$scale
          sc <- exp(fit_results$coefficients[[1]])

          ch<- -pweibull(q = 180, shape = sh, scale = sc, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }



        weib_survprob_TRT1<-lapply(fit_results_TRT1, FUN_surv_prob)

        weib_survprob_TRT0 <- lapply(fit_results_TRT0, FUN_surv_prob)


      }

      ### exponential model
      if("exponential" %in% dist){

        ####Exponential distribution
        ####probabilities

        FUN_survprob_exp <- function(X){
          ra <- (1/exp(X$coefficients[[1]]))
          sh <- -pexp(q = 180, rate = ra, lower = FALSE, log = TRUE)

          surv_prob<- exp(-sh)
          return(surv_prob)
        }



        exp_survprob_TRT1 <- lapply(fit_results_exp_TRT1, FUN_survprob_exp)

        exp_survprob_TRT0 <- lapply(fit_results_exp_TRT0, FUN_survprob_exp)
      }

      ### lognormal model
      if("lognormal" %in% dist){

        FUN_survprob_logn <- function(fit_results){
          mu <- fit_results$coefficients[[1]]
          sigma <- fit_results$scale
          ch <- -plnorm(q = 180, meanlog = mu, sdlog = sigma, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        logn_survprob_TRT1<-lapply(fit_results_logn_TRT1, FUN_survprob_logn)
        logn_survprob_TRT0 <- lapply(fit_results_logn_TRT0, FUN_survprob_logn)
      }

      ### Gompertz model
      if("gompertz" %in% dist){
        FUN_survprob_gom <- function(X){
          sh <- X$coefficients[[1]]
          ra <- exp(X$coefficients[[2]])
          ch <- -flexsurv::pgompertz(q=180, shape = sh, rate = ra, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gom_survprob_TRT1<-lapply(fit_results_gom_TRT1, FUN_survprob_gom)
        gom_survprob_TRT0 <- lapply(fit_results_gom_TRT0, FUN_survprob_gom)
      }

      ### Gamma model
      if("gamma" %in% dist){

        FUN_survprob_gamma <- function(X){
          ## coefficients are returned when X$coefficients are in log scale
          shape <-exp(X$coefficients[[1]])
          rate <-exp(X$coefficients[[2]])


          ch <- -pgamma(q = 180, shape = shape, rate = rate, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gamma_survprob_TRT1 <- lapply(fit_results_gamma_TRT1, FUN_survprob_gamma)
        gamma_survprob_TRT0 <- lapply(fit_results_gamma_TRT0, FUN_survprob_gamma)
      }

      ### Gather all probabilities in one dataframe, and add them together
      All_probs <- tibble::tibble(
        logn_trt1_prob = unlist(logn_survprob_TRT1),
        exp_trt1_prob  = unlist(exp_survprob_TRT1),
        weib_trt1_prob = unlist(weib_survprob_TRT1),
        gom_trt1_prob = unlist(gom_survprob_TRT1),
        gamma_trt1_prob = unlist(gamma_survprob_TRT1),
        exp_trt0_prob  = unlist(exp_survprob_TRT0),
        logn_trt0_prob = unlist(logn_survprob_TRT0),
        weib_trt0_prob = unlist(weib_survprob_TRT0),
        gom_trt0_prob = unlist(gom_survprob_TRT0),
        gamma_trt0_prob = unlist(gamma_survprob_TRT0))

      # compute the probability difference of the best fit model
      prob_difference <- data.frame(Data_set_no = 1:n_boot, prob_diff = NA)
      prob_ratio <- data.frame(Data_set_no = 1:n_boot, prob_ratio = NA)

      for (i in 1:n_boot){
        prob_difference$prob_diff[i] <-ifelse(best_fit_model_index[i,] == 1,All_probs$exp_trt1_prob[i] - All_probs$exp_trt0_prob[i],
                                              ifelse(best_fit_model_index[i,] == 2, All_probs$logn_trt1_prob[i] - All_probs$logn_trt0_prob[i],
                                                     ifelse(best_fit_model_index[i,] == 3,All_probs$weib_trt1_prob[i] - All_probs$weib_trt0_prob[i],
                                                            ifelse(best_fit_model_index[i,] == 4, All_probs$gom_trt1_prob[i] - All_probs$gom_trt0_prob[i], All_probs$gamma_trt1_prob[i] - All_probs$gamma_trt0_prob[i] ))))

        prob_ratio$prob_ratio[i] <-
          switch(best_fit_model_index[i,],
                 '1' = All_probs$exp_trt1_prob[i]/All_probs$exp_trt0_prob[i],
                 '2' = All_probs$logn_trt1_prob[i]/All_probs$logn_trt0_prob[i],
                 '3' = All_probs$weib_trt1_prob[i]/All_probs$weib_trt0_prob[i],
                 '4' = All_probs$gom_trt1_prob[i]/All_probs$gom_trt0_prob[i],
                 '5' = All_probs$gamma_trt1_prob[i]/All_probs$gamma_trt0_prob[i])

      }

      # compute confidence interval for best fit model probability difference
      res_bestfit_ci_6m <- prob_difference |>
        dplyr::reframe(qs=quantile(prob_diff,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))

      res_bestfit_ci_ratio_6m <- prob_ratio |>
        dplyr::reframe(qs=quantile(prob_ratio,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))
      # res_bestfit_ci

      # Statistcal significance of the median difference
      res_bestfit_sig_diff_6m <- res_bestfit_ci_6m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))

      res_bestfit_sig_ratio_6m <- res_bestfit_ci_ratio_6m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))
      # res_bestfit_sig_diff

      p_val_greater_6m <- prob_difference |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_diff<=0,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_greater_ratio_6m <- prob_ratio |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_ratio<=1,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_less_6m <- prob_difference |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_diff>=0,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_less_ratio_6m <- prob_ratio |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_ratio>=1,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_2sid_6m <- min(p_val_greater_6m,p_val_less_6m)*2
      p_val_2sid_ratio_6m <- min(p_val_greater_ratio_6m,p_val_less_ratio_6m)*2

      p_value_6m <-
        switch(alternative,
               two.sided = p_val_2sid_6m,
               one.sided = p_val_greater_6m)

      p_value_ratio_6m <-
        switch(alternative,
               two.sided = p_val_2sid_ratio_6m,
               one.sided = p_val_greater_ratio_6m)

      # best fit model plots
      # p5 <- ggplot(Median_difference, aes(median_diff)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 0,col="red") +
      #   geom_vline(xintercept = res_bestfit_ci$qs)
      # p5

      # p5 <- ggplot(prob_ratio, aes(prob_ratio)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 1,col="red") +
      #   geom_vline(xintercept = res_bestfit_ci_ratio_6m$qs)
      # p5

      #
    }

    # 12mo survival calculations ----
    if("12mo" %in% values){

      ### Weibul model
      if("weibull" %in% dist){

        #### Weibull distribution
        #### survival probability
        FUN_surv_prob <- function(fit_results){

          #   survreg's scale  =    1/(rweibull shape)
          #   survreg's intercept = log(rweibull scale)
          sh <- 1/fit_results$scale
          sc <- exp(fit_results$coefficients[[1]])

          ch<- -pweibull(q = 365, shape = sh, scale = sc, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }



        weib_survprob_TRT1<-lapply(fit_results_TRT1, FUN_surv_prob)

        weib_survprob_TRT0 <- lapply(fit_results_TRT0, FUN_surv_prob)


      }

      ### exponential model
      if("exponential" %in% dist){

        ####Exponential distribution
        ####probabilities

        FUN_survprob_exp <- function(X){
          ra <- (1/exp(X$coefficients[[1]]))
          sh <- -pexp(q = 365, rate = ra, lower = FALSE, log = TRUE)

          surv_prob<- exp(-sh)
          return(surv_prob)
        }



        exp_survprob_TRT1 <- lapply(fit_results_exp_TRT1, FUN_survprob_exp)

        exp_survprob_TRT0 <- lapply(fit_results_exp_TRT0, FUN_survprob_exp)
      }

      ### lognormal model
      if("lognormal" %in% dist){

        FUN_survprob_logn <- function(fit_results){
          mu <- fit_results$coefficients[[1]]
          sigma <- fit_results$scale
          ch <- -plnorm(q = 365, meanlog = mu, sdlog = sigma, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        logn_survprob_TRT1<-lapply(fit_results_logn_TRT1, FUN_survprob_logn)
        logn_survprob_TRT0 <- lapply(fit_results_logn_TRT0, FUN_survprob_logn)
      }

      ### Gompertz model
      if("gompertz" %in% dist){
        FUN_survprob_gom <- function(X){
          sh <- X$coefficients[[1]]
          ra <- exp(X$coefficients[[2]])
          ch <- -flexsurv::pgompertz(q=365, shape = sh, rate = ra, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gom_survprob_TRT1<-lapply(fit_results_gom_TRT1, FUN_survprob_gom)
        gom_survprob_TRT0 <- lapply(fit_results_gom_TRT0, FUN_survprob_gom)
      }

      ### Gamma model
      if("gamma" %in% dist){

        FUN_survprob_gamma <- function(X){
          ## coefficients are returned when X$coefficients are in log scale
          shape <-exp(X$coefficients[[1]])
          rate <-exp(X$coefficients[[2]])


          ch <- -pgamma(q = 365, shape = shape, rate = rate, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gamma_survprob_TRT1 <- lapply(fit_results_gamma_TRT1, FUN_survprob_gamma)
        gamma_survprob_TRT0 <- lapply(fit_results_gamma_TRT0, FUN_survprob_gamma)
      }

      ### Gather all probabilities in one dataframe, and add them together
      All_probs <- tibble::tibble(
        logn_trt1_prob = unlist(logn_survprob_TRT1),
        exp_trt1_prob  = unlist(exp_survprob_TRT1),
        weib_trt1_prob = unlist(weib_survprob_TRT1),
        gom_trt1_prob = unlist(gom_survprob_TRT1),
        gamma_trt1_prob = unlist(gamma_survprob_TRT1),
        exp_trt0_prob  = unlist(exp_survprob_TRT0),
        logn_trt0_prob = unlist(logn_survprob_TRT0),
        weib_trt0_prob = unlist(weib_survprob_TRT0),
        gom_trt0_prob = unlist(gom_survprob_TRT0),
        gamma_trt0_prob = unlist(gamma_survprob_TRT0))

      # compute the probability difference of the best fit model
      prob_difference_12m <- data.frame(Data_set_no = 1:n_boot, prob_diff = NA)
      prob_ratio_12m <- data.frame(Data_set_no = 1:n_boot, prob_ratio = NA)

      for (i in 1:n_boot){
        prob_difference_12m$prob_diff[i] <-ifelse(best_fit_model_index[i,] == 1,All_probs$exp_trt1_prob[i] - All_probs$exp_trt0_prob[i],
                                              ifelse(best_fit_model_index[i,] == 2, All_probs$logn_trt1_prob[i] - All_probs$logn_trt0_prob[i],
                                                     ifelse(best_fit_model_index[i,] == 3,All_probs$weib_trt1_prob[i] - All_probs$weib_trt0_prob[i],
                                                            ifelse(best_fit_model_index[i,] == 4, All_probs$gom_trt1_prob[i] - All_probs$gom_trt0_prob[i], All_probs$gamma_trt1_prob[i] - All_probs$gamma_trt0_prob[i] ))))

        prob_ratio_12m$prob_ratio[i] <-
          switch(best_fit_model_index[i,],
                 '1' = All_probs$exp_trt1_prob[i]/All_probs$exp_trt0_prob[i],
                 '2' = All_probs$logn_trt1_prob[i]/All_probs$logn_trt0_prob[i],
                 '3' = All_probs$weib_trt1_prob[i]/All_probs$weib_trt0_prob[i],
                 '4' = All_probs$gom_trt1_prob[i]/All_probs$gom_trt0_prob[i],
                 '5' = All_probs$gamma_trt1_prob[i]/All_probs$gamma_trt0_prob[i])

      }

      # compute confidence interval for best fit model probability difference
      res_bestfit_ci_12m <- prob_difference_12m |>
        dplyr::reframe(qs=quantile(prob_diff,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))
      # res_bestfit_ci

      res_bestfit_ci_ratio_12m <- prob_ratio_12m |>
        dplyr::reframe(qs=quantile(prob_ratio,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))

      # Statistcal significance of the median difference
      res_bestfit_sig_diff_12m <- res_bestfit_ci_12m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))
      # res_bestfit_sig_diff

      res_bestfit_sig_ratio_12m <- res_bestfit_ci_ratio_12m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))

      p_val_greater_12m <- prob_difference_12m |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_diff<=0,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_greater_ratio_12m <- prob_ratio_12m |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_ratio<=1,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_less_12m <- prob_difference_12m |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_diff>=0,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_less_ratio_12m <- prob_ratio_12m |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_ratio>=1,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_2sid_12m <- min(p_val_greater_12m,p_val_less_12m)*2
      p_val_2sid_ratio_12m <- min(p_val_greater_ratio_12m,p_val_less_ratio_12m)*2

      p_value_12m <-
        switch(alternative,
               two.sided = p_val_2sid_12m,
               one.sided = p_val_greater_12m)

      p_value_ratio_12m <-
        switch(alternative,
               two.sided = p_val_2sid_ratio_12m,
               one.sided = p_val_greater_ratio_12m)

      # best fit model plots
      # p5 <- ggplot(Median_difference, aes(median_diff)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 0,col="red") +
      #   geom_vline(xintercept = res_bestfit_ci$qs)
      # p5
      #

    }

    # 24mo survival calculations ----
    if("24mo" %in% values){

      ## Weibul model
      if("weibull" %in% dist){

        #### Weibull distribution
        #### survival probability
        FUN_surv_prob <- function(fit_results){

          #   survreg's scale  =    1/(rweibull shape)
          #   survreg's intercept = log(rweibull scale)
          sh <- 1/fit_results$scale
          sc <- exp(fit_results$coefficients[[1]])

          ch<- -pweibull(q = 730, shape = sh, scale = sc, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }



        weib_survprob_TRT1<-lapply(fit_results_TRT1, FUN_surv_prob)

        weib_survprob_TRT0 <- lapply(fit_results_TRT0, FUN_surv_prob)


      }

      ### exponential model
      if("exponential" %in% dist){

        ####Exponential distribution
        ####probabilities

        FUN_survprob_exp <- function(X){
          ra <- (1/exp(X$coefficients[[1]]))
          sh <- -pexp(q = 730, rate = ra, lower = FALSE, log = TRUE)

          surv_prob<- exp(-sh)
          return(surv_prob)
        }



        exp_survprob_TRT1 <- lapply(fit_results_exp_TRT1, FUN_survprob_exp)

        exp_survprob_TRT0 <- lapply(fit_results_exp_TRT0, FUN_survprob_exp)
      }

      ### lognormal model
      if("lognormal" %in% dist){

        FUN_survprob_logn <- function(fit_results){
          mu <- fit_results$coefficients[[1]]
          sigma <- fit_results$scale
          ch <- -plnorm(q = 730, meanlog = mu, sdlog = sigma, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        logn_survprob_TRT1<-lapply(fit_results_logn_TRT1, FUN_survprob_logn)
        logn_survprob_TRT0 <- lapply(fit_results_logn_TRT0, FUN_survprob_logn)
      }

      ### Gompertz model
      if("gompertz" %in% dist){
        FUN_survprob_gom <- function(X){
          sh <- X$coefficients[[1]]
          ra <- exp(X$coefficients[[2]])
          ch <- -flexsurv::pgompertz(q=730, shape = sh, rate = ra, lower = FALSE, log = TRUE)
          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gom_survprob_TRT1<-lapply(fit_results_gom_TRT1, FUN_survprob_gom)
        gom_survprob_TRT0 <- lapply(fit_results_gom_TRT0, FUN_survprob_gom)
      }

      ### Gamma model
      if("gamma" %in% dist){

        FUN_survprob_gamma <- function(X){
          ## coefficients are returned when X$coefficients are in log scale
          shape <-exp(X$coefficients[[1]])
          rate <-exp(X$coefficients[[2]])


          ch <- -pgamma(q = 730, shape = shape, rate = rate, lower = FALSE, log = TRUE)

          surv_prob <- exp(-ch)
          return(surv_prob)
        }

        gamma_survprob_TRT1 <- lapply(fit_results_gamma_TRT1, FUN_survprob_gamma)
        gamma_survprob_TRT0 <- lapply(fit_results_gamma_TRT0, FUN_survprob_gamma)
      }

      ### Gather all probabilities in one dataframe, and add them together
      All_probs <- tibble::tibble(
        logn_trt1_prob = unlist(logn_survprob_TRT1),
        exp_trt1_prob  = unlist(exp_survprob_TRT1),
        weib_trt1_prob = unlist(weib_survprob_TRT1),
        gom_trt1_prob = unlist(gom_survprob_TRT1),
        gamma_trt1_prob = unlist(gamma_survprob_TRT1),
        exp_trt0_prob  = unlist(exp_survprob_TRT0),
        logn_trt0_prob = unlist(logn_survprob_TRT0),
        weib_trt0_prob = unlist(weib_survprob_TRT0),
        gom_trt0_prob = unlist(gom_survprob_TRT0),
        gamma_trt0_prob = unlist(gamma_survprob_TRT0))

      # compute the probability difference of the best fit model
      prob_difference_24m <- data.frame(Data_set_no = 1:n_boot, prob_diff = NA)
      prob_ratio_24m <- data.frame(Data_set_no = 1:n_boot, prob_ratio = NA)

      for (i in 1:n_boot){
        prob_difference_24m$prob_diff[i] <-ifelse(best_fit_model_index[i,] == 1,All_probs$exp_trt1_prob[i] - All_probs$exp_trt0_prob[i],
                                                  ifelse(best_fit_model_index[i,] == 2, All_probs$logn_trt1_prob[i] - All_probs$logn_trt0_prob[i],
                                                         ifelse(best_fit_model_index[i,] == 3,All_probs$weib_trt1_prob[i] - All_probs$weib_trt0_prob[i],
                                                                ifelse(best_fit_model_index[i,] == 4, All_probs$gom_trt1_prob[i] - All_probs$gom_trt0_prob[i], All_probs$gamma_trt1_prob[i] - All_probs$gamma_trt0_prob[i] ))))
        prob_ratio_24m$prob_ratio[i] <-
          switch(best_fit_model_index[i,],
                 '1' = All_probs$exp_trt1_prob[i]/All_probs$exp_trt0_prob[i],
                 '2' = All_probs$logn_trt1_prob[i]/All_probs$logn_trt0_prob[i],
                 '3' = All_probs$weib_trt1_prob[i]/All_probs$weib_trt0_prob[i],
                 '4' = All_probs$gom_trt1_prob[i]/All_probs$gom_trt0_prob[i],
                 '5' = All_probs$gamma_trt1_prob[i]/All_probs$gamma_trt0_prob[i])
      }

      # compute confidence interval for best fit model probability difference
      res_bestfit_ci_24m <- prob_difference_24m |>
        dplyr::reframe(qs=quantile(prob_diff,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))
      # res_bestfit_ci

      res_bestfit_ci_ratio_24m <- prob_ratio_24m |>
        dplyr::reframe(qs=quantile(prob_ratio,
                                   c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))

      # Statistcal significance of the median difference
      res_bestfit_sig_diff_24m <- res_bestfit_ci_24m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))
      # res_bestfit_sig_diff

      res_bestfit_sig_ratio_24m <- res_bestfit_ci_ratio_24m |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))

      p_val_greater_24m <- prob_difference_24m |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_diff<=0,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_greater_ratio_24m <- prob_ratio_24m |>
        dplyr::mutate(val_gt_null=dplyr::if_else(prob_ratio<=1,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_less_24m <- prob_difference_24m |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_diff>=0,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_less_ratio_24m <- prob_ratio_24m |>
        dplyr::mutate(val_lt_null=dplyr::if_else(prob_ratio>=1,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_2sid_24m <- min(p_val_greater_24m,p_val_less_24m)*2
      p_val_2sid_ratio_24m <- min(p_val_greater_ratio_24m,p_val_less_ratio_24m)*2

      p_value_24m <-
        switch(alternative,
               two.sided = p_val_2sid_24m,
               one.sided = p_val_greater_24m)

      p_value_ratio_24m <-
        switch(alternative,
               two.sided = p_val_2sid_ratio_24m,
               one.sided = p_val_greater_ratio_24m)

      # best fit model plots
      # p5 <- ggplot(Median_difference, aes(median_diff)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 0,col="red") +
      #   geom_vline(xintercept = res_bestfit_ci$qs)
      # p5
      #

    }

    # AIC weighted difference calculations ----
    if("aic_diff" %in% values){

      ### merge the lists containing bootstraped datasets from arm TRT0 and TRT1
      FUN_merge_lists <- function(X,Y) {
        merged_df<-rbind(X, Y)
        return(merged_df)
      }

      TRT_merged <- mapply(bootstraped_survival_datasets_TRT0,
                           bootstraped_survival_datasets_TRT1,
                           FUN = FUN_merge_lists, SIMPLIFY = FALSE)

      ### fit the merged datasets with model which was previosuly defined as the best fit

      FUN_fit_merged <- function(TRT_merged,best_fit_model_index)
      {
        ### index key: 1 - exponential; 2 - lognormal; 3 - weibull; 4 - gompertz
        distribution <-
          switch(best_fit_model_index,
                 '1' = 'exp',
                 '2'= 'lognormal',
                 '3' = 'weibull',
                 '4' = 'gompertz',
                 '5' = 'gamma')

        ### gompertz and gamma use flexsurv instead of survreg
        fit_results <- if(best_fit_model_index == 4 | best_fit_model_index == 5){
          flexsurv::flexsurvreg(Surv(t, evt)~1,data=TRT_merged,dist=distribution)
        } else {
          survreg(Surv(t, evt)~1,data=TRT_merged,dist=distribution)
        }
        return(fit_results)
      }

      TRT_merged_fit <-
        mapply(TRT_merged,
               as.list(c(best_fit_model_index[[1]])),
               FUN = FUN_fit_merged,
               SIMPLIFY = FALSE)
      AIC_values_merged_TRT<-lapply(TRT_merged_fit, AIC)

      ## Add separate bestfit min AIC values and combined fit AIC values to a single data frame
      AIC_values_sep_n_comb <-
        tibble::tibble(AIC_sep_fit = AIC_values_sumv2$AIC_min,
                       AIC_comb_fit = unlist(AIC_values_merged_TRT))

      ## Obtain the difference of AIC values
      AIC_values_sep_n_comb <- AIC_values_sep_n_comb |>
        dplyr::mutate(AIC_diff = AIC_sep_fit - AIC_comb_fit)

      # compute confidence interval for AIC difference
      res_AIC_diff_ci <- AIC_values_sep_n_comb |>
        dplyr::reframe(qs=quantile(AIC_diff, c((1-level)/2, level+(1-level)/2)),
                       prob = c((1-level)/2, level+(1-level)/2))

      ####statistcal significance of the AIC difference
      res_AIC_diff_sig_diff <- res_AIC_diff_ci |>
        dplyr::mutate(across(everything(),sign)) |>
        dplyr::reframe(sig_diff=sum(qs)) |>
        dplyr::mutate(sig_diff=dplyr::if_else(sig_diff==0,0,1))

       p_val_greater_AIC_diff <- AIC_values_sep_n_comb |>
        dplyr::mutate(val_gt_null=dplyr::if_else(AIC_diff<=0,1,0)) |>
        dplyr::reframe(p_val_greater=(sum(val_gt_null)/dplyr::n()))

      p_val_less_AIC_diff <- AIC_values_sep_n_comb |>
        dplyr::mutate(val_lt_null=dplyr::if_else(AIC_diff>=0,1,0)) |>
        dplyr::reframe(p_val_less=(sum(val_lt_null)/dplyr::n()))

      p_val_2sid_AIC_diff <- min(p_val_greater_AIC_diff,p_val_less_AIC_diff)*2

      p_value_AIC_diff <-
        switch(alternative,
               two.sided = p_val_2sid_AIC_diff,
               one.sided = p_val_less_AIC_diff)

      #### weighted model plots
      # p1 <- ggplot(AIC_values_sep_n_comb, aes(AIC_diff)) +
      #   geom_histogram() +
      #   geom_vline(xintercept = 0,col="red") +
      #   geom_vline(xintercept = res_AIC_diff_ci$qs)
      # p1
    }

    # Final considerations ----
    return(list(
      p = p_value[[1]],

      p_6m_diff = p_value_6m[[1]],
      p_6m_ratio = p_value_ratio_6m[[1]],

      p_12m_diff = p_value_12m[[1]],
      p_12m_ratio = p_value_ratio_12m[[1]],

      p_24m_diff = p_value_24m[[1]],
      p_24m_ratio = p_value_ratio_24m[[1]],

      p_AIC_diff = p_value_AIC_diff[[1]],
      alternative = alternative,
      med_diff_est = median(Median_difference$median_diff),
      med_diff_lower = res_bestfit_ci[1,1],
      med_diff_upper = res_bestfit_ci[2,1],

      surv_6m_diff_est = median(prob_difference$prob_diff),
      surv_6m_diff_lower = res_bestfit_ci_6m[1,1],
      surv_6m_diff_upper = res_bestfit_ci_6m[2,1],

      surv_6m_ratio_est = median(prob_ratio$prob_ratio),
      surv_6m_ratio_lower = res_bestfit_ci_ratio_6m[1,1],
      surv_6m_ratio_upper = res_bestfit_ci_ratio_6m[2,1],

      surv_12m_diff_est = median(prob_difference_12m$prob_diff),
      surv_12m_diff_lower = res_bestfit_ci_12m[1,1],
      surv_12m_diff_upper = res_bestfit_ci_12m[2,1],

      surv_12m_ratio_est = median(prob_ratio_12m$prob_ratio),
      surv_12m_ratio_lower = res_bestfit_ci_ratio_12m[1,1],
      surv_12m_ratio_upper = res_bestfit_ci_ratio_12m[2,1],

      surv_24m_diff_est = median(prob_difference_24m$prob_diff),
      surv_24m_diff_lower = res_bestfit_ci_24m[1,1],
      surv_24m_diff_upper = res_bestfit_ci_24m[2,1],

      surv_24m_ratio_est = median(prob_ratio_24m$prob_ratio),
      surv_24m_ratio_lower = res_bestfit_ci_ratio_24m[1,1],
      surv_24m_ratio_upper = res_bestfit_ci_ratio_24m[2,1],

      AIC_diff_est = median(AIC_values_sep_n_comb$AIC_diff),
      AIC_diff_lower = unname(res_AIC_diff_ci[1,1][[1]]),
      AIC_diff_upper = unname(res_AIC_diff_ci[2,1][[1]]),
      sig_diff = res_bestfit_sig_diff$sig_diff,
      sig_diff_6m_diff = res_bestfit_sig_diff_6m$sig_diff,
      sig_diff_6m_ratio = res_bestfit_sig_ratio_6m$sig_diff,

      sig_diff_12m_diff = res_bestfit_sig_diff_12m$sig_diff,
      sig_diff_12m_ratio = res_bestfit_sig_ratio_12m$sig_diff,

      sig_diff_24m = res_bestfit_sig_diff_24m$sig_diff,
      sig_diff_24m_ratio = res_bestfit_sig_ratio_24m$sig_diff,

      sig_diff_AIC = res_AIC_diff_sig_diff$sig_diff,

      CI_level = level,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    ))
  }
}
