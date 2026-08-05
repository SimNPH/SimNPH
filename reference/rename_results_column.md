# Rename Columns in Simulation Results and Update Attributes

Rename Columns in Simulation Results and Update Attributes

## Usage

``` r
rename_results_column(results, rename)

rename_results_column_pattern(results, pattern, replacement)
```

## Arguments

- results:

  `SimDesign` object

- rename:

  named vector of new names

- pattern:

  regexp pattern as understood by
  [`stringr::str_replace_all`](https://stringr.tidyverse.org/reference/str_replace.html)

- replacement:

  replacement as understood by
  [`stringr::str_replace_all`](https://stringr.tidyverse.org/reference/str_replace.html)

## Value

`SimDesign` object with updated column names

## Functions

- `rename_results_column()`: Rename Columns in Simulation Results

- `rename_results_column_pattern()`: Rename Columns in Simulation
  Results by Pattern

## Examples

``` r
# \donttest{
condition <- merge(
assumptions_delayed_effect(),
design_fixed_followup(),
by=NULL
) |>
  tail(4) |>
  true_summary_statistics_delayed_effect(cutoff_stats = 15)

sim_results <- runSimulation(
  design=condition,
  replications=10,
  generate=generate_delayed_effect,
  analyse=list(
    logrank  = analyse_logrank(alternative = "one.sided"),
    mwlrt = analyse_modelstly_weighted(t_star = m2d(24))
  ),
  summarise = create_summarise_function(
    logrank = summarise_test(0.025),
    mwlrt = summarise_test(0.025)
  )
)
#> 
#> Design: 1/4;   Replications: 10;   RAM Used: 201.3 Mb;   Total Time: 0.00s 
#>  Conditions: delay=121., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=1034, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 2/4;   Replications: 10;   RAM Used: 201.3 Mb;   Total Time: 0.21s 
#>  Conditions: delay=182., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=1004, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 3/4;   Replications: 10;   RAM Used: 201.3 Mb;   Total Time: 0.43s 
#>  Conditions: delay=243., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=974, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 4/4;   Replications: 10;   RAM Used: 201.3 Mb;   Total Time: 0.64s 
#>  Conditions: delay=304., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=943., mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Simulation complete. Total execution time: 0.85s

names(sim_results)
#>  [1] "delay"                   "hazard_ctrl"            
#>  [3] "hazard_trt"              "random_withdrawal"      
#>  [5] "n_trt"                   "n_ctrl"                 
#>  [7] "followup"                "recruitment"            
#>  [9] "median_survival_trt"     "median_survival_ctrl"   
#> [11] "rmst_trt_15"             "rmst_ctrl_15"           
#> [13] "gAHR_15"                 "AHR_15"                 
#> [15] "AHRoc_15"                "AHRoc_robust_15"        
#> [17] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#> [19] "logrank.N"               "logrank.mean_n_pat"     
#> [21] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#> [23] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#> [25] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [27] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [29] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [31] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [33] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> [35] "REPLICATIONS"            "SIM_TIME"               
#> [37] "RAM_USED"                "SEED"                   
#> [39] "COMPLETED"              
attr(sim_results, "design_names")
#> $design
#>  [1] "delay"                "hazard_ctrl"          "hazard_trt"          
#>  [4] "random_withdrawal"    "n_trt"                "n_ctrl"              
#>  [7] "followup"             "recruitment"          "median_survival_trt" 
#> [10] "median_survival_ctrl" "rmst_trt_15"          "rmst_ctrl_15"        
#> [13] "gAHR_15"              "AHR_15"               "AHRoc_15"            
#> [16] "AHRoc_robust_15"     
#> 
#> $sim
#>  [1] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#>  [3] "logrank.N"               "logrank.mean_n_pat"     
#>  [5] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#>  [7] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#>  [9] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [11] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [13] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [15] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [17] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> 
#> $bootCI
#> character(0)
#> 
#> $extra
#> [1] "REPLICATIONS" "SIM_TIME"     "RAM_USED"     "SEED"         "COMPLETED"   
#> 
#> $errors
#> [1] "ERRORS"
#> 
#> $warnings
#> [1] "WARNINGS"
#> 

sim_results <- sim_results |>
  rename_results_column(c("delay"="onset"))

names(sim_results)
#>  [1] "onset"                   "hazard_ctrl"            
#>  [3] "hazard_trt"              "random_withdrawal"      
#>  [5] "n_trt"                   "n_ctrl"                 
#>  [7] "followup"                "recruitment"            
#>  [9] "median_survival_trt"     "median_survival_ctrl"   
#> [11] "rmst_trt_15"             "rmst_ctrl_15"           
#> [13] "gAHR_15"                 "AHR_15"                 
#> [15] "AHRoc_15"                "AHRoc_robust_15"        
#> [17] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#> [19] "logrank.N"               "logrank.mean_n_pat"     
#> [21] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#> [23] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#> [25] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [27] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [29] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [31] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [33] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> [35] "REPLICATIONS"            "SIM_TIME"               
#> [37] "RAM_USED"                "SEED"                   
#> [39] "COMPLETED"              
attr(sim_results, "design_names")
#> $design
#>  [1] "onset"                "hazard_ctrl"          "hazard_trt"          
#>  [4] "random_withdrawal"    "n_trt"                "n_ctrl"              
#>  [7] "followup"             "recruitment"          "median_survival_trt" 
#> [10] "median_survival_ctrl" "rmst_trt_15"          "rmst_ctrl_15"        
#> [13] "gAHR_15"              "AHR_15"               "AHRoc_15"            
#> [16] "AHRoc_robust_15"     
#> 
#> $sim
#>  [1] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#>  [3] "logrank.N"               "logrank.mean_n_pat"     
#>  [5] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#>  [7] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#>  [9] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [11] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [13] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [15] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [17] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> 
#> $bootCI
#> character(0)
#> 
#> $extra
#> [1] "REPLICATIONS" "SIM_TIME"     "RAM_USED"     "SEED"         "COMPLETED"   
#> 
#> $errors
#> [1] "ERRORS"
#> 
#> $warnings
#> [1] "WARNINGS"
#> 
# }
# \donttest{
  condition <- merge(
    assumptions_delayed_effect(),
    design_fixed_followup(),
    by=NULL
  ) |>
    tail(4) |>
    true_summary_statistics_delayed_effect(cutoff_stats = 15)

  sim_results <- runSimulation(
    design=condition,
    replications=10,
    generate=generate_delayed_effect,
    analyse=list(
      logrank  = analyse_logrank(alternative = "one.sided"),
      mwlrt = analyse_modelstly_weighted(t_star = m2d(24))
    ),
    summarise = create_summarise_function(
      logrank = summarise_test(0.025),
      mwlrt = summarise_test(0.025)
    )
  )
#> 
#> Design: 1/4;   Replications: 10;   RAM Used: 201.5 Mb;   Total Time: 0.00s 
#>  Conditions: delay=121., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=1034, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 2/4;   Replications: 10;   RAM Used: 201.4 Mb;   Total Time: 0.21s 
#>  Conditions: delay=182., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=1004, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 3/4;   Replications: 10;   RAM Used: 201.4 Mb;   Total Time: 0.42s 
#>  Conditions: delay=243., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=974, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Design: 4/4;   Replications: 10;   RAM Used: 201.4 Mb;   Total Time: 0.63s 
#>  Conditions: delay=304., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=943., mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Simulation complete. Total execution time: 0.84s

  names(sim_results)
#>  [1] "delay"                   "hazard_ctrl"            
#>  [3] "hazard_trt"              "random_withdrawal"      
#>  [5] "n_trt"                   "n_ctrl"                 
#>  [7] "followup"                "recruitment"            
#>  [9] "median_survival_trt"     "median_survival_ctrl"   
#> [11] "rmst_trt_15"             "rmst_ctrl_15"           
#> [13] "gAHR_15"                 "AHR_15"                 
#> [15] "AHRoc_15"                "AHRoc_robust_15"        
#> [17] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#> [19] "logrank.N"               "logrank.mean_n_pat"     
#> [21] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#> [23] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#> [25] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [27] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [29] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [31] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [33] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> [35] "REPLICATIONS"            "SIM_TIME"               
#> [37] "RAM_USED"                "SEED"                   
#> [39] "COMPLETED"              
  attr(sim_results, "design_names")
#> $design
#>  [1] "delay"                "hazard_ctrl"          "hazard_trt"          
#>  [4] "random_withdrawal"    "n_trt"                "n_ctrl"              
#>  [7] "followup"             "recruitment"          "median_survival_trt" 
#> [10] "median_survival_ctrl" "rmst_trt_15"          "rmst_ctrl_15"        
#> [13] "gAHR_15"              "AHR_15"               "AHRoc_15"            
#> [16] "AHRoc_robust_15"     
#> 
#> $sim
#>  [1] "logrank.rejection_0.025" "logrank.N_missing_0.025"
#>  [3] "logrank.N"               "logrank.mean_n_pat"     
#>  [5] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#>  [7] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#>  [9] "logrank.N_missing_n_evt" "mwlrt.rejection_0.025"  
#> [11] "mwlrt.N_missing_0.025"   "mwlrt.N"                
#> [13] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [15] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [17] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> 
#> $bootCI
#> character(0)
#> 
#> $extra
#> [1] "REPLICATIONS" "SIM_TIME"     "RAM_USED"     "SEED"         "COMPLETED"   
#> 
#> $errors
#> [1] "ERRORS"
#> 
#> $warnings
#> [1] "WARNINGS"
#> 

  sim_results <- sim_results |>
    rename_results_column_pattern(pattern = "_0.025", replacement = "")

  names(sim_results)
#>  [1] "delay"                   "hazard_ctrl"            
#>  [3] "hazard_trt"              "random_withdrawal"      
#>  [5] "n_trt"                   "n_ctrl"                 
#>  [7] "followup"                "recruitment"            
#>  [9] "median_survival_trt"     "median_survival_ctrl"   
#> [11] "rmst_trt_15"             "rmst_ctrl_15"           
#> [13] "gAHR_15"                 "AHR_15"                 
#> [15] "AHRoc_15"                "AHRoc_robust_15"        
#> [17] "logrank.rejection"       "logrank.N_missing"      
#> [19] "logrank.N"               "logrank.mean_n_pat"     
#> [21] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#> [23] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#> [25] "logrank.N_missing_n_evt" "mwlrt.rejection"        
#> [27] "mwlrt.N_missing"         "mwlrt.N"                
#> [29] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [31] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [33] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> [35] "REPLICATIONS"            "SIM_TIME"               
#> [37] "RAM_USED"                "SEED"                   
#> [39] "COMPLETED"              
  attr(sim_results, "design_names")
#> $design
#>  [1] "delay"                "hazard_ctrl"          "hazard_trt"          
#>  [4] "random_withdrawal"    "n_trt"                "n_ctrl"              
#>  [7] "followup"             "recruitment"          "median_survival_trt" 
#> [10] "median_survival_ctrl" "rmst_trt_15"          "rmst_ctrl_15"        
#> [13] "gAHR_15"              "AHR_15"               "AHRoc_15"            
#> [16] "AHRoc_robust_15"     
#> 
#> $sim
#>  [1] "logrank.rejection"       "logrank.N_missing"      
#>  [3] "logrank.N"               "logrank.mean_n_pat"     
#>  [5] "logrank.sd_n_pat"        "logrank.mean_n_evt"     
#>  [7] "logrank.sd_n_evt"        "logrank.N_missing_n_pat"
#>  [9] "logrank.N_missing_n_evt" "mwlrt.rejection"        
#> [11] "mwlrt.N_missing"         "mwlrt.N"                
#> [13] "mwlrt.mean_n_pat"        "mwlrt.sd_n_pat"         
#> [15] "mwlrt.mean_n_evt"        "mwlrt.sd_n_evt"         
#> [17] "mwlrt.N_missing_n_pat"   "mwlrt.N_missing_n_evt"  
#> 
#> $bootCI
#> character(0)
#> 
#> $extra
#> [1] "REPLICATIONS" "SIM_TIME"     "RAM_USED"     "SEED"         "COMPLETED"   
#> 
#> $errors
#> [1] "ERRORS"
#> 
#> $warnings
#> [1] "WARNINGS"
#> 
# }
```
