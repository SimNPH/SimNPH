mod_dir <- "./"
for(i in 1:2500){
  rand_nums <- sample(1:1000000, 3, replace=FALSE)
  new_sim_line <- paste0("$SIMULATION (",rand_nums[1],") (",rand_nums[2]," UNIFORM) (",rand_nums[3]," UNIFORM) ONLYSIM NSUB=1 ;")

  mod_file <- paste0(mod_dir,"run",i,".mod")

  mod_code <- readr::read_lines(mod_file)
  sim_line_num <- grep("\\$SIMULATION",mod_code)

  mod_code[sim_line_num] <- new_sim_line
  readr::write_lines(mod_code,file=mod_file)
}
