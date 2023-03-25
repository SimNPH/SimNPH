tmp_obj <- SIMDESIGN_TEMPFILE_UUC_D6MYHKW75

tmp1 <- tmp_obj$'1'[1:193]
tmp2 <- tmp_obj$'2'[1:193]
tmp3 <- tmp_obj$'3'[1:193]
tmp4 <- tmp_obj$'4'[1:193]
tmp5 <- tmp_obj$'5'[1:193]
tmp6 <- tmp_obj$'6'[1:193]

result <- tibble::as_tibble(rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6))

View(result)
