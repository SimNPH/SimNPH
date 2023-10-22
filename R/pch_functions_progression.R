progression_cdf_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
                               0, -hazard_after,  hazard_after,
                               0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })
  }
}

progression_surv_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      1-as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })
  }
}

progression_pdf_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% Q %*% c(0,0,1))
    })
  }
}

progression_haz_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    S <- sapply(v, \(v_){
      1-as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })

    f <- sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% Q %*% c(0,0,1))
    })

    S/f
  }
}

progression_quant_fun <- function(hazard_before, prog_rate, hazard_after){
  F <- progression_cdf_fun(hazard_before, prog_rate, hazard_after)
  target_fun <- function(q, p){
    F(q) - p
  }

  function(v){
    sapply(v, \(v_){
      uniroot(target_fun, interval=c(0,1), p=v_, extendInt = "upX")$root
    })
  }
}
