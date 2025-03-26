
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, moments, patchwork,
               doParallel, foreach, progressr, pbapply,
               purrr # pacote p/ usar a função pmap no cáuculo das variâncias dos quantis
               )


# IMPORTAR DADOS ----------------------------------------------------------

# path.df.param <- file.choose() # tbl_df/data.frame containing fitted parameters
# df.param <- readRDS(path.df.param)


# CÁLCULO DOS QUANTIS -----------------------------------------------------

# Function for quantile variance -> implemented inside fun.quantile function
fun.var.q <- function(p, par, hess){
  
  tryCatch({
    
    # Parâmetros
    par <- unlist(par)   # parameters column inside df are "listed"
    n.par <- length(par) # number of parameters
    csi <- par[1]
    alpha <- par[2]
    kappa <- par[3]
    
    # Hessiana
    hess <- hess %>% unlist %>% matrix(ncol = 3, nrow = 3, byrow = TRUE)
    
    # Matriz de covariância -> inverte a hessiana
    cov.par <- solve(hess)
    
    # Criar vetor vazio para armazenar as variâncias
    var.q <- rep(NA, length(p))
    
    # Cálculo p/ vetor de probabilidades
    # P/ cada probabilidade (Tr), calcular a matriz gradiente (derivadas parciais) e 
    for(i in 1:length(p)){
      
      # Derivadas parciais do quantil
      d.csi <- 1                                                                # csi
      d.alpha <- 1/kappa*(1 - (-log(p[i]))^kappa)                               # alpha
      yp <- -log(p[i])                                                          # simplificar a expressão
      d.kappa <- -alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp) # kappa
      
      # Vetor gradiente
      grad.q <- c(d.csi, d.alpha, d.kappa) # matriz coluna 3x1
      
      # Variância do quantil
      var.q[i] <- t(grad.q) %*% cov.par %*% grad.q
      
    }
    
  }, error = function(e){
    
    # Usando o tryCatch(), podemos retornar uma mensagem de erro, isso ajuda no futuro a
    # identificar a fonte do problema caso haja algum
    
    message(sprintf("\nErro no cálculo da variância dos quantis:\n", e))
    
  })
  
  return(var.q)
  
}

# Process all stations (only gl.gev model at the moment)
# This function structures a tibble for a single station at a time
fun.quantile <- function(df,               # data.frame with fitted parameters (df.param)
                         tempo.retorno,    # vector of return periods to estimate quantiles
                         alpha.ic = 0.05   # significance level for quantile confidence intervals
                         ){
  
  # Initial information
  estacao <- df$estacao        # gauge codes
  n.estacao <- length(estacao) # number of stations
  p <- 1 - 1/tempo.retorno     # probabillity of non-exceedance -> vector if length(return.period) > 1
  n.p <- length(p)             # lenght of probabilities vector
  
  # GEV-GML parameters
  par <- unlist(df$par.gl.gev) # parameter vector is listed
  n.par <- length(par)         # number of distribution parameters
  csi <- par[1]                # extract GEV location parameter
  alpha <- par[2]              # extract GEV scale parameter
  kappa <- par[3]              # extract GEV shape parameter
  
  # Quantiles
  q <- fun.q.gev(p, param = par)
  
  # Hessian matrix extraction and uncertainty estimation by the Delta Method
  tryCatch({
    
    # Extract hessian and turn it back into square matrix
    hess <- df$hess.gl.gev %>% unlist %>% matrix(ncol = n.par, nrow = n.par, byrow = TRUE)
    
    # Inverse hessian matrix (covariance matrix) -> if error, variance will be NA
    cov.par <- tryCatch(solve(hess), error = function(e){
      
      message(sprintf("\nUnable to inverse Hessian for station %s: %s", estacao, e$message))
      matrix(NA, nrow = n.par, ncol = n.par)
        
      })
    
    var.q <- rep(NA, n.p)  # create empty variance matrix
    
    # Calculate, for each p, partials, gradient matrix, and quantile variance vector
    for(i in 1:n.p){
      
      # Quantile partials in relation to each GEV parameter
      d.csi <- 1                                                                # dq/dξ
      d.alpha <- 1/kappa*(1 - (-log(p[i]))^kappa)                               # dq/dα
      yp <- -log(p[i])                                                          # simplifying dq/dκ expresison
      d.kappa <- -alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp) # dq/dκ
      
      # Vetor gradiente
      grad.q <- c(d.csi, d.alpha, d.kappa) # matriz coluna 3x1
      
      # Variância do quantil
      var.q[i] <- t(grad.q) %*% cov.par %*% grad.q
      
    }
    
  }, error = function(e){
    
    message(sprintf("\nStation %s: Hessian inversion failed - %s", estacao, e$message))
    var.q <- NA
    
  })
  
  # Confidence intervals - CI
  delta.ci <- qnorm(p = 1 - alpha.ic/2, 0, 1) # standard normal quantile
  sd.q <- sqrt(var.q)                         # quantile standard deviation
  ic.l <- q - delta.ci*sd.q                   # CI's lower limit
  ic.u <- q + delta.ci*sd.q                   # CI's upper limit
  
  # Resulting tibble
  res <- tibble(estacao = rep(estacao, n.p), # repeat gauge number for n.p rows
                tr = tempo.retorno,
                hess = list(tibble(hess)),
                q = q,
                sd.q = sd.q,
                ic.l = ic.l,
                ic.u = ic.u)
  
  return(res)
  
}

# This funcion loops through all stations running fun.quantile,
# binding each stations tibble together to strucure the resulting tibble
fun.quantile.process <- function(df,
                                 tempo.retorno,
                                 alpha.ic = 0.05){
  
  # Load packages
  pacman::p_load(pacman, dplyr, beepr)
  
  # Start timer
  start.time <- Sys.time()

  # List gauges
  estacoes <- df$estacao          # vector of all stations
  n.estacoes <- length(estacoes)  # number of stations
  
  tryCatch({
    
    # Start resulting tibble
    df.q <- tibble()
    
    # Run fun.quantile for all gauges to estimate quantiles and confidence intervals
    for(i in seq_along(estacoes)){

      # Select only the current station
      df.estacao <- df %>% filter(estacao == estacoes[i])
      
      # Run fun.quantile for current station
      df.proxy <- fun.quantile(df = df.estacao, tempo.retorno, alpha.ic = alpha.ic)

      # Combine with previous processed stations
      df.q <- bind_rows(df.q, df.proxy)

      cat("\nEstação", estacoes[i], "processada →", round(i/n.estacoes*100, 2), "%")

    }
    
  }, error = function(e){
    
    message(e$message)
    
  })
  
  # Processing time message
  procss.time <- round(difftime(Sys.time(), start.time, units = "mins"), 1)
  message(sprintf("\nProcessamento concluído!\n%d quantis calculados de %d\nDuração: %s min",
                  nrow(df.q), n.estacoes*length(tempo.retorno), procss.time))
  beep(sound = 10)
  
  return(df.q)
  
}

tempo.retorno <- c(10.0, 20.0, 50.0, 100.0, 200.0, 500.0)

# Testing fun.quantile for a single station
teste.quantil <- 
  df.param %>% 
  filter(estacao == 2147001) %>% 
  fun.quantile(tempo.retorno = tempo.retorno, alpha.ic = 0.05)

# Aplicação da função fun.quantile.process
df.quantiles <- 
  df.param %>% 
  # filter(estacao %in% c(57000, 47002, 2147001, 147010)) %>% 
  fun.quantile.process(tempo.retorno = tempo.retorno, alpha.ic = 0.05)