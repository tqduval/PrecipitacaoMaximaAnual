
# PACOTES -----------------------------------------------------------------

# Muitos não estão sendo utilizados → organizar depois
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


# VERSÃO P/ TODOS OS MODELOS ----------------------------------------------------------

# Turn df.param into long table
# Ainda não está sendo usado p/ nada
df.param_long <- ({
  df.param %>% 
    select(estacao, n.serie, bacia, sub.bacia, starts_with("par"), starts_with("hess")) %>% 
    pivot_longer(
      cols = c(starts_with("par."), starts_with("hess.")),
      names_to = c("type", "model"),
      names_pattern = "(par|hess)\\.(.*)",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = type,
      values_from = value
    )
})

# Testar depois fun.q.process.all olhando pra df.param normal
# a função calcula procura os nomes padronizados das colunas em df.param p/ saber qual modelo usar
fun.q.process.all <- function(df,
                              tempo.retorno,
                              alpha.ic = 0.05,
                              model = c("lmom.gu", "lmom.gev", "l.gu", "l.gev", "gl.gev")
                              ){
  
  # Timer
  start.time <- Sys.time()
  
  # Initial information
  estacoes <- unique(df$estacao)
  n.estacoes <- length(estacoes)
  p <- 1 - 1/tempo.retorno
  n.p <- length(p)
  
  # Start resulting table
  df.q <- tibble()
  
  # Loop through stations
  for(est in estacoes){
    
    cat("\nEstação", est, "em processamento...")
    
    # Set current station
    df.current <- df %>% filter(estacao == est)
    series.length <- df.current$n.serie
    
    # Loop through models to determine which should be calculated
    for(mod in model){
      
      # Extract parameter vector for the current model
      par.col <- paste0("par.", mod)
      
      # Check if the column exists
      if(!par.col %in% names(df.current)){
        message(sprintf("Station %s: column %s not found. Skipping model %s.", est, par.col, mod))
        next
      }
      
      # Unlist parameters if the column exists and extract them
      par <- unlist(df.current[[par.col]])
      n.par <- length(par) # number of parameters for the current distribution
      
      # Select appropriate quantile function and extract parameters
      if(grepl("gu", mod)){ # check for the distribution suffix in the model name
        
        # Extract parameters
        csi <- par[1]
        alpha <- par[2]
        
        # Quantiles
        q <- fun.q.gu(p, param = par)
        
      } else {               # if it's not Gumbel, then use the GEV quantile function
        
        # Extract parameters
        csi <- par[1]
        alpha <- par[2]
        kappa <- par[3]
        
        # Quantiles
        q <- fun.q.gev(p, param = par)
        
      }
      
      # Estimate uncertainty using the Delta Method -> only for Maximum Likelihood models
      # Set NA for lmom models
      if(mod %in% c("l.gu", "l.gev", "gl.gev")){
        
        # Extract hessian
        hess.col <- paste0("hess.", mod)
        if(!hess.col %in% names(df.current)){
          
          message(sprintf("Station %s: hessian column %s not found. Skipping uncertainty estimation for model %s.", est, hess.col, mod))
          
          # Set columns to NA
          sd.q <- rep(NA, n.p)
          ic.l <- rep(NA, n.p)
          ic.u <- rep(NA, n.p)
          
        } else {
          
          # If the hessian column exists, then estimate quantile variance
          # Extract and invert hessian matrix
          hess <- df.current[[hess.col]] %>% unlist %>% matrix(ncol = n.par, nrow = n.par, byrow = TRUE) # extract hessian for current model
          cov.q <- tryCatch(
            solve(hess),
            error = function(e){
              message(sprintf("\nStation %s: hessian is null as the station was not optimized in model %s → %s", est, mod, e$message))
              matrix(NA, nrow = length(par), ncol = length(par))
            })
          
          # Calculate variance via the Delta Method
          var.q <- rep(NA, n.p) # create empty vector
          for(prob in seq_along(p)){
            
            # Structure gradient matrix depending on the distribution
            if(grepl("gev", mod)){
              
              d.csi <- 1                                                                # dq/dξ
              d.alpha <- 1/kappa*(1 - (-log(p[prob]))^kappa)                               # dq/dα
              yp <- -log(p[prob])                                                       # to simplify dq/dκ expression
              d.kappa <- -alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp) # dq/dκ
              
              grad.q <- c(d.csi, d.alpha, d.kappa)
              
            }
            
            if(grepl("gu", mod)){
              
              d.csi <- 1                     # dq/dξ
              d.alpha <- -log(-log(p[prob])) # dq/dα
              
              grad.q <- c(d.csi, d.alpha)
              
            }
            
            # Variance
            var.q[prob] <- t(grad.q) %*% cov.q %*% grad.q
            
          }
          
          # Estimate confidence intervals
          sd.q <- sqrt(var.q)               # quantile standard deviation
          delta.ci <- qnorm(1 - alpha.ic/2) # standard normal quantile for given confidence level
          ic.l <- q - delta.ci*sd.q         # lower limit of CI
          ic.u <- q + delta.ci*sd.q         # upper limit of CI
          
        }
        
      } else {
        
        # For lmom models: uncertainty values are not estimated
        sd.q <- rep(NA, n.p)
        ic.l <- rep(NA, n.p)
        ic.u <- rep(NA, n.p)
        
      }
      
      # Organize results for current station and model
      df.q.aux <- tibble(estacao = rep(est, n.p),
                         n.serie = rep(series.length, n.p),
                         model = mod,
                         tempo.retorno = tempo.retorno,
                         q = q,
                         sd.q = sd.q,
                         ic.l = ic.l,
                         ic.u = ic.u)
      
      # Combine results
      df.q <- rbind(df.q, df.q.aux)
      
    } # end model loop
    
  } # end station loop
  
  # Processing time and message
  procss.time <- round(difftime(Sys.time(), start.time, units = "mins"), 1)
  message(sprintf("\nProcessing complete!\nStations processed: %d\nDuration: %s min", n.estacoes, procss.time))
  beep(sound = 10) # sound alert
  gc()             # clean memory
  
  return(df.q)
  
}

# Tempos de retorno
tempo.retorno <- c(10.0, 20.0, 50.0, 100.0, 200.0, 500.0)

teste.quantil <- 
  df.param %>% 
  # filter(estacao == 57000) %>%
  # filter(estacao %in% c(47002, 47003, 57000)) %>%
  # head(20) %>% 
  fun.q.process.all(tempo.retorno = tempo.retorno,
                    alpha.ic = 0.05)
