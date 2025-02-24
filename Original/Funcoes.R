
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# QUANTIL -----------------------------------------------------------------

# GEV
fun.q.gev <- function(p,     # probabilidade de não excedência
                      param  # vetor com parâmetros c(csi, alpha, kappa)
){                 
  
  csi <- param[1]; alpha <- param[2]; kappa <- param[3]
  
  # Quantil
  xp <- csi + alpha/kappa*(1 - (-log(p))^kappa)
  
  return(xp)  # quantil explícito em função dos parâmetros
  # e da probabilidade de não excedência
  
}

# Gumbel
fun.q.gu <- function(p,     # probabilidade de não excedência
                     param  # vetor com parâmetros c(csi, alpha)
){                 
  
  csi <- param[1]; alpha <- param[2]
  
  # Quantil
  xp <- csi - alpha*log(-log(p))
  
  return(xp)  # quantil explícito em função dos parâmetros
  # e da probabilidade de não excedência
  
}


# MOM-L -------------------------------------------------------------------

# Calcula os parâmetros iniciais que serão otimizados pelo
# Estimador de Máxima Verossimilhança (MLE) e
# Máxima Verossimilhança Generalizada (GMLE)

fun.lmom <- function(xp,          # coluna do dataframe c/ máximos anuais
                     dist = "gev" # "gumbel" ou "gev"
){
  
  # Carrega o pacote lmom
  pacman::p_load(lmom, pacman)
  
  # Estimativas dos momentos-L não-tendenciosas da amostra
  est.lmom <- lmom::samlmu(x = xp, sort.data = TRUE)
  
  # testar κ
  if(dist == "gumbel"){
    
    param <- lmom::pelgum(est.lmom)
    param <- data.frame(csi = param[1], alpha = param[2], row.names = NULL)
  }
  else{
    
    param <- lmom::pelgev(est.lmom)
    param <- data.frame(csi = param[1], alpha = param[2], kappa = param[3], row.names = NULL)
    
  }
  
  return(param)
  
}


# ML ----------------------------------------------------------------------

# Estimador de Máxima Verossimilhança - GEV
fun.l.gev <- function(param0, # vetor com os parâmetros, 1. csi 2. alpha, 3. kappa
                      xp){
  
  # Paramentros Iniciais
  
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  n <- length(xp)
  
  # Maxima Verossimilhanca de GEV 
  
  if(kappa < 0){
    if(min(xp) < (csi + alpha/kappa)){
      ln.L <- 1e6
    }else{
      ln.L <- -(- n*log(alpha)
                + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi))
                      - (1 - kappa/alpha*(xp - csi))^(1/kappa)))
    }
  }else{
    if(max(xp) > (csi+alpha/kappa)){
      ln.L <- 1e6
    }else{
      ln.L <- -(- n*log(alpha)
                + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi))
                      - (1 - kappa/alpha*(xp - csi))^(1/kappa)))
    }
  }
  
  if(alpha < 0){
    ln.L <- 1e6
  }
  
  return(ln.L)
  
}

# Estimador de Máxima Verossimilhança - Gumbel
fun.l.gu <- function(param0,
                     xp){
  
  # Parâmetros iniciaiss
  csi <- param0[1]
  alpha <- param0[2]
  
  n <- length(xp)
  
  # Função de verossimilhança
  if (alpha < 0){
    
    ln.L <- 1e6
    
  } else{
    
    y <- (xp - csi)/alpha
    ln.L <- -(-n*log(alpha)
              - sum(y + exp(-y)))
    
  }
  
  return(ln.L) # naturalmente, o sinal do ln é negativo, então no momento q inverte o sinal pra minimizar a função,
               # precisa tirar reverter o sinal do máximo depois
  
}

# GML ---------------------------------------------------------------------

# Estimador de Máxima Verossimilhança Generalizada - GEV
# GMLE - priori Beta
fun.gl.beta.gev <- function(param0,        # vetor com os parâmetros, 1. csi 2. alpha, 3. kappa
                            media.kappa,   # -0,10
                            desvpad.kappa, # 0,122
                            xp){
  
  # Paramentros Iniciais
  
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  # Distribuicao Beta
  
  n <- length(xp)
  
  var.kappa <- desvpad.kappa^2
  
  a <- 0.5
  b <- media.kappa+a
  p <- b^2*(1-b)/var.kappa-b
  q <- p*(1/b-1)
  
  # Maxima Verossimilhanca de GEV 
  
  if(kappa < 0){
    if(min(xp) < (csi + alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha)
                  + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa))
                  + (p - 1)*log(a + kappa)+(q - 1)*log(a - kappa)))
    }
  }else{
    if(max(xp) > (csi+alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha)
                  + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa))
                  + (p - 1)*log(a + kappa)+(q - 1)*log(a - kappa)))
    }
  }
  
  if(abs(kappa) > a){
    ln.GL <- 1e6
  }
  
  if(alpha < 0){
    ln.GL <- 1e6
  }
  
  return(ln.GL) # naturalmente, o sinal do ln é negativo, então no momento q inverte o sinal pra minimizar a função,
                # precisa tirar reverter o sinal do máximo depois
  
}

# GMLE - priori Normal (está com problema, dando NaN)
fun.gl.normal.gev <- function(param0,        # vetor com os parâmetros, 1. csi 2. alpha, 3. kappa
                              media.kappa,   # -0,092
                              desvpad.kappa, # 0,12
                              xp){
  
  # Paramentros Iniciais
  
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  n <- length(xp)
  
  # Maxima Verossimilhanca de GEV 
  
  if(kappa < 0){
    if(min(xp) < (csi + alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha)
                  + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa))
                  + log(1/(desvpad.kappa*sqrt(2*pi))*exp(-0.5*((kappa - media.kappa)/desvpad.kappa)^2))))
    }
  }else{
    if(max(xp) > (csi+alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha)
                  + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa))
                  + log(1/(desvpad.kappa*sqrt(2*pi))*exp(-0.5*((kappa - media.kappa)/desvpad.kappa)^2))))
    }
  }
  
  return(ln.GL) # naturalmente, o sinal do ln é negativo, então no momento q inverte o sinal pra minimizar a função,
                # precisa tirar reverter o sinal do máximo depois
  
}


# PRIORIS -----------------------------------------------------------------

# Distribuições a priori p/ serem plotadas c/ ggplot2
# Priori Beta
fun.priori.beta <- function(x, p, q){
  (0.5 + x)^(p-1)*(0.5 - x)^(q-1)/(beta(p,q))
} # sinais normais

# Priori Normal
fun.priori.normal <- function(x, mu, sd){
  
  1/(sd*sqrt(2*pi))*exp(-1/2*((-x-mu)/sd)^2)
  
} # expressão da fdp normal com sinal de kappa invertido


# OPTIMX ------------------------------------------------------------------

# Máxima verossimilhança 
fun.max.l <- function(param = c(NULL, NULL, NULL),
                      xp,
                      dist = "gev"){
  
  p_load(pacman, optimx, dplyr)
  
  if(dist == "gev"){
    
    
    max.l <- optimx::optimx(param, fun.l.gev, method = c("Nelder-Mead", "BFGS", "L-BFGS-B"), xp = xp) %>%
      slice_min(value) %>%             # pega somente a linha com a menor (maior) verossimilhança
      select(csi, alpha, kappa, value) # pega só os parâmetros e o valor do máximo
    
    
  }
  else {
    if(dist == "gumbel"){
      
      
      max.l <- optimx::optimx(param, fun.l.gu, method = c("Nelder-Mead", "BFGS", "L-BFGS-B"), xp = xp) %>%
        slice_min(value) %>%      # pega somente a linha com a menor (maior) verossimilhança
        select(csi, alpha, value) # pega só os parâmetros e o valor do máximo
      
      
    }
    else{errorCondition(message = "gumbel ou gev")}
  }
  
  max.l$value <- -max.l$value # corrige sinal do máximo
  
  return(max.l)
  
}

# Máxima verossimilhança generalizada
fun.max.gl <- function(param = c(NULL, NULL, NULL),
                       xp,
                       priori = "beta",
                       media.kappa,
                       desvpad.kappa){
  
  # Pacotes
  p_load(pacman, optimx, dplyr)
  
  if(priori == "beta"){
    
    # Otimizador p/ Beta
    
    
    max.optimx <- optimx::optimx(param, 
                                 fun.gl.beta.gev,
                                 method = c("Nelder-Mead", "BFGS", "L-BFGS-B"),
                                 xp = xp,
                                 media.kappa = media.kappa,
                                 desvpad.kappa = desvpad.kappa)
    max.gl <- data.frame(
      csi = max.optimx$csi[which(max.optimx$value == min(max.optimx$value))],
      alpha = max.optimx$alpha[which(max.optimx$value == min(max.optimx$value))],
      kappa = max.optimx$kappa[which(max.optimx$value == min(max.optimx$value))],
      value = min(max.optimx$value),
      hess = nest(as.data.frame(
        attributes(max.optimx)$details[which(max.optimx[, 4] == min(max.optimx[, 4])), "nhatend"][[1]]))) %>% 
      rename("hess" = "data")
    
  }
  else {
    
    if(priori == "normal"){
      
      # Otimizador p/ Normal
      max.gl <- optimx::optimx(param, 
                               fun.gl.normal.gev,
                               method = c("Nelder-Mead", "BFGS", "L-BFGS-B"),
                               xp = xp,
                               media.kappa = media.kappa,
                               desvpad.kappa = desvpad.kappa) %>%
        slice_min(value) %>%             # pega somente a linha com a menor (maior) verossimilhança
        select(csi, alpha, kappa, value) # pega só os parâmetros e o valor do máximo
    }
    else {errorCondition(message = "beta ou normal")}
  }
  
  max.gl$value <- -max.gl$value # corrige sinal do máximo
  
  return(max.gl) # corrige o sinal
  
} # retorna um dataframe com 5 colunas: csi, alpha, kappa, value e hess

# L-RATIO -----------------------------------------------------------------

fun.tao3 <- function(kappa){
  
  tao3 <- 2*(1 - 3^(-kappa))/(1 - 2^(-kappa)) - 3
  
  return(tao3)
  
}

fun.tao4 <- function(kappa){
  
  tao4 <- (1 - 6*2^(-kappa) + 10*3^(-kappa) - 5*4^(-kappa))/(1 - 2^(-kappa))
  
  return(tao4)
  
}

# DERIVADAS PARCIAIS ------------------------------------------------------

# Derivada parcial de yp em relação a csi: ∂yp/∂ξ = 1

# Derivada parcial de yp em relação a alpha
fun.deriv.alpha <- function(p, kappa){
  
  deriv.alpha <- 1/kappa*(1 - (-log(p))^kappa)
  return(deriv.alpha)
  
}

# Derivada parcial do quantil em relação a kappa
fun.deriv.kappa <- function(p, alpha, kappa){
  
  yp <- -log(p)
  deriv.kappa <- - alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp)
  return(deriv.kappa)
  
}

# RESUMO DAS ESTATÍSTICAS BÁSICAS -----------------------------------------

# Função que calcula a tabela resumo igual Papalexiou e Koutsoyiannis (2013)
fun.estat.bas <- function(x){
  c(min = min(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    media = mean(x, na.rm = TRUE),
    assimetria = skewness(x, na.rm = TRUE))}

fun.estat.bas2 <- function(x){
  c(min = min(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE), # adiciona o desvio
    media = mean(x, na.rm = TRUE),
    assimetria = skewness(x, na.rm = TRUE))}

