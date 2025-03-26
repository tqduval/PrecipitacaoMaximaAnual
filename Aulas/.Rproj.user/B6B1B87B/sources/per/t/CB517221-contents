
# IMPORTAR DADOS ----------------------------------------------------------

path <- file.choose()

df.47002 <- read.table(file = path,
                       sep = ",",
                       dec = ".",
                       header = TRUE,
                       fileEncoding = "UTF-8")


# PACOTES -----------------------------------------------------------------

install.packages("pacman")
library(pacman)

pacman::p_load(pacman, tidyverse, ggplot2)


# FUNÇÕES -----------------------------------------------------------------

# Função de verossimilhança - GEV
fun.verossimilhanca.gev <- function(param = c(NA, NA, NA), # parâmetros da GEV: csi, alpha, kappa
                                    x){                    # série de observações pd.max
  
  # Parâmetros iniciais
  csi <- param[1]
  alpha <- param[2]
  kappa <- param[3]
  
  n <- length(x)
  
  # Testar limites da distribuição GEV
  if(kappa < 0){
    
    # Testar limite inferior
    if(min(x) < (csi + alpha/kappa)){
      
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(x - csi)) - (1 - kappa/alpha*(x - csi))^(1/kappa)))
      
    }
    
  } else{ # kappa positivo (Weibull)
    
    # Testar limite superior
    if(max(x) > (csi + alpha/kappa)){
      
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(x - csi)) - (1 - kappa/alpha*(x - csi))^(1/kappa)))
      
    }
  }
  
  # Testar parâmetro de escala (deve ser maior que zero)
  if(alpha < 0) ln.L <- 1e6
  
  return(ln.L)
  
}

# Função de verossimilhança Generalizadaa - GEV
fun.verossimilhanca.generalizada.gev <- function(param = c(NA, NA, NA), # parâmetros da GEV: csi, alpha, kappa
                                                 mu.kappa,              # -0.10
                                                 sd.kappa,              # 0.122
                                                 x){                    # série de observações pd.max
  
  # Parâmetros iniciais
  csi <- param[1]
  alpha <- param[2]
  kappa <- param[3]
  
  n <- length(x)
  
  # Distribuição Beta
  var.kappa <- sd.kappa^2
  a <- 0.5
  b <- mu.kappa + a
  p <- b^2*(1 - b)/var.kappa - b
  q <- p*(1/b - 1)
  
  # Testar limites da distribuição GEV
  if(kappa < 0){
    
    # Testar limite inferior
    if(min(x) < (csi + alpha/kappa)){
      
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      if(abs(kappa) > 0.5){
        ln.L <- 1e6
      } else{
        ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(x - csi)) - (1 - kappa/alpha*(x - csi))^(1/kappa)) +
                    (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa))
      }
      
    }
    
  } else{ # kappa positivo (Weibull)
    
    # Testar limite superior
    if(max(x) > (csi + alpha/kappa)){
      
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      if(abs(kappa) > 0.5){
        ln.L <- 1e6
      } else{
        ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(x - csi)) - (1 - kappa/alpha*(x - csi))^(1/kappa)) +
                    (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa))
      }
      
    }

  }
  
  # Testar parâmetro de escala (deve ser maior que zero)
  if(alpha < 0) ln.L <- 1e6
  
  return(ln.L)
  
}

# Função de maximização da verossimilhança
fun.max.verossimilhanca.gev <- function(param = c(NA, NA, NA),
                                        x){
  
  max.l <-
    optim(par = param,
          fn = fun.verossimilhanca.gev,
          hessian = TRUE,
          x = x) %>% {
            # Organizar resultados finais
            tibble(csi = .$par[1],
                   alpha = .$par[2],
                   kappa = .$par[3],
                   value = -(.$value),
                   hess = list(as_tibble(.$hessian)))}
  
  return(max.l)
  
}

# Função de maximização da verossimilhança generalizada
fun.max.verossimilhanca.generalizada.gev <- function(param = c(NA, NA, NA),
                                                     mu.kappa,
                                                     sd.kappa,
                                                     x){
  
  max.l <-
    optim(par = param,
          fn = fun.verossimilhanca.generalizada.gev,
          hessian = TRUE,
          mu.kappa = mu.kappa,
          sd.kappa = sd.kappa,
          x = x) %>% {
            # Organizar resultados finais
            tibble(csi = .$par[1],
                   alpha = .$par[2],
                   kappa = .$par[3],
                   value = -(.$value),
                   hess = list(as_tibble(.$hessian)))}
  
  return(max.l)
  
}

# AJUSTE ------------------------------------------------------------------

# Vetor com parâmetros iniciais (chute)
param <- c(80, 20, -0.1)

# Teste da função de verossimilhança
fun.verossimilhanca.gev(param = c(80, 20, -0.1),
                        x = df.47002$Pdmax)

# Teste da função de verossimilhança generalizada
fun.verossimilhanca.generalizada.gev(param = c(80, 20, -0.1),
                                     mu.kappa = -0.10,
                                     sd.kappa = 0.122,
                                     x = df.47002$Pdmax)

# Maximização da função de verossimilhança
df.par.ml <- fun.max.verossimilhanca.gev(param = param, x = df.47002$Pdmax)

# Maximização da função de verossimilhança generalizada
df.par.gml <- fun.max.verossimilhanca.generalizada.gev(param = param, mu.kappa = -0.10, sd.kappa = 0.122, x = df.47002$Pdmax)

# Organizar data.frame dos resultados
df.par <- tibble(model = c("ml", "gml"), # tibble() é equivalente à data.frame(), mas mais rígido (e seguro)
                 rbind(df.par.ml,        # juntar os dois resultados
                       df.par.gml)) %>% 
  select(-value)                         # remover a coluna de máximo e manter somente os parâmetros
