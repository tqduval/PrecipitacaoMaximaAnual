
# PACOTES -----------------------------------------------------------------

pacman::p_load(dplyr, matricks)

# FUNÇÃO ------------------------------------------------------------------

# Ler parâmetros "populacionais"
# Input do tamanho da série

fun.serie.sintetica <- function(par = c(NA, NA, NA), # parâmetros populacionais
                                n.serie = 50,        # tamanho da série
                                n.repl = 1000){
  
  # Probabilidades retiradas de uma distribuição Uniforma definida entre 0 e 1
  p <- matricks::runifm(nrow = n.serie, ncol = n.repl, min = 0, max = 1)
  
  # Calcular "série sintética"
  x.sintetico <- fun.q.gev(p = p, param = par)
  
  return(x.sintetico)
  
}


# Teste
param.pop <- c(80, 20, -0.10)
p <- 1 - 1/tempo.retorno
x.sintetico <- fun.serie.sintetica(par = param, n.serie = 50) # série "observada"
xp.pop <- fun.q.gev(p = p, param = param.pop)                 # quantis populacionais
n.serie <- 50
n.repl <- 1000

# Estimar os parâmetros
df.param.obs <- fun.max.verossimilhanca.generalizada.gev(param = c(120, 50, 0.3),
                                                         mu.kappa = -0.10,
                                                         sd.kappa = 0.122,
                                                         x = x.sintetico[,1])
for(j in 2:n.repl){
  
  df.proxy <- fun.max.verossimilhanca.generalizada.gev(param = c(120, 50, 0.3),
                                                              mu.kappa = -0.10,
                                                              sd.kappa = 0.122,
                                                              x = x.sintetico[,j])
  
  df.param.obs <- rbind(df.param.obs, df.proxy)
  
}

