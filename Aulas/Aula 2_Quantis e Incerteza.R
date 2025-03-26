
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, ggplot2)


# IMPORTAR DADOS ----------------------------------------------------------

# Escolher arquivo no explorador
path <- file.choose()

# Importar arquivo
df.47002 <- read.table(file = path,
                       sep = ",",
                       dec = ".",
                       header = TRUE,
                       fileEncoding = "UTF-8")

# QUANTIS -----------------------------------------------------------------

# Tempos de retorno
tempo.retorno <- c(10, 20, 50, 100) # vetor c/ tempos de retorno desejados
p <- 1 - 1/tempo.retorno            # vetor c/ probabilidade de não excedência

# Extrair os parâmetros desejados (filtrar somente resultado do método GML)
par.gml <- 
  df.par %>% 
  filter(model == "gml") %>% 
  select(c(csi, alpha, kappa)) %>% 
  as.vector # retorna uma lista

# Função quantil GEV
fun.q.gev <- function(p,     # probabilidade de não excedência
                      param  # vetor c/ parâmetros c(csi, alpha, kappa)
                      ){                 
  
  csi <- param[1]; alpha <- param[2]; kappa <- param[3]
  
  # Quantil
  xp <- csi + alpha/kappa*(1 - (-log(p))^kappa)
  
  return(xp)  # quantil explícito em função dos parâmetros
  # e da probabilidade de não excedência
  
}

# Calcular os quantis p/ tempos de retorno desejados
q.gml <- fun.q.gev(p, param = c(par.gml$csi, par.gml$alpha, par.gml$kappa))


# INCERTEZA DOS PARÂMETROS ------------------------------------------------

# Extrair hessiana do método GML
hess.gml <- 
  df.par %>% 
  filter(model == "gml") %>% 
  select(hess) %>% 
  unlist %>%                               # transforma em uma matriz 1x9
  matrix(ncol = 3, nrow = 3, byrow = TRUE) # cria novamente a matriz 3x3 linha a linha

# Tira a inversa da hessiana -> Matriz de Covariância
cov.gml <- solve(hess.gml)

# Extrair diagonais principais
# Variância nos parâmetros
var.par <- c(cov.gml[1,1],
             cov.gml[2,2],
             cov.gml[3,3])

# Desvio padrão dos parâmetros
sd.par <- sqrt(var.par)

# Intervalos de confiança dos parâmetros
alpha.ic <- 0.05                            # nível de confiança
delta.ic <- qnorm(p = 1 - alpha.ic/2, 0, 1) # quantil da normal padrão

# Criar colunas no data.frame -> colunas de 1 a 3 (parâmetros)
df.par.gml$ic.l <- df.par.gml[,1:3] - delta.ic*sd.par # limite inferior
df.par.gml$ic.u <- df.par.gml[,1:3] + delta.ic*sd.par # limite superior


# INCERTEZA DOS QUANTIS ---------------------------------------------------

# Função p/ calcular a variância dos quantis
fun.var.q <- function(p, par, hess){
  
  tryCatch({
    
    # Parâmetros
    par <- unlist(par) # os parâmetros no df estão "listed"
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

# Testar aplicação da função
var.q <- fun.var.q(p,                                                           # vetor c/ probabilidades
                   par = c(df.par.gml$csi, df.par.gml$alpha, df.par.gml$kappa), # parâmetros do método GML
                   hess = df.par.gml$hess)                                      # hessiana do método GML

# Organizando os resultados
df.q <- tibble(                                 # cria um df c/ tibble (posteriormente o modelo do df c/ todas as estações terá essa forma)
  estacao = 47002,                              # código da estação
  par = list(tibble(csi = df.par.gml$csi,       # coluna que armazena os parâmetros dentro de um df listado (como no df.param da aula 1)
                    alpha = df.par.gml$alpha,   
                    kappa = df.par.gml$kappa)),
  hess = df.par.gml$hess,                       # hessiana
  tempo.retorno = tempo.retorno,                # tempos de retorno
  prob = p,                                     # probabilidades calculadas
  q = q.gml) %>%                                # vetor c/ quantil calculado
  mutate(var.q = purrr::pmap_dbl(list(prob, par, hess), # calcula as variâncias linha a linha usando pmap
                                 function(prob, par, hess){
                                   fun.var.q(p = prob, par = par, hess = hess)}),
         ic.l = q - delta.ic*sqrt(var.q),
         ic.u = q + delta.ic*sqrt(var.q))


# GRÁFICO DE FREQUÊNCIAS --------------------------------------------------

pacman::p_load(pacman, dplyr, ggplot2, purrr) # pacote p/ plotagens no R

# Calcular probabilidadee empírica
df.47002.obs <- ({
  
  df.47002 %>% 
    select(Estacao_codigo, Pdmax) %>%  # selecionando colunas
    arrange(desc(Pdmax)) %>%           # ordenando na ordem decrescente
    mutate(n = length(Estacao_codigo), # número total de observações
           i = row_number(),           # índice de cada observação ordenada
           p.weibull = i/(n + 1),      # posição de plotagem de Weibull (prob empírica)
           Tr = 1/p.weibull)           # tempo de retorno calculado
  
})

# Calcular quantis
tempo.retorno <- c(1.5, 2, 3, 5, 7, 10, 20, 50, 100, 200, 500) # vetor c/ tempos de retorno p/ plotagem 

# Data.frame de quantis estimados p/ plotagem
df.47002.est <- ({
  
  tibble(
    estacao = 47002,
    par = list(tibble(csi = df.par.gml$csi,
                      alpha = df.par.gml$alpha,
                      kappa = df.par.gml$kappa)),
    hess = df.par.gml$hess,
    tempo.retorno = tempo.retorno,
    prob = 1 - 1/tempo.retorno,
    q = fun.q.gev(p = prob, param = unlist(par))) %>% 
    mutate(var.q = purrr::pmap_dbl(list(prob, par, hess),
                                        function(prob, par, hess){
                                          fun.var.q(p = prob, par = par, hess = hess)}),
           ic.l = q - delta.ic*sqrt(var.q),
           ic.u = q + delta.ic*sqrt(var.q))
  
})

# Gráfico de frequência
ggplot(data = df.47002.est, aes(x = tempo.retorno, y = q)) +
  geom_line()

plot.freq <- ({
  ggplot() +
    # Precipitação estimada (quantis)
    geom_line(data = df.47002.est, aes(x = tempo.retorno, y = q, linetype = "Estimado"), linewidth = 0.6) +
    # Intervalos de confiança
    geom_line(data = df.47002.est, aes(x = tempo.retorno, y = ic.l, linetype = "IC 95%"), linewidth = 0.4) +
    geom_line(data = df.47002.est, aes(x = tempo.retorno, y = ic.u, linetype = "IC 95%"), linewidth = 0.4) +
    # Precipitação observada
    geom_point(data = df.47002.obs, aes(x = Tr, y = Pdmax, shape = "Observado", color = "Observado"), size = 1.5) +
    scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500),
                  labels = c("10", "20", "50", "100", "200", "500")) + 
    scale_linetype_manual(values = c("Estimado" = "solid",
                                     "IC 95%" = "dashed")) +
    scale_color_manual(values = c("Observado" = "red")) +
    scale_shape_manual(values = c("Observado" = 8)) + 
    labs(x = "TR [anos]", y = "Precipitação diária máxima anual [mm]", title = "Análise de frequências",
         linetype = "", shape = "", colour = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"))
})

ggsave(filename = "Plotagens/Análise de Frequência.png",
       plot = plot.freq,
       width = 12,
       height = 10,
       units = "cm",
       dpi = 300)
