
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, moments, patchwork,
               doParallel, foreach, progressr, pbapply)


# IMPORTAR DADOS ----------------------------------------------------------

# Importar dados caso não estejam no 'Global Environment'
path.df.precip.max.anual <- file.choose()
df.precip.max.anual <- readRDS(path.df.precip.max.anual)


# AJUSTE DAS ESTAÇÕES -----------------------------------------------------

# Definir nível de significância e limite do teste de Razão de Verossimilhança
l.ratio.alpha <- 0.05
l.ratio.lim <- qchisq(l.ratio.alpha, df = 1, lower.tail = FALSE)

# Número de estações
n.estacoes <- df.precip.max.anual$Estacao_codigo %>% unique %>% length

# Estações
estacoes <- df.precip.max.anual$Estacao_codigo %>% unique

# Função para processar uma estação
# Esta função está sem as avaliações dos testes de AIC, BIC e LRT
fun.ajuste.pmax <- function(df,                  # data.frame contendo séries AM
                            l.ratio.alpha = 0.05 # nível de significância para o teste da Razão de Verossimilhança
                            ){
  
  tryCatch({
    
    # Informações básicas
    estacao <- unique(df$Estacao_codigo)  # código da estação
    pd.max <- df$Pdmax                    # série AM
    p0 <- sum(pd.max == 0)/length(pd.max) # probabilidade de determinada série ter pd.max = 0
    pd.max.pos <- pd.max[pd.max > 0]      # séries AM sem anos em que pd.max = 0
    n.serie <- length(pd.max)             # tamanho da série
    
    # Ajuste com L-MOM
    par.lmom.gev <- fun.lmom(xp = pd.max.pos, dist = "gev")   # GEV
    par.lmom.gu <- fun.lmom(xp = pd.max.pos, dist = "gumbel") # Gumbel
    
    # Ajuste com ML
    par.l.gev <- fun.max.l(param = unlist(par.lmom.gev),
                           xp = pd.max.pos,
                           dist = "gev")  # GEV
    par.l.gu <- fun.max.l(param = unlist(par.lmom.gu),
                          xp = pd.max.pos,
                          dist = "gumbel") # Gumbel
    
    # Ajuste com GML
    par.gl.gev <- fun.max.gl(param = unlist(par.lmom.gev),
                             xp = pd.max.pos,
                             media.kappa = -0.10,
                             desvpad.kappa = 0.122)
    
    # Extrair máximos
    max.l.gu <- par.l.gu$value
    max.l.gev <- par.l.gev$value
    max.gl.gev <- par.gl.gev$value
    
    # Testes da Razão de Verossimilhança
    l.ratio.lim <- qchisq(l.ratio.alpha, df = 1, lower.tail = FALSE) # limite da distribuição X²
    rv.l <- 2*(max.l.gev - max.l.gu)                                 # GEV_ML e Gumbel_ML
    rv.gl <- 2*(max.gl.gev - max.l.gu)                               # GEV_GML e Gumbel_ML
    l.ratio.l <- tibble(rv = rv.l,
                        p.value = pchisq(q = rv.l, df = 1, lower.tail = FALSE), # p-valor
                        l.ratio = ifelse(rv.l > l.ratio.lim, "GEV", "Gumbel"))  # teste de hipótese
    l.ratio.gl <- tibble(rv = rv.gl,
                         p.value = pchisq(q = rv.gl, df = 1, lower.tail = FALSE), # p-valor
                         l.ratio = ifelse(rv.gl > l.ratio.lim, "GEV", "Gumbel"))  # teste de hipótese
    
    # # Critérios de Informação (AIC e BIC)
    # aic.bic <- tibble(estat = c("AIC", "BIC"),
    #                   l.gu = c(-2*max.l.gu + 2*2, -2*max.l.gu + 2*log(n.serie)),    # gumbel tem 2 parâmetros
    #                   l.gev = c(-2*max.l.gev + 3*2, -2*max.l.gev + 3*log(n.serie)), # gev completa tem 3 parâmetros
    #                   gl.gev = c(-2*max.gl.gev.beta + 3*2, -2*max.gl.gev.beta + 3*log(n.serie)))
    
    # Matriz hessiana
    hess.l.gu <- par.l.gu$hess
    hess.l.gev <- par.l.gev$hess
    hess.gl.gev <- par.gl.gev$hess
    
    # Latitude e Longitude
    lat <- df$Latitude[1]
    long <- df$Longitude[1]
    
    # Bacia e Sub-bacia
    bacia <- df$BaciaCodigo[1]
    sub.bacia <- df$SubBaciaCodigo[1]
    
    # Lista de resultados
    tibble(estacao = estacao,
           pd.max = list(tibble(Pdmax = pd.max)),
           p0 = p0,
           n.serie = n.serie,
           par.lmom.gev = list(tibble(par.lmom.gev)),
           par.lmom.gu = list(tibble(par.lmom.gu)),
           par.l.gev = list(tibble(par.l.gev[,-c(4,5)])),   # tirando coluna 'value' e 'hess'
           par.l.gu = list(tibble(par.l.gu[,-c(3,4)])),     # tirando coluna 'value' e 'hess'
           par.gl.gev = list(tibble(par.gl.gev[,-c(4,5)])), # tirando coluna 'value' e 'hess'
           max.l.gu = max.l.gu,
           max.l.gev = max.l.gev,
           max.gl.gev = max.gl.gev,
           hess.l.gu = hess.l.gu,
           hess.l.gev = hess.l.gev,
           hess.gl.gev = hess.gl.gev,
           l.ratio.l = list(l.ratio.l),
           l.ratio.gl = list(l.ratio.gl),
           # aic.bic = list(aic.bic),
           lat = lat,
           long = long,
           bacia = bacia,
           sub.bacia = sub.bacia,
           .rows = 1)
  
  }, error = function(e){
    message(sprintf("\nErro na estação %s: %s", estacao, e))
    NULL
  })
  
}

# Função de processamento com for loop convencional
fun.process <- function(df, l.ratio.alpha = 0.05){
  
  # Pacotes
  if(!require(pacman)) install.packages("pacman")
  pacman::p_load(pacman, dplyr, beepr)
  
  # Iniciar "cronômetro"
  tempo.inicio <- Sys.time()
  
  # Lista de estações
  estacoes <- unique(df$Estacao_codigo)
  n.estacoes <- length(estacoes)
  
  # Redefinir df.param
  df.param <- tibble()
  
  # Processar estações restantes
  for(i in seq_along(estacoes)){
    
    # Process current station
    df.proxy <-
      df %>% 
      filter(Estacao_codigo == estacoes[i]) %>% 
      fun.ajuste.pmax()
    
    # Juntar resultados
    df.param <- bind_rows(df.param, df.proxy)
    
    cat("\nEstação", estacoes[i], "processada →", round(i/n.estacoes*100, 2), "%")
    
  }
  
  # Mensagem final e duração
  duracao <- round(difftime(Sys.time(), tempo.inicio, units = "mins"), 1)
  message(sprintf("\nProcessamento concluído!\nEstações processadas: %d/%d\nDuração: %s min",
                  nrow(df.param), n.estacoes, duracao))
  beep(sound = 10)
  
  return(df.param)
  
}

# Testar resultados de uma estação
teste <- 
  df.precip.max.anual %>% 
  filter(Estacao_codigo == 47003) %>% 
  fun.ajuste.pmax()

# Rodar processamento completo
df.param <- 
  df.precip.max.anual %>%
  # filter(Estacao_codigo %in% c(57000, 47002, 147010)) %>%
  fun.process()

# Substituir versão anterior
df.param.032525 <- df.param


# AVALIAÇÃO DOS RESULTADOS ------------------------------------------------

# Qtde estações que retornaram 1e6 (7 resultados)
# Estações 2147001, 2854003, 8361005, 2447075, 950001, 968001, 1846002
df.param %>% 
  filter(max.gl.gev == -1e6 | max.l.gu == -1e6 | max.l.gev == -1e6) %>% 
  nrow

# Diferença entre máximos
# Versão de 031125 contém os resultados do optimx usando os métodos c("Nelder-Mead", "BFGS", "L-BFGS-B")
df.diff <- 
  tibble(
    estacao = df.param$estacao,
    dif.max.gl = (df.param$max.gl.gev - df.param.031125$max.gl.gev)/df.param$max.gl.gev)