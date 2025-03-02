
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

# Configurações do processamento paralelo
parallel.cores <- parallel::detectCores() - 2 # deixar 2 cores livres
parallel.cl <- parallel::makeCluster(parallel.cores)
doParallel::registerDoParallel(parallel.cl)

# Função para processar uma estação
# Esta função está sem as avaliações dos testes de AIC, BIC e LRT
fun.ajuste.pmax <- function(df){
  
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
    par.l.gev <- fun.max.l(param = unlist(par.lmom.gev), xp = pd.max.pos, dist = "gev")  # GEV
    par.l.gu <- fun.max.l(param = unlist(par.lmom.gu), xp = pd.max.pos, dist = "gumbel") # Gumbel
    
    # Ajuste com GML
    par.gl.gev.beta <- fun.max.gl(param = unlist(par.lmom.gev),
                                  xp = pd.max.pos,
                                  priori = "beta",
                                  media.kappa = -0.10,
                                  desvpad.kappa = 0.122)
    
    # Extrair máximos
    max.l.gev <- par.l.gev$value
    max.l.gu <- par.l.gu$value
    max.gl.gev.beta <- par.gl.gev.beta$value
    
    # Matriz hessiana
    hess.gl <- par.gl.gev.beta$hess
    
    # Lista de resultados
    tibble(estacao = estacao,
           pd.max = list(tibble(Pdmax = pd.max)),
           p0 = p0,
           n.serie = n.serie,
           par.lmom.gev = list(tibble(par.lmom.gev)),
           par.lmom.gu = list(tibble(par.lmom.gu)),
           par.l.gev = list(tibble(par.l.gev[,-4])),
           par.l.gu = list(tibble(par.l.gu[,-4])),
           par.gl.gev.beta = list(tibble(par.gl.gev.beta[,-c(4,5)])),
           max.l.gu = max.l.gu,
           max.l.gev = max.l.gev,
           max.gl.gev = max.gl.gev.beta,
           hess.gl.beta = hess.gl,
           .rows = 1) %>% unique
  
  }, error = function(e){
    message(sprintf("Erro na estação %s: %s", estacao, e))
    NULL
  })
  
}

# Função p/ rodar a função fun.ajuste.pmax p/ cada estação em processamento paralelo usando foreach
fun.process <- function(df){
  
  # Iniciar "cronômetro"
  tempo.inicio <- Sys.time()
  
  estacoes <- df %>% 
    group_by(Estacao_codigo) %>% 
    group_split()
  
  total <- length(estacoes)
  
  # Configurar processamento paralelo
  parallel.cores <- parallel::detectCores() - 2 # deixar 2 cores livres
  parallel.cl <- parallel::makeCluster(parallel.cores)
  doParallel::registerDoParallel(parallel.cl)
  
  # Exportar funções críticas
  parallel::clusterExport(cl = parallel.cl, 
                          varlist = c("fun.ajuste.pmax", "fun.lmom", "fun.max.l", "fun.max.gl"),
                          envir = .GlobalEnv)
  
  progressr::handlers(progressr::handler_progress(format = "Processando :current/:total :percent [:bar] :eta"))
  
  # Combina cada estação processada em paralelo usando row_bind em um data.frame usando tibble
  resultados <- progressr::with_progress({
    p <- progressor(total)
    
    foreach(i = 1:total,
            .packages = c("lmom", "optimx", "dplyr", "purrr", "tidyr"),
            .export = c("fun.ajuste.pmax", "fun.ajuste.pmax", "fun.lmom", "fun.max.l", "fun.max.gl"),
            .combine = function(...) {
              res <- suppressWarnings(dplyr::bind_rows(...))
              if(is.null(res) || nrow(res) == 0) {
                return(tibble(estacao = numeric()))
              }
              res
            }) %dopar% {
              tryCatch({
                res <- fun.ajuste.pmax(estacoes[[i]])
                p()
                res
              }, error = function(e) {
                message(sprintf("Erro grave na estação %s: %s", 
                                estacoes[[i]]$Estacao_codigo[1], e$message))
                NULL
              })
            }
  })
  
  # # Filtragem segura
  # if("estacao" %in% names(resultados)) {
  #   resultados <- resultados %>% 
  #     filter(!is.na(estacao), 
  #            map_lgl(pd.max, ~!is.null(.x) && nrow(.x) > 0))
  # } else {
  #   resultados <- tibble()  # Fallback para dataframe vazio
  # }
  
  # Mensagem final e alerta
  duracao <- round(difftime(Sys.time(), tempo.inicio, units = "mins"), 1)
  message(sprintf("\nProcessamento concluído!\nEstações processadas: %d/%d\nDuração: %s minutos",
                  nrow(resultados), total, duracao))
  
  if(nrow(resultados) > 0) {
    beepr::beep(sound = 10)
  } else {
    message("Nenhuma estação processada com sucesso")
  }
  
  stopCluster(parallel.cl)
  
  return(resultados)
  
}

# Função de processamento com for loop convencional
fun.process <- function(df){
  
  # Iniciar "cronômetro"
  tempo.inicio <- Sys.time()
  
  # Lista de estações
  estacoes <- unique(df.precip.max.anual$Estacao_codigo)
  n.estacoes <- length(estacoes)
  
  # Redefinir df.param
  df.param <- tibble()
    
  # Processar primeira estação
  suppressWarnings({
    
    df.param <-
      df.precip.max.anual %>% 
      filter(Estacao_codigo == estacoes[1]) %>% 
      fun.ajuste.pmax()
    
  })
  
  # Processar estações restantes
  for(i in 2:n.estacoes){
    
    estacao <- estacoes[i]
    
    # Processar estações
    suppressWarnings({
      
      df.proxy <-
        df.precip.max.anual %>% 
        filter(Estacao_codigo == estacao) %>% 
        fun.ajuste.pmax()
      
    })
    
    # Juntar resultados
    df.param <- bind_rows(df.param, df.proxy)
    
    cat("\nEstação", estacao, "processada.")
    
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
  filter(Estacao_codigo %in% c(57000)) %>% 
  fun.ajuste.pmax()

# Rodar processamento completo
df.param <- 
  df.precip.max.anual %>%
  fun.process()

# Teste com for loop simples
for(i in 1:n.estacoes){
  
  estacao <- estacoes[i]
  
  suppressWarnings({
    
    df.proxy <-
      df.precip.max.anual %>% 
      filter(Estacao_codigo == estacao) %>% 
      fun.ajuste.pmax()
    
  })
  
  if(i == 1){
    
    df.param <- df.proxy
    
  } else{
    
    df.param <- bind_rows(df.param, df.proxy)
    
    }
  
  cat("\nEstação", estacao, "processada.")
  
}; beep(sound = 10)


# AVALIAÇÃO DOS RESULTADOS ------------------------------------------------

# Qtde estações que retornaram 1e6
df.param %>% 
  filter(max.gl.gev == -1e6) %>% 
  nrow

rowsum(df.param[, c("max.l.gu", "max.l.gev", "max.gl.gev")] == -1e6)
