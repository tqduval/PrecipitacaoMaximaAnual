
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ----------------------------------------------------------

# Importar dados caso não estejam no 'Global Environment'
path.df.precip.max.anual <- file.choose()
df.precip.max.anual <- readRDS(path.df.precip.max.anual)

# AJUSTE DAS ESTAÇÕES -----------------------------------------------------

# Nível de significância e teste Ratio Likelihood
alfa <- 0.05
lim.l.ratio <- qchisq(alfa, df = 1, lower.tail = FALSE)

# Criar dataframe que resume os dados por e estação e faz o ajuste das distribuições
{ 
  tempo.inicio <- Sys.time()
  
  df.param.completo <-               # corrige séries de precipitação com pd.max = 0
    df.precip.max.anual %>%
    # filter(Estacao_codigo %in% c(47002, 2849002, 3052007)) %>% # filtrar somente 3 estações p/ testar)
    group_by(estacao = Estacao_codigo) %>%
    summarise(
      
      # puxa a série de máximos de cada estação
      pd.max = list(Pdmax %>% as.data.frame()),
      
      # p.0 = probabilidade de determinada série histórica ter pd.max = 0
      p.0 = length(Pdmax[Pdmax == 0])/length(Pdmax),
      
      # tamanho da série
      n.serie = length(pd.max %>% unlist()),
      
      # estatísticas básicas
      estat.basica = list(data.frame( # min, # max
        media = mean(pd.max %>% unlist),
        mediana = median(pd.max %>% unlist),
        desv.pad = sd(pd.max %>% unlist),
        assimetria = moments::skewness(pd.max %>% unlist),
        l.assimetria = lmom::samlmu(x = unlist(pd.max))[[3]],
        curtose = moments::kurtosis(pd.max %>% unlist),
        l.curtose = lmom::samlmu(x = unlist(pd.max))[[4]])), 
      
      # coordenadas
      lat = list(Latitude[1]),
      long = list(Longitude[1]),
      
      # parâmetros iniciais usando LMOM
      param.lmom.gev = list(fun.lmom(xp = Pdmax[Pdmax > 0], dist = "gev")),
      param.lmom.gu = list(fun.lmom(xp = Pdmax[Pdmax > 0], dist = "gumbel")),
      
      # máxima verossimilhança (otimizador)
      max.l.gev = list(fun.max.l(param = unlist(param.lmom.gev), xp = Pdmax[Pdmax > 0], dist = "gev")),  # gev
      max.l.gu = list(fun.max.l(param = unlist(param.lmom.gu), xp = Pdmax[Pdmax > 0], dist = "gumbel")), # gumbel
      max.gl.gev.beta = list(fun.max.gl(param = unlist(param.lmom.gev),                                # gev e priori beta
                                        xp = Pdmax[Pdmax > 0],
                                        priori = "beta",
                                        media.kappa = -0.10,
                                        desvpad.kappa = 0.122)),
      
      # pega o máximo do resultado do optimx
      max.gev = list(max.l.gev %>% as.data.frame() %>% select(value)) %>% unlist() %>% as.numeric(),            # máximo gev
      max.gu = list(max.l.gu %>% as.data.frame() %>% select(value)) %>% unlist() %>% as.numeric(),              # máximo gumbel
      max.gev.beta = list(max.gl.gev.beta %>% as.data.frame() %>% select(value)) %>% unlist() %>% as.numeric(), # máximo gev beta
      
      # hessiana verossimilhança generalizada
      hess.gl = list(max.gl.gev.beta %>% as.data.frame %>% select(hess) %>% unlist %>% as.numeric),
      
      # tirar coluna ´value´ da lista e deixa só os parâmetros resultantes da GML (deixar bonito)
      max.l.gev = list(max.l.gev %>% as.data.frame() %>% select(-value)),
      max.l.gu = list(max.l.gu %>% as.data.frame() %>% select(-value)),
      max.gl.gev.beta = list(max.gl.gev.beta %>% as.data.frame() %>% select(-value)),
      
      # calcular a razão de verossimilhança
      dist.l = ifelse(2*(max.gev - max.gu) > lim.l.ratio, "GEV", "Gumbel"), # entre gumbel likelihood e gev likelihood
      
      # p-value l.gu e l.gev
      p.value.l = pchisq(2*(max.gev - max.gu), df = 1, lower.tail = FALSE),
      
      # Aikake's Information Criterion (AIC) e Bayesian Information Criterion (BIC)
      criterios = list(data.frame(estat = c("AIC", "BIC"),
                                  l.gu = c(-2*max.gu + 2*2, -2*max.gu + 2*log(n.serie)),      # gumbel tem 2 parâmetros
                                  l.gev = c(-2*max.gev + 3*2, -2*max.gev + 3*log(n.serie)),   # gev completa tem 3 parâmetros
                                  gl.gev = c(-2*max.gev.beta + 3*2, -2*max.gev.beta + 3*log(n.serie)))), 
      
      # informações complementares
      bacia.codigo = list(BaciaCodigo[1]),       # pega somente primeiro item da lista
      subbacia.codigo = list(SubBaciaCodigo[1]), # pq são todos iguais
      altitude = list(Altitude[1]),              # info indisponível
      area.drenagem = list(AreaDrenagem[1]),     # info indisponível
      nome.rio = list(NomedoRio[1])              # info indisponível
      
    ) %>%
    unique() # remove linhas duplicadas
    
  cat(paste(df.param[df.param$max.gev.beta == -1e6,] %>% nrow, "estações não foram otimizadas e serão removidas do data.frame final"), "\n")
  # df.param <- df.param %>% filter(max.gev.beta != -1e6)  # remove linhas que não foram otimizadas (49 estações)
  
  tempo.fim <- Sys.time() - tempo.inicio
  cat(paste("Tempo de processamento:", round(tempo.fim, 3), "\n"))
  
  beep(sound = 10)
  
} # final bloco df.param



# Exportar RDS
saveRDS(object = df.param, file = "Dados Gerados/df.param.rds")

# Exportar em CSV
write.table(
  x = df.param %>%
    select(estacao, n.serie, param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta) %>% 
    mutate(csi.lmom.gu = sapply(param.lmom.gu, function(x) x[[1]][[1]]),
           csi.ml.gu = sapply(max.l.gu, function(x) x[[1]][[1]]),
           alpha.lmom.gu = sapply(param.lmom.gu, function(x) x[[2]][[1]]),
           alpha.ml.gu = sapply(max.l.gu, function(x) x[[2]][[1]]),
           csi.lmom.gev = sapply(param.lmom.gev, function(x) x[[1]][[1]]),
           csi.ml.gev = sapply(max.l.gev, function(x) x[[1]][[1]]),
           csi.gml.gev = sapply(max.gl.gev.beta, function(x) x[[1]][[1]]),
           alpha.lmom.gev = sapply(param.lmom.gev, function(x) x[[2]][[1]]),
           alpha.ml.gev = sapply(max.l.gev, function(x) x[[2]][[1]]),
           alpha.gml.gev = sapply(max.gl.gev.beta, function(x) x[[2]][[1]]),
           kappa.lmom.gev = sapply(param.lmom.gev, function(x) x[[3]][[1]]),
           kappa.ml.gev = sapply(max.l.gev, function(x) x[[3]][[1]]),
           kappa.gml.gev = sapply(max.gl.gev.beta, function(x) x[[3]][[1]])) %>% 
    select(-c(param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta)),
  file = "Dados Gerados/df.param.csv",
  sep = ";",
  dec = ",",
  row.names = FALSE,
  fileEncoding = "UTF-8")

# Estações não otimizadas
df.kappa.out <- df.param.completo %>% 
  select(estacao, n.serie, param.lmom.gev, max.gev.beta) %>% 
  mutate(kappa.out = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]])) %>% 
  # filter(!between(kappa.out, -0.5, 0.5))
  filter(max.gev.beta == -1e6)

df.precip.max.anual %>% 
  filter(Estacao_codigo %in% c(2147001, 2854003)) %>% 
  ggplot(aes(x = as_date(Date, format = "%d/%m/%Y"))) +
    facet_wrap(~Estacao_codigo) +
    geom_line(aes(y = Pdmax))