
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ----------------------------------------------------------

# Importar dados caso não estejam no 'Global Environment'
path.df.param <- file.choose()
df.param <- readRDS(path.df.param)

# ESTATÍSTICAS BÁSICAS ----------------------------------------------------

# Dataframe somente c/ estatísticas básicas para cada estação
df.estat.bas <- ({
  data.frame(
    estacao = df.param$estacao,
    n.serie = df.param$n.serie,
    min = sapply(df.param$pd.max, min, na.rm = TRUE),
    max = sapply(df.param$pd.max, max, na.rm = TRUE),
    media = sapply(df.param$estat.basica, function(estacao) estacao[[1]][[1]]),
    mediana = sapply(df.param$estat.basica, function(estacao) estacao[[2]][[1]]),
    desv.pad = sapply(df.param$estat.basica, function(estacao) estacao[[3]][[1]]),
    assimetria = sapply(df.param$estat.basica, function(estacao) estacao[[4]][[1]]),
    l.assimetria = sapply(df.param$estat.basica, function(estacao) estacao[[5]][[1]]),
    curtose = sapply(df.param$estat.basica, function(estacao) estacao[[6]][[1]]),
    l.curtose = sapply(df.param$estat.basica, function(estacao) estacao[[7]][[1]])) %>% 
    mutate(intervalo = cut(x = n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55")))})

# Tabela resumo com as estatísticas básicas igual em Papalexiou e Koutsoyiannis (2013)
df.estat.bas_resumo <- df.estat.bas %>% 
  select(-c(intervalo, estacao)) %>% 
  summarise(across(.cols = where(is.numeric), .fns = fun.estat.bas))

# Salvar tabela em CSV
write.table(x = df.estat.bas_resumo, file = "Dados Gerados/df.estat.bas_resumo.csv", sep = ";", dec = ",", fileEncoding = "UTF-8")

# AIC, BIC, RV ---------------------------------------------------------------

# Dataframe somente c/ resultado das avaliações de melhor ajuste
df.aic.bic.rv <- ({
  data.frame(
    estacao = df.param$estacao,
    n.serie = df.param$n.serie,
    l.ratio = df.param$dist.l,
    aic.gu = sapply(df.param$criterios, function(estacao) estacao[[2]][[1]]),
    aic.gev.l = sapply(df.param$criterios, function(estacao) estacao[[3]][[1]]),
    bic.gu = sapply(df.param$criterios, function(estacao) estacao[[2]][[2]]),
    bic.gev = sapply(df.param$criterios, function(estacao) estacao[[3]][[2]])) %>% 
    mutate(    aic = ifelse(aic.gu < aic.gev, "Gumbel", "GEV"), # AIC pega mais Gumbel
               bic = ifelse(bic.gu < bic.gev, "Gumbel", "GEV"), # BIC pega GEV quase sempre
               comp.aic.l.ratio = ifelse(aic == l.ratio, "Converge", "Diverge"),
               comp.bic.l.ratio = ifelse(bic == l.ratio, "Converge", "Diverge"),
               comp.aic.bic = ifelse(aic == bic, "Converge", "Diverge"),
               intervalo = cut(x = n.serie,           # cria uma nova coluna separando os dados
                               breaks = c(15, 25, 35, 45, 55, Inf), # nos intervalos definiddos
                               right = FALSE,
                               labels = c("15-25", "25-35", "35-45", "45-55", "≥55")))
})

# Dataframe c/ distribuição mais frequente c/ melhor ajuste entre RV, AIC e BIC
df.ajuste <- ({
  data.frame(
    dist = c("GEV", "Gumbel"),
    rv = c(df.aic.bic.rv %>% filter(l.ratio == "GEV") %>% nrow(),      # RV GEV: 563
           df.aic.bic.rv %>% filter(l.ratio == "Gumbel") %>% nrow()),  # RV Gumbel: 3237
    aic = c(df.aic.bic.rv %>% filter(aic == "GEV") %>% nrow(),         # AIC GEV: 1054
            df.aic.bic.rv %>% filter(aic == "Gumbel") %>% nrow()),     # AIC Gumbel: 2746
    bic = c(df.aic.bic.rv %>% filter(bic == "GEV") %>% nrow(),         # BIC GEV: 3748
            df.aic.bic.rv %>% filter(bic == "Gumbel") %>% nrow())) %>% # BIC Gumbel: 32
    pivot_longer(cols = c(rv, aic, bic),
                 names_to = "teste",
                 values_to = "valor") %>% 
    mutate(teste = toupper(teste))
})

# Dataframe c/ distribuição mais frequente por intervalo
df.ajuste_intervalo <- ({
  df.aic.bic.rv %>%
  group_by(intervalo) %>%
  summarise(
    rv.gev = sum(l.ratio == "GEV"),
    rv.gumbel = sum(l.ratio == "Gumbel"),
    aic.gev = sum(aic == "GEV"),
    aic.gumbel = sum(aic == "Gumbel"),
    bic.gev = sum(bic == "GEV"),
    bic.gumbel = sum(bic == "Gumbel")
  ) %>%
  ungroup()
  
df.ajuste <- df.aic.bic.rv %>%
  group_by(intervalo) %>%
  summarise(rv.gev = sum(l.ratio == "GEV"),
            rv.gumbel = sum(l.ratio == "Gumbel"),
            aic.gev = sum(aic == "GEV"),
            aic.gumbel = sum(aic == "Gumbel"),
            bic.gev = sum(bic == "GEV"),
            bic.gumbel = sum(bic == "Gumbel")) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("rv") | starts_with("aic") | starts_with("bic"), 
               names_to = c("teste", "dist"), 
               names_pattern = "(.*)\\.(.*)", 
               values_to = "valor") %>%
  mutate(teste = toupper(teste)) %>%
  select(dist, intervalo, teste, valor)})

# Dataframe c/ média dos parâmetros obtidos por cada método p/ cada distribuição por intervalo
df.param_intervalo <- ({
  df.param %>% select(estacao, n.serie, param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta) %>% 
  mutate(csi.lmom.gu = sapply(param.lmom.gu, function(estacao) estacao[[1]][[1]]),
         csi.ml.gu = sapply(max.l.gu, function(estacao) estacao[[1]][[1]]),
         alpha.lmom.gu = sapply(param.lmom.gu, function(estacao) estacao[[2]][[1]]),
         alpha.ml.gu = sapply(max.l.gu, function(estacao) estacao[[2]][[1]]),
         csi.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[1]][[1]]),
         csi.ml.gev = sapply(max.l.gev, function(estacao) estacao[[1]][[1]]),
         csi.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[1]][[1]]),
         alpha.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[2]][[1]]),
         alpha.ml.gev = sapply(max.l.gev, function(estacao) estacao[[2]][[1]]),
         alpha.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[2]][[1]]),
         kappa.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]]),
         kappa.ml.gev = sapply(max.l.gev, function(estacao) estacao[[3]][[1]]),
         kappa.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[3]][[1]]),
         intervalo = cut(x = n.serie,           # cria uma nova coluna separando os dados
                         breaks = c(15, 25, 35, 45, 55, Inf), # nos intervalos definiddos
                         right = FALSE,
                         labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    select(-c(param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta)) %>% 
    group_by(intervalo) %>% 
    summarise(csi.lmom.gu = mean(csi.lmom.gu),
              csi.ml.gu = mean(csi.ml.gu),
              alpha.lmom.gu = mean(alpha.lmom.gu),
              alpha.ml.gu = mean(alpha.ml.gu),
              csi.lmom.gev = mean(csi.lmom.gev),
              csi.ml.gev = mean(csi.ml.gev),
              csi.gml.gev = mean(csi.gml.gev),
              alpha.lmom.gev = mean(alpha.lmom.gev),
              alpha.ml.gev = mean(alpha.ml.gev),
              alpha.gml.gev = mean(alpha.gml.gev),
              kappa.lmom.gev = mean(kappa.lmom.gev),
              kappa.ml.gev = mean(kappa.ml.gev),
              kappa.gml.gev = mean(kappa.gml.gev)) %>% 
    pivot_longer(cols = c(csi.lmom.gu, csi.ml.gu, csi.lmom.gev, csi.ml.gev, csi.gml.gev, 
                          alpha.lmom.gu, alpha.ml.gu, alpha.lmom.gev, alpha.ml.gev, alpha.gml.gev, 
                          kappa.lmom.gev, kappa.ml.gev, kappa.gml.gev),
                 names_to = c("parametro", "metodo", "distribuicao"),
                 names_sep = "\\.") %>%
    pivot_wider(names_from = c(parametro, distribuicao),
                values_from = value) %>%
    mutate(across(c(csi_gu, alpha_gu), ~ if_else(metodo == "gml", NA_real_, .))) %>%
    mutate(metodo = factor(metodo, levels = c("lmom", "ml", "gml"))) %>%
    arrange(intervalo, metodo)
})

write.table(x = df.param_intervalo,
            file = "Dados Gerados/df.param_intervalo.csv",
            dec = ",",
            sep = ";",
            row.names = FALSE,
            fileEncoding = "UTF-8")

# Dataframe c/ média dos parâmetros obtidos por cada método p/ cada dsitribui por RH
df.param_rh <- ({
  df.param %>% select(estacao, bacia.codigo, p.0, param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta) %>% 
    mutate(bacia.codigo = as.integer(bacia.codigo),
           csi.lmom.gu = sapply(param.lmom.gu, function(estacao) estacao[[1]][[1]]),
           csi.ml.gu = sapply(max.l.gu, function(estacao) estacao[[1]][[1]]),
           alpha.lmom.gu = sapply(param.lmom.gu, function(estacao) estacao[[2]][[1]]),
           alpha.ml.gu = sapply(max.l.gu, function(estacao) estacao[[2]][[1]]),
           csi.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[1]][[1]]),
           csi.ml.gev = sapply(max.l.gev, function(estacao) estacao[[1]][[1]]),
           csi.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[1]][[1]]),
           alpha.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[2]][[1]]),
           alpha.ml.gev = sapply(max.l.gev, function(estacao) estacao[[2]][[1]]),
           alpha.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[2]][[1]]),
           kappa.lmom.gev = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]]),
           kappa.ml.gev = sapply(max.l.gev, function(estacao) estacao[[3]][[1]]),
           kappa.gml.gev = sapply(max.gl.gev.beta, function(estacao) estacao[[3]][[1]])) %>% 
    select(-c(param.lmom.gu, max.l.gu, param.lmom.gev, max.l.gev, max.gl.gev.beta)) %>% 
    group_by(bacia.codigo) %>% 
    summarise(p.0 = mean(p.0),
              n.estacoes = n(),
              csi.lmom.gu = mean(csi.lmom.gu),
              csi.ml.gu = mean(csi.ml.gu),
              alpha.lmom.gu = mean(alpha.lmom.gu),
              alpha.ml.gu = mean(alpha.ml.gu),
              csi.lmom.gev = mean(csi.lmom.gev),
              csi.ml.gev = mean(csi.ml.gev),
              csi.gml.gev = mean(csi.gml.gev),
              alpha.lmom.gev = mean(alpha.lmom.gev),
              alpha.ml.gev = mean(alpha.ml.gev),
              alpha.gml.gev = mean(alpha.gml.gev),
              kappa.lmom.gev = mean(kappa.lmom.gev),
              kappa.ml.gev = mean(kappa.ml.gev),
              kappa.gml.gev = mean(kappa.gml.gev)) %>% 
    pivot_longer(cols = c(csi.lmom.gu, csi.ml.gu, csi.lmom.gev, csi.ml.gev, csi.gml.gev, 
                          alpha.lmom.gu, alpha.ml.gu, alpha.lmom.gev, alpha.ml.gev, alpha.gml.gev, 
                          kappa.lmom.gev, kappa.ml.gev, kappa.gml.gev),
                 names_to = c("parametro", "metodo", "distribuicao"),
                 names_sep = "\\.") %>%
    pivot_wider(names_from = c(parametro, distribuicao),
                values_from = value) %>%
    mutate(across(c(csi_gu, alpha_gu), ~ if_else(metodo == "gml", NA_real_, .))) %>%
    mutate(metodo = factor(metodo, levels = c("lmom", "ml", "gml"))) %>%
    arrange(bacia.codigo, metodo)
})

write.table(x = df.param_rh,
            file = "Dados Gerados/df.param_rh.csv",
            dec = ",",
            sep = ";",
            row.names = FALSE,
            fileEncoding = "UTF-8")


# GRÁFICOS ----------------------------------------------------------------

# Distribuição de melhor ajuste conforme os AIC, BIC e RV (gráfico de barras)
plot.ajuste <- ({
  ggplot(df.ajuste, aes(x = teste, y = valor, fill = dist)) +
    geom_bar(stat = "identity",
             position = "stack",
             width = 0.4) +
    geom_text(aes(label = valor, family = "serif", fontface = "bold"),
              size = 4,
              position = position_stack(),
              hjust = 1.2) +
    ylim(c(0, 4000)) +
    coord_flip() + # inverter eixos do gráfico (deixar deitado)
    labs(x = "", y = "", fill = "",) +
    scale_fill_manual(values = c("GEV" = "lightblue", "Gumbel" = "steelblue")) +
    # scale_x_discrete(expand = c(0.1, 1)) +
    theme(panel.background = NULL,
          text = element_text(family = "serif", size = 15, color = "black"),
          axis.text = element_text(color = "black"),
          legend.position = "bottom")
}); plot.ajuste # barras horizontais

plot.ajuste <- ({
  ggplot(df.ajuste, aes(x = teste, y = valor, fill = dist)) +
    geom_bar(stat = "identity",
             position = "stack",
             width = 0.4) +
    geom_text(aes(label = valor, family = "serif", fontface = "bold"),
              size = 4,
              position = position_stack(),
              vjust = -0.5) +
    ylim(c(0, 4000)) +
    # coord_flip() + # inverter eixos do gráfico (deixar deitado)
    labs(x = "", y = "", fill = "",) +
    scale_fill_manual(values = c("GEV" = "lightblue", "Gumbel" = "steelblue")) +
    # scale_x_discrete(expand = c(0.1, 1)) +
    theme(panel.background = NULL,
          text = element_text(family = "serif", size = 15, color = "black"),
          axis.text = element_text(color = "black"),
          legend.position = "bottom")
}); plot.ajuste # barras verticais

ggsave(filename = "Plotagens/Análise Exploratória/Ajuste AIC BIC RV.png",
       plot = plot.ajuste,
       height = 8, # 5 p/ horizontal
       width = 5,  # 8 p/ horizontal
       dpi = 300)

ggplot(df.ajuste_intervalo, aes(x = teste, y = valor, fill = dist)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(intervalo ~ ., scales = "free_y") +
  coord_flip() +
  labs(x = "Intervalo de duração da série histórica", y = "Contagem", fill = "Distribuição") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom")


# Histogramas intervalos c(15, 25, 35, 45, 55, >55)
plot.n.serie.intervalo <- ({
  ggplot(data = df.aic.bic.rv, aes(x = intervalo)) +
    geom_bar(color = "steelblue", fill = "lightblue") +
    ylim(c(0, 1050)) +
    geom_text(stat = "count",                                        # adiciona os rótulos em cima das barras
              aes(label = ..count..), 
              vjust = -0.5, family = "serif") +
    labs(x = "Intervalo de Duração *N* (anos)", y = "Nº Estações") +
    theme(panel.background = NULL,
          text = element_text(family = "serif", size = 13.5, color = "black"),
          axis.title.x = ggtext::element_markdown())
}); plot.n.serie.intervalo
