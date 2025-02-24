
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ---------------------------------------------------------------

# Ler dados de precipitação máxima anual gerados pelo Saulo
df.precip.completa <- read.table(file = "Fonte de Dados/df.precip.diaria.txt",
                                 header = TRUE,
                                 sep = "\t",
                                 dec = ".",
                                 fileEncoding = "UTF-8") %>% 
  select(c(Date, Pdmax, Estacao_codigo, BaciaCodigo, SubBaciaCodigo, Altitude, AreaDrenagem, NomedoRio, Latitude, Longitude))

saveRDS(object = df.precip.completa, file = "Dados Gerados/df.precip.completa.RDS")

# Dataframe resumo
df.precip.completa_resumo <- 
  df.precip.completa %>% 
  group_by(estacao = Estacao_codigo) %>% 
  summarise(n.serie = n()) %>% 
  mutate(intervalo = cut(x = n.serie,
                         breaks = c(15, 25, 35, 45, 55, Inf),
                         right = FALSE,
                         labels = c("15-25", "25-35", "35-45", "45-55", "≥55")))

# OUTLIERS ----------------------------------------------------------------

# Identificar estações que apresentem Pdmax > 700 mm
outlier.limiar <- 700

df.outlier <- df.precip.completa %>%
  select(Estacao_codigo, Date, Pdmax) %>% 
  rename(estacao = Estacao_codigo, data = Date, pd.max = Pdmax) %>% 
  filter(estacao %in% unique(df.precip.completa[which(df.precip.completa$Pdmax > outlier.limiar), 3])) %>%
  mutate(data = as_date(data, format = "%d/%m/%Y")) %>% 
  arrange(estacao, data) %>% 
  group_by(estacao) %>%
  mutate(dif.anos = as.numeric(difftime(data, lag(data), units = "days")) / 365.25) %>% # diferença entre anos consecutivos
  ungroup()

# REMOÇÃO DE OUTLIERS -----------------------------------------------------

# Número de estações em df.outlier (40)
df.outlier$estacao %>% unique %>% length

# Número de estações do df.outlier c/ dif.anos > 5 (11)
df.outlier$estacao[df.outlier$dif.anos >= 5] %>% unique %>% length

# Número de estações do df.param c/ max.gml.gev = 1e+6
df.param.completo$estacao[df.param.completo$max.gev.beta == -1e+6] %>% length

# Verificar sobreposição entre Estações c/ pd.max > 700 e dif.anos > 5 e κ fora de ±0.5 (54 estações)
estacoes.removidas <- rbind(
  # estações não otimizadas pelo gmle (max.gev.beta == -1e6)
  data.frame(filtro = 'gmle', estacao = df.param.completo$estacao[df.param.240817$max.gev.beta == -1e+6]),
  # estações com mais de 5 anos consecutivos sem dados
  data.frame(filtro = 'year.nodata', estacao = df.outlier$estacao[which(df.outlier$dif.anos >= 5)]),
  # estações com kappa fora do intervalo da priori beta
  data.frame(filtro = 'kappa.out', estacao = df.quantis.240817$estacao[(df.quantis.240817$kappa.gml >= 0.5 | df.quantis.240817$kappa.gml <= -0.5)]))

estacoes.removidas <- 
  estacoes.removidas %>% 
  arrange(estacao, desc(filtro == "kappa.out")) %>%
  distinct(estacao, .keep_all = TRUE)  # Remove duplicatas, mantendo kappa.out se presente

estacoes.removidas$estacao %>% unique %>% length # 56 estações

# Juntar todos os dados removidos
df.removido <- ({
  df.precip.completa %>% 
    rename("estacao" = "Estacao_codigo", "pd.max" = "Pdmax", "data" = "Date") %>% 
    mutate(data = as_date(data, format = "%d/%m/%Y")) %>% 
    filter(estacao %in% estacoes.removidas$estacao) %>% 
    mutate(filtro = estacoes.removidas$filtro[match(estacao, estacoes.removidas$estacao)]) %>% 
    bind_rows(df.outlier %>% 
                filter(pd.max > 700) %>% 
                mutate(filtro = "pd.max")) %>% 
    select(estacao, data, pd.max, filtro)
})

write.table(x = df.removido,
            file = "Dados Gerados/df.removido.csv",
            sep = ";",
            dec = ",",
            row.names = FALSE,
            fileEncoding = "UTF-8")

saveRDS(object = df.removido, file = "Dados Gerados/df.removido.RDS")

# Remove os outliers conforme df.outlier e mais de 5 anos consecutivos sem dados
df.precip.max.anual <- df.precip.completa %>% 
  filter(Pdmax < outlier.limiar,                                            # remove anos que tem Pdmax > 700 mm
         !Estacao_codigo %in% df.outlier$estacao[df.outlier$dif.anos >= 5]) # remove as estações c/ Pdmax > 700 mm e +5 anos s/ dados

# Exportar RDS
saveRDS(object = df.precip.max.anual, file = "Dados Gerados/df.precip.max.anual.rds")

# GRÁFICOS ----------------------------------------------------------------

# Plotar estações com pdmax > 700 mm
plot.outlier <- ({
  df.outlier %>% 
    mutate(ano = year(data)) %>% 
    group_by(estacao) %>% 
    complete(ano = seq(min(ano), max(ano), by = 1)) %>% 
    mutate(data = make_date(ano)) %>% 
    ungroup() %>% 
    # Gráfico
    ggplot(aes(x = ano, y = pd.max)) +
    geom_line(linewidth = 0.35, color = "magenta") +
    geom_point(shape = 3, size = 0.2, color = "purple4", alpha = 0.6) +
    facet_wrap(~estacao, scale = "free_y", ncol = 8) +
    labs(x = "", y = "Precipitação [mm]") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = NULL,
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(vjust = 0.5, hjust = 1),
          text = element_text(family = "serif", size = 10, color = "black"))
}); plot.outlier

ggsave(filename = "Plotagens/Análise Exploratória/Outliers 700 mm.png", plot = plot.outlier, 
       height = 190 , width = 350, units = "mm",
       dpi = 300)

# Histograma tamanho das séries
plot.n.serie <- ({
  ggplot(data = df.precip.completa_resumo, aes(x = n.serie)) +
    geom_histogram(binwidth = 10,
                   boundary = 112,
                   breaks = c(15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 112), # colocar de 10 em 10 anos
                   color = "steelblue",
                   fill = "lightblue") +
    geom_vline(xintercept = max(df.aic.bic.rv$n.serie),
               color = "orange",
               linewidth = 0.8) +
    labs(y = "Nº estações\n", x = "N [anos]") +
    geom_text(stat = "bin",
              aes(label = after_stat(count)),
              vjust = -0.8,
              breaks = seq(15, 115, by = 10),
              family = "serif") +
    annotate(geom = "text",
             x = 109.8, y = 100, angle = 90,
             label = "Estação mais longa do conjunto (N = 112)",
             hjust = "left",
             size = 4,
             family = "serif") +
    scale_x_continuous(breaks = seq(15, 115, 20)) +
    ylim(c(0,1010)) +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(color = "white", fill = "white"),
          text = element_text(family = "serif", size = 13.5, color = "black"),
          # axis.title.x = ggtext::element_markdown(), # deixa a legenda do eixo x em itálico ou negrito usando *
          axis.text = element_text(color = "black"),
          axis.ticks = NULL)
}); plot.n.serie

ggsave(filename = "Plotagens/Análise Exploratória/Histograma Tamanho Séries.png",
       plot = plot.n.serie,
       height = 6,
       width = 8,
       dpi = 300)


# DIAGRAMA DE MOMENTOS-L --------------------------------------------------

# Dataframe c/ l-ratios amostrais c/ n.serie >= 30 anos
df.diagrama.amostral <- ({
  df.param %>% 
    select(estacao, n.serie, pd.max, param.lmom.gev) %>% 
    mutate(kappa.gev.lmom = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]])) %>%  # extrair κ do df listado
    select(-c(param.lmom.gev)) %>% 
    mutate(origem = "Pontos amostrais",
           tao3 = sapply(pd.max, function(x) {lmom::samlmu(x = unlist(x))[[3]]}),
           tao4 = sapply(pd.max, function(x) {lmom::samlmu(x = unlist(x))[[4]]}))
})

# Dataframe c/ l-ratios sintéticos da GEV p/ plotar uma curva
df.diagrama.gev <- ({
  data.frame(kappa = seq(from = min(df.diagrama.amostral$kappa.gev.lmom),
                         to = max(df.diagrama.amostral$kappa.gev.lmom),
                         by = 0.001)) %>% 
    mutate(origem = "GEV",
           tao3 = sapply(kappa, fun.tao3),
           tao4 = sapply(kappa, fun.tao4))
})

# Dataframe c/ l-ratios sintéticos da Pareto Generalizada
df.diagrama.gp <- ({
  data.frame(origem = "GP",
             tao3 = (1 - seq(from = -0.66, to = 17, by = 0.001)) / (3 + seq(from = -0.66, to = 17, by = 0.001)),
             tao4 = (1 - seq(from = -0.66, to = 17, by = 0.001)) * (2 - seq(from = -0.66, to = 17, by = 0.001))/((3 + seq(from = -0.66, to = 17, by = 0.001)) * (4 + seq(from = -0.66, to = 17, by = 0.001)))) %>% 
    filter(between(tao3, min(df.diagrama.gev$tao3), max(df.diagrama.gev$tao3)))
})

# Dataframe c/ l-ratios sintéticos para LN2 e LN3
df.diagrama.ln23 <- ({ # as estimativas são válidas p/ |tao3| < 0.9, valores da GEV são válidos
  data.frame(origem = "LN2 e LN3",
             tao3 = seq(from = min(df.diagrama.gev$tao3), to = max(df.diagrama.gev$tao3), by = 0.001)) %>%
    mutate(tao4 = 0.12282 + 0.77518*tao3^2 + 0.12279*tao3^2 - 0.13683*tao3^6 + 0.11368*tao3^8)
})

# Dataframe c/ l-ratios sintéticos p/ Pearson Tipo III
df.diagrama.p3 <- ({ # as estimativas são válidas p/ |tao3| < 0.9, valores da GEV são válidos
  data.frame(origem = "P3",
             tao3 = seq(from = min(df.diagrama.gev$tao3), to = max(df.diagrama.gev$tao3), by = 0.001)) %>%
    mutate(tao4 = 0.1224 + 0.30115*tao3^2 + 0.95812*tao3^2 - 0.57488*tao3^6 + 0.19388*tao3^8)
})

# Dataframe c/ l-ratios Gumbel
df.diagrama.gumbel <- ({
  data.frame(origem = "Gumbel",
             tao3 = 0.1699,
             tao4 = 0.1504)
})

# Juntar data.frames
df.diagrama.lmom <- ({
  rbind(df.diagrama.amostral[
  df.diagrama.amostral$n.serie >= 30,
  c(5,6,7)], # estações c/ N ≥ 30 anos
                          df.diagrama.gev[,-1],                                               # GEV
                          df.diagrama.gumbel,                                                 # Gumbel
                          df.diagrama.gp,                                                     # Pareto Generalizada
                          df.diagrama.ln23,                                                   # LN2 e LN3
                          data.frame(origem = "Ponto médio",                                  # adiciona o ponto médio
                                     tao3 = mean(df.diagrama.amostral$tao3),
                                     tao4 = mean(df.diagrama.amostral$tao4)))})

# Plotar Diagrama de Momentos-L
plot.diagrama.lmom <- ({
  ggplot(data = df.diagrama.lmom, aes(x = tao3, y = tao4, shape = origem, color = origem, linetype = origem)) +
    # Pontos amostrais
    geom_point(data = subset(df.diagrama.lmom, origem == "Pontos amostrais"), size = 1.5, alpha = 0.5) +
    # Ponto médio
    geom_point(data = subset(df.diagrama.lmom, origem == "Ponto médio"), size = 2.5) +
    # Curva GEV
    geom_line(data = subset(df.diagrama.lmom, origem == "GEV"), linewidth = 0.8) +
    # Curva GP
    geom_line(data = subset(df.diagrama.lmom, origem == "GP"), linewidth = 0.8) +
    # Curva LN2 e LN3
    geom_line(data = subset(df.diagrama.ln23, origem == "LN2 e LN3"), linewidth = 0.8) +
    # Curva P3
    geom_line(data = subset(df.diagrama.p3, origem == "P3"), linewidth = 1) +
    # Gumbel
    geom_point(data = subset(df.diagrama.lmom, origem == "Gumbel"), size = 2.5) +
    scale_shape_manual(name = "Legenda",
                       values = c("Pontos amostrais" = 16, "Ponto médio" = 15, "Gumbel" = 17, "P3" = NA, "GP" = NA, "GEV" = NA, "LN2 e LN3" = NA)) +
    scale_color_manual(name = "Legenda",
                       values = c("Pontos amostrais" = "grey80", "GEV" = "#201E43", "Gumbel" = "darkorange", "Ponto médio" = "red", "GP" = "darkgreen", "LN2 e LN3" = "blue", "P3" = "purple")) +
    guides(shape = guide_legend(override.aes = list(linetype = c("solid", "dashed", "blank", "dashed", "dotted", "blank","blank")))) +
    labs(x = "L-assimetria", y = "L-curtose", color = "Legenda", shape = "Legenda", linetype = "Legenda") +
    theme_minimal() +
    theme(panel.background = NULL,
          text = element_text(family = "serif", size = 13.5, color = "black"),
          axis.text = element_text(color = "black"),
          legend.title = element_blank(), # tira o título da legenda
          legend.key = element_blank(),   # tira o fundo dos símbolos
          # legend.position = "bottom",
          legend.position = c(0.01, 0.99),  # posição da legenda (canto superior esquerdo)
          legend.justification = c(0, 1), # justifica a legenda no canto superior esquerdo
          legend.direction = "vertical",  # coloca as legendas uma em cima da outra
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # ajusta as margens da legenda
          legend.spacing.y = unit(0.1, 'cm'))                 # ajusta o espaçamento vertical entre os itens 
}); plot.diagrama.lmom


ggsave(filename = "Plotagens/Análise Exploratória/Diagrama de Momentos-L.png",
       plot = plot.diagrama.lmom,
       height = 6,
       width = 8,
       dpi = 300)

# Quantidade de valores de pontos amostrais que estão à esquerda do ponto de Gumbel (correspondendo à Weibull)
df.diagrama.amostral %>% filter(tao3 < df.diagrama.gumbel$tao3, n.serie >= 30) %>% nrow # 1881 estações
df.diagrama.amostral %>% filter(tao3 > df.diagrama.gumbel$tao3, n.serie >= 30) %>% nrow # 1860 estações

# Número de anos com y = 0 (precipitações nulas) 114
df.precip.max.anual$Estacao_codigo[df.precip.max.anual$Pdmax == 0] %>% length /
  df.precip.max.anual %>% nrow * 100

# Número de estações com y = 0
df.precip.max.anual$Estacao_codigo[df.precip.max.anual$Pdmax == 0] %>% unique %>% length # 25

# Número de estações com y = 0 por bacia
df.precip.max.anual %>% 
  select(Estacao_codigo, BaciaCodigo, Pdmax) %>% 
  # filter(Pdmax == 0) %>% 
  group_by(BaciaCodigo) %>% 
  summarise(NumEstacoes = n_distinct(Estacao_codigo), # n_distinct() pega o número distindo de estações
            NumPdmax0 = n())                          # n() pega o número de observações agrupadas em cada bacia