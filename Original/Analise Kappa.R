
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ----------------------------------------------------------

# Importar dados caso não estejam no 'Global Environment'
path.df.quantis <- file.choose()
df.quantis <- readRDS(path.df.quantis)

# Importar SF Brasil
path.sf.brasil <- file.choose()
sf.brasil <- st_read(dsn = path.sf.brasil)

# GERAR DADOS MAPA --------------------------------------------------------

# Dataframe relacionaal Regiões Hidrográficas DNAEE e nomes
df.rh <- ({
  data.frame(id = seq(1, 8),
             rh = c("Amazonas", "Tocantins", "Atlântico\nNorte/Nordeste", "São Francisco",
                    "Atlântico\nLeste", "Paraná", "Uruguai", "Atlântico\nSudeste"))})

# Cria um objeto espacial contento as estações, tamanho das séries e os valores de kappa ajustados por gml e lmom
df.kappa.mapa <- ({
  df.param %>% 
    select(estacao, n.serie, bacia.codigo, param.lmom.gev, max.gl.gev.beta, lat, long) %>% 
    mutate(kappa.lmom = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]]),
           kappa.gml = sapply(max.gl.gev.beta, function(estacao) estacao[[3]][[1]]),
           bacia.codigo = as.integer(bacia.codigo),
           intervalo = cut(x = kappa.gml,
                           breaks = c(-Inf ,-0.2, -0.1, 0, 0.1, Inf), # intervalos de kappa
                           right = FALSE,                             # se FALSE, inclui o limite inferior, mas não o superior
                           labels = c("1", "2", "3", "4", "5"))) %>% 
    select(-c(param.lmom.gev, max.gl.gev.beta))
})

# Transforma em um objeto vetorial SF
sf.kappa <- sf::st_as_sf(x = df.kappa.mapa, coords = c("long", "lat"), crs = 4326) %>%  # WGS84
  st_transform(st_crs(sf.brasil))
st_write(obj = sf.kappa, dsn = "Dados Gerados/sf.kappa.gpkg")

# RESUMO KAPPA ------------------------------------------------------------

# Tabela resumo ajuste de kappa
df.kappa.gml_resumo2 <- ({
  df.param %>% 
    select(estacao, n.serie, bacia.codigo, max.gl.gev.beta) %>% 
    mutate(kappa.gml = sapply(max.gl.gev.beta, function(estacao) estacao[[3]][[1]]),
           intervalo = cut(n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    select(-max.gl.gev.beta) %>% 
    group_by(intervalo) %>%
    summarise(estatisticas = list(fun.estat.bas2(kappa.gml))) %>%
    unnest_wider(estatisticas) %>% # expande a lista de estatísticas em colunas separadas
    pivot_longer(cols = -intervalo, names_to = "estatistica", values_to = "valor") %>%
    pivot_wider(names_from = intervalo, values_from = valor)
})

write.table(x = df.kappa.gml_resumo,
            file = "Dados Gerados/df.kappa.gml_resumo.csv",
            sep = ";",
            dec = ",",
            row.names = FALSE,
            fileEncoding = "UTF-8")

# Conferência do resumo
df.param %>% 
  select(estacao, n.serie, max.gl.gev.beta) %>% 
  mutate(kappa.gml = sapply(max.gl.gev.beta, function(estacao) estacao[[3]][[1]]),
         intervalo = cut(n.serie,
                         breaks = c(15, 25, 35, 45, 55, Inf),
                         right = FALSE,
                         labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
  filter(intervalo == "15-25") %>% 
  group_by(intervalo) %>% 
  summarise(sd = sd(kappa.gml))

# DISTRIBUIÇÃO DE KAPPA -----------------------------------------------------

# Histograma de κ GMLE p/ N > 40 anos
plot.hist.kappa.gml <- ({
  
  # N > 40
    ggplot(data = df.quantis %>% 
             filter(kappa.gml >= -0.5 & kappa.gml <= 0.5,
                    n.serie >= 40),
           aes(x = kappa.gml)) +
    
    # Histograma c/ densidade
    geom_histogram(aes(y = after_stat(density)),
                   color = "#c15123", fill = "orange",
                   bins = 40) +
    
    # Distribuição normal teórica dos valores de kappa
    stat_function(fun = dnorm, n = 1000, args = list(
      mean = mean(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40]),
      sd = sd(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40])),
      linewidth = 0.7,
      aes(color = "Normal teórica ajustada", linetype = "Normal teórica ajustada")) +
    
    annotate("text", x = 0.5, y = 5, hjust = 1, vjust = 1, family = "serif", size = 5,
             label = paste0("N(", round(mean(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40]), 3),
                            ", ", round(sd(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40]), 3), "²)")) +
    
    # Priori Beta
    stat_function(fun = fun.priori.beta, n = 1000,
                  args = list(p = 6, q = 9), linewidth = 0.7,
                  aes(color = "Beta (Martins e Stedinger, 2000)", linetype = "Beta (Martins e Stedinger, 2000)")) +
    
    # Priori Normal
    stat_function(fun = fun.priori.normal, n = 1000,
                  args = list(mu = 0.092, sd = 0.12), linewidth = 0.7,
                  aes(color = "Normal (Papalexiou e Koutsoyiannis, 2013)", linetype = "Normal (Papalexiou e Koutsoyiannis, 2013)")) +
    
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(name = "Legenda",
                       values = c("Normal teórica ajustada" = "black",
                                  "Beta (Martins e Stedinger, 2000)" = "purple",
                                  "Normal (Papalexiou e Koutsoyiannis, 2013)" = "blue")) +
    scale_linetype_manual(name = "Legenda",
                          values = c("Normal teórica ajustada" = "solid",
                                     "Beta (Martins e Stedinger, 2000)" = "dashed",
                                     "Normal (Papalexiou e Koutsoyiannis, 2013)" = "twodash")) +
    scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, 5.1)) +
    scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1), limits = c(-0.5, 0.5)) +
    labs(x = "Parâmetro de forma κ da GEV", y = "Densidade", color = "", linetype = "") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white", color = "white"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(family = "serif",
                                     size = 14,
                                     color = "black"),
          text = element_text(family = "serif",
                              size = 13,
                              color = "black"))
}); plot.hist.kappa.gml

# Histograma de κ MOM-L p/ N > 40 anos
plot.hist.kappa.lmom <- ({
  
  # N > 40
  ggplot(data = df.quantis %>% 
           filter(kappa.lmom >= -0.5 & kappa.lmom <= 0.5,
                  n.serie >= 40),
         aes(x = kappa.lmom)) +
    
    # Histograma c/ densidade
    geom_histogram(aes(y = after_stat(density)),
                   color = "#c15123", fill = "orange",
                   bins = 40) +
    
    # Distribuição normal teórica dos valores de kappa (MOM-L)
    stat_function(fun = dnorm, n = 1000, args = list(
      mean = mean(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40]),
      sd = sd(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40])),
      linewidth = 0.7,
      aes(color = "Normal teórica ajustada (MOM-L)", linetype = "Normal teórica ajustada (MOM-L)")) +
    
    # Distribuição normal teórica dos valores de kappa (GMLE)
    stat_function(fun = dnorm, n = 1000, args = list(
      mean = mean(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40]),
      sd = sd(df.quantis$kappa.gml[df.quantis$kappa.gml >= -0.5 & df.quantis$kappa.gml <= 0.5 & df.quantis$n.serie >= 40])),
      linewidth = 0.7,
      aes(color = "Normal teórica ajustada (GMLE)", linetype = "Normal teórica ajustada (GMLE)")) +
    
    annotate("text", x = 0.5, y = 5, hjust = 1, vjust = 1, family = "serif", size = 5,
             label = paste0("N(", round(mean(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40]), 3),
                            ", ", round(sd(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40]), 3), "²)")) +
    
    # Priori Beta
    stat_function(fun = fun.priori.beta, n = 1000,
                  args = list(p = 6, q = 9), linewidth = 0.7,
                  aes(color = "Beta (Martins e Stedinger, 2000)", linetype = "Beta (Martins e Stedinger, 2000)")) +
    
    # Priori Normal
    stat_function(fun = fun.priori.normal, n = 1000,
                  args = list(mu = 0.092, sd = 0.12), linewidth = 0.7,
                  aes(color = "Normal (Papalexiou e Koutsoyiannis, 2013)", linetype = "Normal (Papalexiou e Koutsoyiannis, 2013)")) +
    
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(name = "Legenda",
                       values = c("Normal teórica ajustada (MOM-L)" = "black",
                                  "Normal teórica ajustada (GMLE)" = "red",
                                  "Beta (Martins e Stedinger, 2000)" = "purple",
                                  "Normal (Papalexiou e Koutsoyiannis, 2013)" = "blue")) +
    scale_linetype_manual(name = "Legenda",
                          values = c("Normal teórica ajustada (MOM-L)" = "solid",
                                     "Normal teórica ajustada (GMLE)" = "solid",
                                     "Beta (Martins e Stedinger, 2000)" = "dashed",
                                     "Normal (Papalexiou e Koutsoyiannis, 2013)" = "twodash")) +
    scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, 5.1)) +
    scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1), limits = c(-0.5, 0.5)) +
    labs(x = "Parâmetro de forma κ da GEV", y = "Densidade", color = "", linetype = "") +
    guides(linetype = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white", color = "white"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(family = "serif",
                                     size = 14,
                                     color = "black"),
          text = element_text(family = "serif",
                              size = 13,
                              color = "black"))
}); plot.hist.kappa.lmom

ggsave(filename = "Plotagens/Ajuste das Estações/Histograma Kappa Tamanho da Série 40 anos (MOM-L).png",
       plot = plot.hist.kappa.lmom,
       height = 15, width = 25, units = "cm",
       dpi = 300)

# Gráfico de barras κ médio por intervalo de duração das séries (GMLE)
plot.kappa.n.serie.gml <- ({ # usando kappa ajustado por GML
  df.quantis %>% 
    mutate(intervalo = cut(x = n.serie,           # cria uma nova coluna separando os dados
                           breaks = c(15, 25, 35, 45, 55, Inf), # nos intervalos definiddos
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    group_by(intervalo) %>% 
    summarise(media.kappa = mean(kappa.gml),
              prct.kappa.neg = sum(kappa.gml < 0)/n()*100) %>%
    ggplot(aes(x = intervalo, y = media.kappa)) +
    geom_bar(stat = "identity", color = "#c15123", fill = "#FFD550") +
    geom_text(aes(label = format(media.kappa, digits = 2), family = "serif"), # rótulo do valor médio de κ por intervalo
              size = 3.5,
              position = position_stack(),
              vjust = -0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prct.kappa.neg),
                  family = "serif", fontface = "bold"), # rótulo da % de valores negativos
              size = 3.5,
              position = position_stack(vjust = 0.5)) +
    scale_y_continuous(limits = c(0, -0.08), trans = "reverse") +
    # scale_y_continuous(limits = c(0.04, -0.03), trans = "reverse") + # usar esse aqui pra plotar kappa.lmom
    labs(y = "Valor médio de κ", x = "N [anos]", title = "a)") +
    theme_minimal() +
    theme(panel.background = NULL,
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.kappa.n.serie.gml

# Gráfico de barras κ médio por intervalo de duração das séries (MOM-L)
plot.kappa.n.serie.lmom <- ({ # usando kappa ajustado por GML
  df.quantis %>% 
    mutate(intervalo = cut(x = n.serie,           # cria uma nova coluna separando os dados
                           breaks = c(15, 25, 35, 45, 55, Inf), # nos intervalos definiddos
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    group_by(intervalo) %>% 
    summarise(media.kappa = mean(kappa.lmom),
              prct.kappa.neg = sum(kappa.lmom < 0)/n()*100) %>%
    ggplot(aes(x = intervalo, y = media.kappa)) +
    geom_bar(stat = "identity", color = "#c15123", fill = "#FFD550") +
    geom_text(aes(label = format(media.kappa, digits = 2), family = "serif"), # rótulo do valor médio de κ por intervalo
              size = 3.5,
              position = position_stack(),
              vjust = -0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prct.kappa.neg),
                  family = "serif", fontface = "bold"), # rótulo da % de valores negativos
              size = 3.5,
              position = position_stack(vjust = 0.5)) +
    scale_y_continuous(limits = c(0.03, -0.025), trans = "reverse") +
    # scale_y_continuous(limits = c(0.04, -0.03), trans = "reverse") + # usar esse aqui pra plotar kappa.lmom
    labs(y = "Valor médio de κ", x = "N [anos]", title = "a)") +
    theme_minimal() +
    theme(panel.background = NULL,
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.kappa.n.serie.lmom

ggsave(filename = "Plotagens/Ajuste das Estações/Kappa médio por intervalo de tamanho das séries.png",
       plot = plot.kappa.n.serie,
       height = 11, width = 9.5, units = "cm",
       dpi = 300)

# MAPA ESTAÇÕES --------------------------------------------------------------

plot.estacoes.br <- ({
  ggplot() +
    geom_sf(data = sf.brasil, fill = "white", color = "black", linewidth = 0.3) + # Brasil
    geom_sf(data = sf.kappa, color = "blue3", size = 0.5, alpha = 0.5) +      # estações
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(st_bbox(sf.kappa)["xmin"], st_bbox(sf.kappa)["xmax"]), # restringe plotagem p/ extensão de sf.kappa
             ylim = c(st_bbox(sf.kappa)["ymin"], st_bbox(sf.kappa)["ymax"])) + 
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
}); plot.estacoes.br

ggsave(filename = "Plotagens/Análise Exploratória/Mapa Estações Brasil.png",
       plot = plot.estacoes.br,
       height = 12, width = 12, units = "cm",
       dpi = 300)

# BOXPLOT KAPPA -----------------------------------------------------------

# Boxplot usando kappa ajustado por GML
plot.boxplot.kappa.n.serie.gml <- ({
  df.quantis %>%
    mutate(intervalo = cut(x = n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    ggplot(aes(x = intervalo, y = kappa.gml)) +
    geom_jitter(size = 0.3, alpha = 0.3, color = "grey70") +
    stat_boxplot(geom ='errorbar', linewidth = 0.5) +
    geom_boxplot(linewidth = 0.5, fill = "#FFD550", outlier.shape = 4, alpha = 0.6) +
    # scale_y_continuous(limits = c(min(df.quantis$kappa.gml) - 0.1, max(df.quantis$kappa.gml) + 0.1)) +
    ylim(c(-0.6, 0.6)) +
    labs(x = "N [anos]", y = "κ", title = "b)") +
    theme_minimal() +
    theme(panel.background = NULL,
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.boxplot.kappa.n.serie.gml

plot.kappa.gml <- ({ plot.kappa.n.serie.gml / plot.boxplot.kappa.n.serie.gml }); plot.kappa.gml

ggsave(filename = "Plotagens/Ajuste das Estações/Kappa (GMLE).png",
       plot = plot.kappa.gml,
       width = 10, height = 14, units = "cm",
       dpi = 300)

# Boxplot usando kappa ajustado por GML
plot.boxplot.kappa.n.serie.lmom <- ({
  df.quantis %>%
    mutate(intervalo = cut(x = n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    ggplot(aes(x = intervalo, y = kappa.lmom)) +
    geom_jitter(size = 0.3, alpha = 0.3, color = "grey70") +
    stat_boxplot(geom ='errorbar', linewidth = 0.5) +
    geom_boxplot(linewidth = 0.5, fill = "#FFD550", outlier.shape = 4, alpha = 0.6) +
    scale_y_continuous(limits = c(min(df.quantis$kappa.lmom) - 0.1, max(df.quantis$kappa.lmom) + 0.1)) +
    labs(x = "N [anos]", y = "κ", title = "b)") +
    theme_minimal() +
    theme(panel.background = NULL,
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.boxplot.kappa.n.serie.lmom

plot.kappa.lmom <- ({ plot.kappa.n.serie.lmom / plot.boxplot.kappa.n.serie.lmom }); plot.kappa.lmom

ggsave(filename = "Plotagens/Ajuste das Estações/Kappa (MOM-L).png",
       plot = plot.kappa.lmom,
       width = 10, height = 15, units = "cm",
       dpi = 300)

# BOXPLOT KAPPA RH --------------------------------------------------------

n <- 30

# Gráfico auxilizar p/ ordenar as bacias por kappa.gml
median_gml <- df.quantis %>% 
  filter(n.serie >= n) %>% 
  group_by(bacia.codigo) %>% 
  summarise(mediana.kappa.gml = median(kappa.gml, na.rm = TRUE)) %>%
  arrange(mediana.kappa.gml)

# Gráfico de barras apresentando o percentual de estações com kappa < 0 por estação
plot.hist.kappa.rh <- ({
  df.quantis %>% 
    filter(n.serie >= n) %>% 
    pivot_longer(cols = starts_with("kappa"), # transforma em formato long
                 names_to = "metodo",
                 names_prefix = "kappa.",
                 values_to = "kappa") %>% 
    mutate(metodo = factor(toupper(metodo), levels = c("LMOM", "ML", "GML")),
           bacia.codigo = as.integer(bacia.codigo)) %>% 
    group_by(bacia.codigo, metodo) %>% 
    summarise(media.kappa = mean(kappa, na.rm = TRUE),               # média de kappa p/ cada estimador
              prct.kappa.neg = 100*sum(kappa < 0, na.rm = TRUE)/n(), # % kappa < 0 p/ cada estimador
              n.estacoes = n()) %>%                                  # nro estações em cada bacia  
    ungroup() %>% 
    ggplot(aes(x = factor(bacia.codigo, levels = median_gml$bacia.codigo), 
               y = media.kappa, fill = metodo)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3, alpha = 0.6) +
    # Rótulo do valor médio de kappa
    geom_text(aes(label = format(round(media.kappa, digits = 3), digits = 3), family = "serif"), 
              size = 2.8, position = position_dodge(width = 0.9), vjust = -0.7) +
    # Rótulo da % de kappa negativos
    geom_text(aes(label = sprintf("%.0f%%", prct.kappa.neg), family = "serif", fontface = "bold"), 
              size = 3, position = position_dodge(width = 0.9), vjust = 1.6) +
    # Rótulo do número de estações
    geom_text(data = . %>% filter(metodo == "ML"), # filtrando somente ML p/ rotulo ficar no meio
              aes(label = n.estacoes[metodo == "ML"], y = 0.06), family = "serif",
              size = 3.8, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(-0.07, 0.06)) + # p/ inverter usar trans = "reverse"
    scale_x_discrete(labels = setNames(df.rh$rh, df.rh$id)) +
    scale_fill_manual(values = c("GML" = "steelblue", "LMOM" = "#FFD550", "ML" = "grey60")) +
    labs(y = "Valor médio de κ", x = "", title = "a)", fill = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(color = "white", fill = "white"),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          text = element_text(family = "serif", size = 11, color = "black"))
}); plot.hist.kappa.rh

# Plot comparando os resutados de LMOM e GML por região
plot.boxplot.kappa.rh <- ({
  df.quantis %>% 
    filter(n.serie >= n) %>% 
    select(estacao, n.serie, bacia.codigo, kappa.lmom, kappa.gml, kappa.ml) %>% 
    pivot_longer(cols = starts_with("kappa"), # transforma em formato long
                 names_to = "metodo",
                 names_prefix = "kappa.",
                 values_to = "kappa") %>% 
    mutate(metodo = factor(toupper(metodo), levels = c("LMOM", "ML", "GML"))) %>% 
    ungroup() %>%
    mutate(bacia.codigo = factor(bacia.codigo, levels = median_gml$bacia.codigo)) %>% # ordena pela mediana
    ggplot(aes(x = bacia.codigo, y = kappa, fill = metodo)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.3) +
    # geom_jitter(aes(color = metodo), size = 0.3, shape = 21, alpha = 0.6,) + # diferenciar cores 
    stat_boxplot(aes(group = interaction(bacia.codigo, metodo)),
                 geom = "errorbar", width = 0.4, position = position_dodge(width = 0.75)) + # tava um pouco deslocado
    geom_boxplot(aes(group = interaction(bacia.codigo, metodo)),                            # interaction p/ plotar lado a lado
                 linewidth = 0.4, alpha = 0.6,
                 outlier.size = 1.5, outlier.alpha = 0.60, outlier.shape = 4) +
    scale_x_discrete(labels = setNames(df.rh$rh, df.rh$id)) +
    labs(x = "", y = "κ", fill = "", title = "b)") +
    scale_fill_manual(values = c(GML = "steelblue", ML = "grey60", LMOM = "#FFD550")) +
    # coord_flip() +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(color = "white", fill = "white"),
          legend.position = "bottom",
          axis.text = element_text(color = "black"),
          # axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          text = element_text(family = "serif", size = 12, color = "black"))
}); plot.boxplot.kappa.rh

plot.kappa.rh <- (plot.hist.kappa.rh + plot_spacer() + plot.boxplot.kappa.rh) + plot_layout(ncol = 1, heights = c(0.8, -0.1, 1.2)); plot.kappa.rh

ggsave(filename = "Plotagens/Ajuste das Estações/Boxplot Kappa Estimadores RH.png",
       plot = plot.boxplot.kappa.rh,
       width = 220, height = 12, units = "cm",
       dpi = 300)

ggsave(filename = "Plotagens/Ajuste das Estações/Kappa Estimadores RH.png",
       plot = plot.kappa.rh,
       width = 24, height = 20, units = "cm",
       dpi = 300)
