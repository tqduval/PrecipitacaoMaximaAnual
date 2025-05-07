
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, readxl, sf,  ggplot2, patchwork, beepr)


# LER DADOS ---------------------------------------------------------------

# Rio Grande do Sul - Daily
df_rs_completo_path <- file.choose()
df_rs_completo <- read_excel(df_rs_completo_path)
df_rs <- ({
  df_rs_completo %>% 
    mutate(
      Date = year,
      Estacao_codigo = gauge_code,
      Pdmax = peak_24h,
      BaciaCodigo = NA,
      SubBaciaCodigo = NA,
      Altitude = elevation,
      AreaDrenagem = NA,
      NomedoRio = NA,
      Latitude = lat,
      Longitude = long) %>%  # reestruturar connforme df.precip.max.anual, depois padronizar os nomes conforme
    select(c(Date, Estacao_codigo, Pdmax, BaciaCodigo, SubBaciaCodigo, Altitude, NomedoRio, Latitude, Longitude))
})
sf_rs_completo <-
  df_rs_completo %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4674) %>% 
  st_write(dsn = "Dados Gerados/GIS/sf_rs_daily.gpkg")

# Rio Grande do Sul - State shapefile
sf_rs_path <- file.choose()
sf_rs <- st_read(dsn = sf.rs.path) %>% 
  filter(SIGLA_UF == "RS")

# AJUSTAR DISTRIBUIÇÕES ---------------------------------------------------

# Ajustar parâmetros (erro no l.gev estação 03152016, parâmetros iniciais)
df_rs_param <- fun.process(df = df_rs,
                           l.ratio.alpha = 0.05)

# Calcualr quantis (somente gl.gev)
df_rs_q <- fun.q.process.all(df = df_rs_param,
                             tempo.retorno = c(10, 20, 50, 100, 200, 500),
                             model = c("gl.gev", "lmom.gev"))


# VISUALIZAÇÃO ------------------------------------------------------------

# Shape parameter distribution (histogram) for gl.gev model
df_rs_param <- 
  df_rs_param %>% 
  mutate(kappa.gl = sapply(par.gl.gev, function(par) par[[3]][[1]]),
         kappa.lmom = sapply(par.lmom.gev, function(par) par[[3]][[1]])) # extract shape parameter in separate col

mu_rs_gl <- mean(df_rs_param$kappa.gl)
sd_rs_gl <- sd(df_rs_param$kappa.gl)

mu_rs_lmom <- mean(df_rs_param$kappa.lmom)
sd_rs_lmom <- sd(df_rs_param$kappa.lmom)


# HISTOGRAMA DE KAPPA -----------------------------------------------------

# LMOM
plot_hist_kappa_rs_lmom <- ({
  
  ggplot(data = df_rs_param, aes(x = kappa.lmom)) +
    geom_histogram(aes(y = after_stat(density)),
                   color = "#c15123", fill = "orange", alpha = 0.8,
                   bins = 50) +
    # priori Beta Martins e Stedinger (2000)
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = -0.10, sd = 0.122), 
                  linewidth = 0.5,
                  aes(color = "Beta geofísica", linetype = "Beta geofísica")) +
    # normal empírica BR
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = -0.038, sd = 0.084),
                  linewidth = 0.5,
                  aes(color = "Beta empírica BR", linetype = "Beta empírica BR")) + 
    # beta empírica RS (GML)
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = mu_rs_gl, sd = sd_rs_gl), 
                  linewidth = 0.7,
                  aes(color = "Beta empírica RS (GML)", linetype = "Beta empírica RS (GML)")) +
    # beta empírica RS (LMOM)
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = mu_rs_lmom, sd = sd_rs_lmom), 
                  linewidth = 0.7,
                  aes(color = "Beta empírica RS (LMOM)", linetype = "Beta empírica RS (LMOM)")) +
    annotate("text", x = 0.5, y = 7, hjust = 1, vjust = 1, size = 4, family = "serif",
             label = paste0("Beta(", round(mu_rs_lmom, 3), ", ", round(sd_rs_lmom, 3), "²)")) +
    scale_color_manual(values = c("Beta empírica RS (GML)" = "black",
                                  "Beta empírica RS (LMOM)" = "orange4",
                                  "Beta empírica BR" = "red3",
                                  "Beta geofísica" = "magenta")) +
    scale_linetype_manual(values = c("Beta empírica RS (GML)" = "solid",
                                     "Beta empírica RS (LMOM)" = "solid",
                                     "Beta empírica BR" = "solid",
                                     "Beta geofísica" = "solid")) + 
    scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1), limits = c(-0.5, 0.5)) +
    scale_y_continuous(breaks = seq(0, 7, 2), limits = c(0, 7)) +
    labs(x = "κ", y = "Densidade", color = "", linetype = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "bottom",
          legend.text = element_text(family = "serif", size = 12, color = "black"),
          text = element_text(family = "serif", size = 12, color = "black")) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), # quebrar legenda em duas linhas
           linetypee = guide_legend(nrow = 2, byrow = TRUE))
  
})

ggsave(filename = "Plotagens/RS/lmom/Distribuição de Kappa no Rio Grande do Sul (24h).png",
       plot = plot_hist_kappa_rs_lmom,
       height = 10, width = 15, units = "cm",
       dpi = 300)

# GML
plot_hist_kappa_rs_gl <- ({
  
  ggplot(data = df_rs_param, aes(x = kappa.gl)) +
    geom_histogram(aes(y = after_stat(density)),
                   color = "#c15123", fill = "orange", alpha = 0.8,
                   bins = 50) +
    # priori Beta Martins e Stedinger (2000)
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = -0.10, sd = 0.122), 
                  linewidth = 0.5,
                  aes(color = "Beta geofísica", linetype = "Beta geofísica")) +
    # normal empírica BR
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = -0.038, sd = 0.084),
                  linewidth = 0.5,
                  aes(color = "Beta empírica BR", linetype = "Beta empírica BR")) + 
    # beta empírica RS
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = mu_rs_gl, sd = sd_rs_gl), 
                  linewidth = 0.7,
                  aes(color = "Beta empírica RS (GML)", linetype = "Beta empírica RS (GML)")) +
    annotate("text", x = 0.5, y = 7, hjust = 1, vjust = 1, size = 4, family = "serif",
             label = paste0("Beta(", round(mu.rs, 3), ", ", round(sd.rs, 3), "²)")) +
    scale_color_manual(values = c("Beta empírica RS (GML)" = "black",
                                  "Beta empírica BR" = "red3",
                                  "Beta geofísica" = "magenta")) +
    scale_linetype_manual(values = c("Beta empírica RS (GML)" = "solid",
                                     "Beta empírica BR" = "solid",
                                     "Beta geofísica" = "solid")) + 
    scale_x_continuous(breaks = seq(-0.5, 0.5, 0.1), limits = c(-0.5, 0.5)) +
    scale_y_continuous(breaks = seq(0, 7, 2), limits = c(0, 7)) +
    labs(x = "κ", y = "Densidade", color = "", linetype = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "bottom",
          legend.text = element_text(family = "serif", size = 12, color = "black"),
          text = element_text(family = "serif", size = 12, color = "black")) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE),
           linetypee = guide_legend(nrow = 2, byrow = TRUE))
  
})

ggsave(filename = "Plotagens/RS/gml/Distribuição de Kappa no Rio Grande do Sul (24h).png",
       plot = plot_hist_kappa_rs_gl,
       height = 10, width = 15, units = "cm",
       dpi = 300)


# MAPA DISTRIBUIÇÃO DE KAPPA ----------------------------------------------

# Shape parameter values accross Rio Grande do Sul
sf_rs_param <- st_as_sf(x = df_rs_param,
                        coords = c("long", "lat"),
                        crs = 4674) # conferir depois

st_write(obj = sf_rs_param,
         dsn = "Dados Gerados/SIG/sf_rs_param.gpkg")

# GML
plot_map_kappa_rs_gl <- ({
  
  ggplot() +
    geom_sf(data = sf.rs, fill = "white", color = "black", linewidth = 0.5) + # RS
    geom_sf(data = sf.rs.param,
            aes(color = kappa.gl), # col with kappa values
            size = 2.5, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", name = "κ", limits = c(-0.3, 0.5)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(st_bbox(sf.rs.param)["xmin"], st_bbox(sf.rs.param)["xmax"]), # restringe plotagem p/ extensão de sf.kappa
             ylim = c(st_bbox(sf.rs.param)["ymin"], st_bbox(sf.rs.param)["ymax"])) + 
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
  
})

ggsave(filename = "Plotagens/RS/gml/Distribuição Regional de Kappa no Rio Grande do Sul (24h).png",
       plot = plot_map_kappa_rs_gl,
       height = 12, width = 15, units = "cm",
       dpi = 300)

# LMOM (como ele troca os intervalos, talvez seja melhor dividir entre positivo e negativo)
plot_map_kappa_rs_lmom <- ({
  
  ggplot() +
    geom_sf(data = sf.rs, fill = "white", color = "black", linewidth = 0.5) + # RS
    geom_sf(data = sf.rs.param,
            aes(color = kappa.lmom), # col with kappa values
            size = 2.5, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", name = "κ", limits = c(-0.3, 0.5)) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude") +
    coord_sf(xlim = c(st_bbox(sf.rs.param)["xmin"], st_bbox(sf.rs.param)["xmax"]), # restringe plotagem p/ extensão de sf.kappa
             ylim = c(st_bbox(sf.rs.param)["ymin"], st_bbox(sf.rs.param)["ymax"])) + 
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
  
})

ggsave(filename = "Plotagens/RS/lmom/Distribuição Regional de Kappa no Rio Grande do Sul (24h).png",
       plot = plot_map_kappa_rs_lmom,
       height = 12, width = 15, units = "cm",
       dpi = 300)


# KAPPA POR TAMANHO DAS SÉRIES --------------------------------------------

# Shape parameter by time series length - Papalexiou and Koutsoyiannis (2013) (GML)
plot_kappa_by_length_gl <- ({
  df_rs_param %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    group_by(interval) %>% 
    summarise(mean.kappa = mean(kappa.gl),
              prct.kappa.neg = sum(kappa.gl < 0)/n()*100) %>% 
    ungroup %>% 
    ggplot(aes(x = interval, y = mean.kappa)) +
    geom_bar(stat = "identity", color = "#c15123", fill = "#FFD550") +
    geom_text(aes(label = format(mean.kappa, digits = 2), family = "serif"),
              size = 3.5, position = position_stack(), vjust = -0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prct.kappa.neg), family = "serif", fontface = "bold"),
              size = 3.5, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(limits = c(0.16, -0.18), trans = "reverse") +
    labs(x = "", y = "Valor médio de κ") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          # axis.text.x = element_blank(),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
  
})

plot_boxplot_kappa_by_length_gl <- ({
  df_rs_param %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    group_by(interval) %>% 
    ggplot(aes(x = interval, y = kappa.gl)) +
    geom_jitter(size = 0.3, color = "grey70") +
    stat_boxplot(geom = "errorbar", linewidth = 0.5) +
    geom_boxplot(linewidth = 0.5, fill = "#FFD550", outlier.shape = 4, alpha = 0.6) +
    ylim(c(-0.2, 0.2)) +
    labs(x = "N[anos]", y = "κ") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
})

plot_kappa_gl_rs <- ({plot_kappa_by_length_gl / plot_boxplot_kappa_by_length_gl})

ggsave(filename = "Plotagens/RS/Kappa médio por tamanho da série (24h).png",
       plot = plot_kappa_gl_rs,
       height = 20, width = 15, units = "cm",
       dpi = 300)

# Shape parameter by time series length - Papalexiou and Koutsoyiannis (2013) (LMOM)
plot_kappa_by_length_lmom <- ({
  df_rs_param %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    group_by(interval) %>% 
    summarise(mean.kappa = mean(kappa.lmom),
              prct.kappa.neg = sum(kappa.lmom < 0)/n()*100) %>% 
    ungroup %>% 
    ggplot(aes(x = interval, y = mean.kappa)) +
    geom_bar(stat = "identity", color = "#c15123", fill = "#FFD550") +
    geom_text(aes(label = format(mean.kappa, digits = 2), family = "serif"),
              size = 3.5, position = position_stack(), vjust = -0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prct.kappa.neg), family = "serif", fontface = "bold"),
              size = 3.5, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(limits = c(0.16, -0.18), trans = "reverse") +
    labs(x = "", y = "Valor médio de κ") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          # axis.text.x = element_blank(),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
  
})

plot_boxplot_kappa_by_length_lmom <- ({
  df_rs_param %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    group_by(interval) %>% 
    ggplot(aes(x = interval, y = kappa.lmom)) +
    geom_jitter(size = 0.3, color = "grey70") +
    stat_boxplot(geom = "errorbar", linewidth = 0.5) +
    geom_boxplot(linewidth = 0.5, fill = "#FFD550", outlier.shape = 4, alpha = 0.6) +
    ylim(c(-0.2, 0.2)) +
    labs(x = "N[anos]", y = "κ") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 12,
                              color = "black"))
})

plot_kappa_lmom_rs <- ({plot_kappa_by_length_lmom / plot_boxplot_kappa_by_length_lmom})

ggsave(filename = "Plotagens/RS/lmom/Kappa médio por tamanho da série (24h).png",
       plot = plot_kappa_lmom_rs,
       height = 20, width = 15, units = "cm",
       dpi = 300)


# INCERTEZA ---------------------------------------------------------------

# Relative CI width
df_rs_q <- 
  df_rs_q %>% 
  mutate(relative_ci_width = (ic.u - ic.l)/q)

plot_boxplot_ci_width <- ({
  df_rs_q %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    ggplot() +
    facet_wrap(~tempo.retorno, labeller = labeller(tempo.retorno = function(tr) paste(tr, "anos"))) +
    stat_boxplot(aes(x = interval, y = relative_ci_width), geom = "errorbar", width = 0.6) +
    geom_jitter(aes(x = interval, y = relative_ci_width), alpha = 0.7, size = 0.3, color = "grey70") +
    geom_boxplot(aes(x = interval, y = relative_ci_width), alpha = 0.6, linewidth = 0.4, fill = "#a086d2", outlier.shape = 4) +
    ylim(c(0, 1.15)) +
    labs(x = "N[anos]", y = "Largura relativa IC 95%") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 10, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
})

ggsave(filename = "Plotagens/RS/Largura relativa do intervalo de confiança (24h).png",
       plot = plot_boxplot_ci_width,
       height = 15, width = 20, units = "cm",
       dpi = 300)


# DISTRIBUIÇÃO TAMANHO SÉRIES ---------------------------------------------

# Time series length
plot_n_serie_by_length <- ({
  df_rs_param %>% 
    mutate(interval = cut(x = n.serie,
                          breaks = c(30, 40, 50, 60, 70, Inf),
                          right = FALSE,
                          labels = c("30-40", "40-50", "50-60", "60-70", ">70"))) %>% 
    ggplot(aes(x = interval)) +
    geom_bar(color = "steelblue", fill = "lightblue") +
    # ylim(c()) +
    geom_text(stat = "count",
              aes(label = ..count..),
              vjust = -0.5, family = "serif") +
    labs(x = "N[anos]", y = "N° estações") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 10, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
})

ggsave(filename = "Plotagens/RS/Número de estações por duração das séries (24h).png",
       plot = plot_n_serie_by_length,
       height = 12, width = 15, units = "cm",
       dpi = 300)
