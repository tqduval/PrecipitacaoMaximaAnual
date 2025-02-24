
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ---------------------------------------------------------

# Importar df.param caso não esteja no 'Global Environment'
path.df.param <- file.choose()
df.param <- readRDS(path.df.param)

# CÁLCULO DOS QUANTIS -----------------------------------------------------

# Novo dataframe c/ parâmetros ajustados
df.quantis <- ({
  df.param %>% 
  select(estacao, n.serie, p.0,
         bacia.codigo, subbacia.codigo,
         hess.gl,
         param.lmom.gu, max.l.gu,
         param.lmom.gev, max.l.gev, max.gl.gev.beta) %>% 
  rename(param.ml.gu = max.l.gu, param.ml.gev = max.l.gev, param.gml.gev = max.gl.gev.beta) %>% 
  mutate(param.ml.gu = map_if(param.ml.gu, ~ nrow(.) > 1, ~ slice(., 1)),           # algumas colunas de parâmetros têm mais de uma linha
         param.ml.gev = map_if(param.ml.gev, ~ nrow(.) > 1, ~ slice(., 1)),         # no dataframe final pq o otimizador não retornou somente
         param.gml.gev = map_if(param.gml.gev, ~ nrow(.) > 1, ~ slice(., 1)),       # então pega só a primeira linha desses c/ mais de um
         csi.gml.gev = sapply(param.gml.gev, function(estacao) estacao[[1]][[1]]),  # csi GEV (GML)
         alpha.gml.gev = sapply(param.gml.gev, function(estacao) estacao[[2]][[1]]),# alpha GEV (GML)
         kappa.lmom = sapply(param.lmom.gev, function(estacao) estacao[[3]][[1]]),  # kappa GEV (MOM-L)
         kappa.ml = sapply(param.ml.gev, function(estacao) estacao[[3]][[1]]),      # kappa GEV (ML)
         kappa.gml = sapply(param.gml.gev, function(estacao) estacao[[3]][[1]]))})  # kappa GEV (GML)

# Calcular quantis p/ Tr = 10, 20, 50, 100, 200, 500 anos
for (Tr in c(10.0, 20.0, 50.0, 100.0, 200.0, 500.0)){
  
  # Probabilidades de excedência
  p <- 1 - 1 / Tr
  p.ajustado <- (p - df.quantis$p.0)/(1 - df.quantis$p.0) # ajusta a probabilidade conforma a probabilidade de pd.max = 0 (p.0)
  
  # Gumbel LMOM
  gu.lmom <- paste0("q", Tr, ".lmom.gu")
  df.quantis[[gu.lmom]] <- sapply(df.quantis$param.lmom.gu, function(param) fun.q.gu(p.ajustado, param) %>% unlist)
  
  # Gumbel ML
  gu.ml <- paste0("q", Tr, ".ml.gu")
  df.quantis[[gu.ml]] <- sapply(df.quantis$param.ml.gu, function(param) fun.q.gu(p.ajustado, param) %>% unlist)
  
  # GEV LMOM
  gev.lmom <- paste0("q", Tr, ".lmom.gev")
  df.quantis[[gev.lmom]] <- sapply(df.quantis$param.lmom.gev, function(param) fun.q.gev(p.ajustado, param) %>% unlist)
  
  # GEV LMOM
  gev.ml <- paste0("q", Tr, ".ml.gev")
  df.quantis[[gev.ml]] <- sapply(df.quantis$param.ml.gev, function(param) fun.q.gev(p.ajustado, param) %>% unlist)
  
  # GEV GML
  gev.gml <- paste0("q", Tr, ".gml.gev")
  df.quantis[[gev.gml]] <- sapply(df.quantis$param.gml.gev, function(param) fun.q.gev(p.ajustado, param) %>% unlist)
  
}

# Fazer dataframe no formato 'long'
df.quantis_long <- ({
  df.quantis %>% 
  pivot_longer(cols = starts_with("q"),                 # fazer um pivot longer e deixar somente as colunas do tempo de retorno,
               names_to = c("tempo.retorno", "modelo"), # do modelo (combinação de distribuição e método de estimação) e do
               names_pattern = "q(\\d+).(.*)",          # quantil respectivo
               values_to = "quantil") %>% 
  select(-starts_with("param")) %>%                     # remove tudo que começa com "param"
  pivot_wider(names_from = modelo,                      # cria uma coluna para cada modelo e uma linha para cada TR
              values_from = quantil) %>% 
  mutate(tempo.retorno = as.numeric(tempo.retorno),     # transformar Tr de caractere para numérico
         erro.lmom = (lmom.gu - lmom.gev)/lmom.gev, 
         erro.ml = (ml.gu - ml.gev)/ml.gev,
         erro.gml = (ml.gu - gml.gev)/gml.gev)})

# Exportar RDS
saveRDS(object = df.quantis, file = "Dados Gerados/df.quantis.rds")
saveRDS(object = df.quantis_long, file = "Dados Gerados/df.quantis_long.rds")


# CONFERÊNCIA KAPPA ± 0.5 -------------------------------------------------

estacoes.0.5 <- df.quantis.240815$estacao[(df.quantis.240815$kappa.gml >= 0.5 | df.quantis.240815$kappa.gml <= -0.5)]

# GRÁFICOS ----------------------------------------------------------------

# Gráfico Quantil-Quantil com κ < 0 e κ > 0 de cores diferentes e histogramas (Gumbel ML e GEV GML)
{
  limite.n.serie <- 30
  plot.qq.gml.kappa <- list()
  plot.hists.gml.kappa <- list()
  
  for(Tr in c(10, 20, 50 , 100, 200, 500)){
    
    plot.qq.gml <- ggplot(data = df.quantis_long %>%
                            filter(tempo.retorno == Tr,
                                   n.serie >= limite.n.serie)) + # filtrando so estações c/ mais de 30 anos
      # gráfico quantil-quantil
      geom_point(aes(x = gml.gev, y = ml.gu, color = factor(sign(kappa.gml))),
                 alpha = 0.3, size = 2.5, stroke = 0, shape = 20) +
      geom_abline(color = "black", linetype = "dashed", linewidth = 0.4) +
      xlim(c(0, max(c(df.quantis_long$ml.gu[which(df.quantis_long$tempo.retorno == Tr)],
                      df.quantis_long$ml.gu[which(df.quantis_long$tempo.retorno == Tr)]), na.rm = TRUE))) +
      ylim(c(0, max(c(df.quantis_long$ml.gu[which(df.quantis_long$tempo.retorno == Tr)],
                      df.quantis_long$ml.gu[which(df.quantis_long$tempo.retorno == Tr)]), na.rm = TRUE))) +
      scale_color_manual(values = c("-1" = "steelblue", "1" = "darkorange"),
                         labels = c("-1" = "κ < 0", "1" = "κ > 0")) +
      coord_fixed() +
      labs(x = "GEV [mm]",
           y = "Gumbel [mm]",
           title = paste(Tr, "anos")) +
      theme_minimal() +
      theme(panel.background = NULL,
            legend.position = "none",
            text = element_text(family = "serif",
                                size = 15,
                                color = "black")) # especifica uma fonte
    
    # histogramas de erro
    plot.hist.erro.gml <- ggplot(data = df.quantis_long %>%
                                   filter(tempo.retorno == Tr,
                                          n.serie >= limite.n.serie)) +
      geom_histogram(aes(x = erro.gml,
                         y = after_stat(count/sum(count))),
                     fill = "orange",
                     color = "#c15123") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1L),
                         limits = c(0, 0.4)) +
      geom_vline(xintercept = 0,
                 color = "red",
                 linewidth = 0.6,
                 linetype = "dashed") +
      labs(x = "Diferença relativa",
           y = "Porcentagem de estações [%]") +
      theme_minimal() +
      theme(panel.background = NULL,
            legend.position = "none",
            text = element_text(family = "serif",
                                size = 15,
                                color = "black"))
    
    plot.qq.gml.kappa[[Tr]] <- plot.qq.gml
    plot.hists.gml.kappa[[Tr]] <- plot.hist.erro.gml
    
  }
  
  plot.qq.erro.gml <- # combinar gráficos usando patchwork
    (plot.qq.gml.kappa[[10]] | plot.hists.gml.kappa[[10]] | plot.qq.gml.kappa[[100]] | plot.hists.gml.kappa[[100]])/
    (plot.qq.gml.kappa[[20]] | plot.hists.gml.kappa[[20]] | plot.qq.gml.kappa[[200]] | plot.hists.gml.kappa[[200]])/
    (plot.qq.gml.kappa[[50]] | plot.hists.gml.kappa[[50]] | plot.qq.gml.kappa[[500]] | plot.hists.gml.kappa[[500]])
  
  
  print(plot.qq.erro.gml)
  
  ggsave(file = "Plotagens/Quantis/Gráficos Quantil-Quantil e Erro GML Kappa 30 anos-teste.png",
         plot = plot.qq.erro.gml,
         width = 14,
         height = 12,
         dpi = 300)
  
} # final bloco GML

# Distribuição acumulada diferença relativa
plot.dr.acum <- ({
  ggplot(data = df.quantis_long, aes(x = erro.gml, color = as.factor(tempo.retorno))) +
    stat_ecdf() +
    scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.1)) +
    scale_y_continuous(labels = scales::label_percent(scale = 100), breaks = seq(0, 1, 0.1)) +
    # geom_vline(xintercept = -0.15) +
    labs(color = "", x = "Diferença Relativa", y = "Estações [%]") +
    scale_color_manual(values = c("10" = "magenta", "20" = "purple", "50" = "blue", "100" = "steelblue", "200" = "#469990", "500" = "green4"),
                       labels = c("10" = "TR = 10 anos", "20" = "TR = 20 anos", "50" = "TR = 50 anos", "100" = "TR = 100 anos", "200" = "TR = 200 anos", "500" = "TR = 500 anos")) +
    guides(color = guide_legend(nrow = 2)) +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill="white", color = "white"),
          legend.position = "bottom",
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.dr.acum

ggsave(filename = "Plotagens/Quantis/Diferença Relativa TR.png",
       plot = plot.dr.acum,
       width = 15, height = 12, units = "cm",
       dpi = 300)

# GRÁFICO QUANTIS TR10 RH -------------------------------------------------

library(grid) # setas

# Estações c/ quantil de 10 anos >= 300 mm
est.q10.300 <- 
  df.quantis_long %>% 
  select(estacao, bacia.codigo, tempo.retorno, gml.gev) %>% 
  mutate(bacia.codigo = as.integer(bacia.codigo)) %>% 
  filter(tempo.retorno == 10, gml.gev >= 300)

# Texto outliers
txt.out.rh1 <- paste0("[", paste0(est.q10.300 %>% filter(bacia.codigo == 1) %>% pull(gml.gev) %>% round, collapse = "; "), "]")
txt.out.rh8 <- paste0("[", est.q10.300 %>% filter(bacia.codigo == 8) %>% pull(gml.gev) %>% round, "]")

# Boxplot dos quantis obtidos por região p/ determinado Tr
Tr <- 10
plot.boxplot.q10 <- ({
  df.quantis_long %>% 
    filter(tempo.retorno == Tr) %>% 
    group_by(bacia.codigo) %>% 
    mutate(mediana.q = median(gml.gev),
           bacia.codigo = as.integer(bacia.codigo)) %>%
    ungroup() %>%
    mutate(bacia.codigo = factor(bacia.codigo, levels = unique(bacia.codigo[order(mediana.q)]))) %>% # ordena pela mediana
    ggplot(aes(x = bacia.codigo)) +
    geom_jitter(aes(group = bacia.codigo, y = gml.gev), size = 0.3, shape = 20, alpha = 0.3, color = "grey60") +
    stat_boxplot(aes(group = bacia.codigo, y = gml.gev), geom = "errorbar", width = 0.4) +
    geom_boxplot(aes(group = bacia.codigo, y = gml.gev), linewidth = 0.4, fill = "#FFD550", alpha = 0.75,
                 outlier.size = 1.5, outlier.alpha = 0.60, outlier.shape = 4) +
    scale_x_discrete(breaks = seq(1, 8, 1), labels = setNames(df.rh$rh, df.rh$id)) +
    # coord_flip() + # colocar
    ylim(c(0, 300)) + # tentar indicar quantas estações ficaram de fora
    labs(x = "", y = "Precipitação Diária Máxima Anual (Tr = 10 anos) [mm]") +
    # Indicar outliers fora do plot
    annotate("segment", x = 5.2, xend = 6, y = 250, yend = 300, arrow = arrow(type = "open", length = unit(0.18, "cm")), color = "grey20", alpha = 0.7) +
    annotate("text", x = 5.2, y = 235, label = txt.out.rh1, size = 4, hjust = 0.5, vjust = 0, family = "serif") +
    annotate("segment", x = 8, xend = 8, y = 255, yend = 300, arrow = arrow(type = "open", length = unit(0.18, "cm")), color = "grey20", alpha = 0.7) +
    annotate("text", x = 8, y = 240, label = txt.out.rh8, size = 4, hjust = 0.5, vjust = 0, family) +
    theme_minimal() +
    theme(panel.background = NULL,
          title = element_text(size = 10),
          plot.background = element_rect(color = "white", fill = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"))
}); plot.boxplot.q10

ggsave(filename = "Plotagens/Quantis/Boxplot Quantil 10 anos RH.png",
       plot = plot.boxplot.q10,
       width = 20, height = 12, units = "cm",
       dpi = 300)
