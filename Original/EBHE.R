
# GRÁFICOS EBHE INGLÊS ----------------------------------------------------

# Figure 3. Comparison K values between basins and estimators ####

# Gráfico de barras apresentando o percentual de estações com kappa < 0 por estação
plot.hist.kappa.rh.en <- ({
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
    labs(y = "Average value of κ estimates", x = "", title = "a)", fill = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(color = "white", fill = "white"),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          text = element_text(family = "serif", size = 11, color = "black"))
}); plot.hist.kappa.rh.en

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

plot.kappa.rh.en <- (plot.hist.kappa.rh.en + plot_spacer() + plot.boxplot.kappa.rh) + plot_layout(ncol = 1, heights = c(0.8, -0.1, 1.2)); plot.kappa.rh.en

ggsave(filename = "Plotagens/EBHE/Kappa Estimadores RH.png",
       plot = plot.kappa.rh.en,
       width = 24, height = 20, units = "cm",
       dpi = 300)

# Figure 4. Geophysical prior ####
plot.hist.kappa.en <- ({
  
  # N > 40
  ggplot(data = df.quantis %>% 
           filter(kappa.lmom >= -0.5 & kappa.lmom <= 0.5,
                  n.serie >= 40),
         aes(x = kappa.lmom)) +
    
    # # Histograma c/ densidade
    # geom_histogram(aes(y = after_stat(density)),
    #                color = "#c15123", fill = "orange",
    #                bins = 40) +
    
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
    
    # annotate("text", x = 0.5, y = 5, hjust = 1, vjust = 1, family = "serif", size = 5,
    #          label = paste0("N(", round(mean(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40]), 3),
    #                         ", ", round(sd(df.quantis$kappa.lmom[df.quantis$kappa.lmom >= -0.5 & df.quantis$kappa.lmom <= 0.5 & df.quantis$n.serie >= 40]), 3), "²)")) +
    
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
    labs(x = "GEV shape parameter κ", y = "Density function", color = "", linetype = "") +
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
}); plot.hist.kappa.en

ggsave(filename = "Plotagens/EBHE/Histograma Kappa Tamanho da Série 40 anos.png",
       plot = plot.hist.kappa.en,
       height = 15, width = 25, units = "cm",
       dpi = 300)

# Figure 5. Accumulated relativa quantile difference ####
plot.dr.acum.en <- ({
  ggplot(data = df.quantis_long, aes(x = erro.gml, color = as.factor(tempo.retorno))) +
    stat_ecdf() +
    scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.1)) +
    scale_y_continuous(labels = scales::label_percent(scale = 100), breaks = seq(0, 1, 0.1)) +
    # geom_vline(xintercept = -0.15) +
    labs(color = "", x = "Relative quantile difference", y = "Gauges [%]") +
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
}); plot.dr.acum.en

ggsave(filename = "Plotagens/EBHE/Diferença Relativa TR.png",
       plot = plot.dr.acum.en,
       width = 15, height = 12, units = "cm",
       dpi = 300)

# Figure 6. 95% CI's relative width ####
# Boxplot largura relativa dos intervalos de confiança
plot.boxplot.lr.en <- ({
  df.ic.quantis %>% 
    # filter(!between(kappa.gml, -0.001, 0.001)) %>%
    mutate(intervalo = cut(n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    ggplot() +
    facet_wrap(~tempo.retorno, labeller = labeller(tempo.retorno = function(x) paste(x, "years"))) +
    stat_boxplot(aes(x = intervalo, y = l.r), geom ='errorbar', width = 0.6) +
    geom_jitter(aes(x = intervalo, y = l.r), alpha = 0.2, size = 0.3, color = "grey70") +
    geom_boxplot(aes(x = intervalo, y = l.r), alpha = 0.6, linewidth = 0.4, fill = "#a086d2",
                 outlier.shape = 4) +
    ylim(c(0, 1.6)) +
    labs(x = "N [years]", y = "Relative width of 95% CI") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
}); plot.boxplot.lr.en

ggsave(filename = "Plotagens/EBHE/Boxplot Largura Relativa IC por TR.png",
       plot = plot.boxplot.lr.en,
       width = 16, height = 16, units = "cm",
       dpi = 300)
