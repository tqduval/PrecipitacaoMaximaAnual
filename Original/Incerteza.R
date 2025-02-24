
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# IMPORTAR DADOS ----------------------------------------------------------

# Importar dados caso não estejam no 'Global Environment'
path.df.quantis_long <- file.choose()
df.quantis_long <- readRDS(path.df.quantis_long)

# INTERVALO DE CONFIANÇA PARÂMETROS --------------------------------------------------

alfa.ic <- 0.05
df.ic.param <- ({
  df.quantis %>% 
    select(-starts_with(c("erro.", "param.", "q")), -starts_with("param."), -c(kappa.ml, kappa.lmom)) %>% 
    mutate(
      
      # Separar estações por tamanho das séries
      intervalo = cut(x = n.serie,
                      breaks = c(15, 25, 35, 45, 55, Inf),
                      right = FALSE,
                      labels = c("15-25", "25-35", "35-45", "45-55", "≥55")),
      
      # Transforma de 1x6 p/ 3x3 e inverte a matriz
      hess.gl = lapply(hess.gl, function(vetor) matrix(vetor, ncol = 3, nrow = 3, byrow = TRUE)),
      cov.gl = lapply(hess.gl, solve),
      
      # Extrai os elementos da diagonal principal da matriz inversa
      inv.fisher.csi = sapply(cov.gl, function(estacao) estacao[[1,1]]),   # diagonal principal: csi
      inv.fisher.alpha = sapply(cov.gl, function(estacao) estacao[[2,2]]), # diagonal principal: alpha
      inv.fisher.kappa = sapply(cov.gl, function(estacao) estacao[[3,3]]), # diagonal principal: kappa
      
      # Intervalos de confiança dos parâmetros p/ nível de confiança de 95%
      ic.e.csi = csi.gml.gev - qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.csi),       # intervalo esquerdo IC p/ csi
      ic.d.csi = csi.gml.gev + qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.csi),       # intervalo direito IC p/ csi
      ic.e.alpha = alpha.gml.gev - qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.alpha), # intervalo esquerdo IC p/ alpha
      ic.d.alpha = alpha.gml.gev + qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.alpha), # intervalo direito IC p/ alpha
      ic.e.kappa = kappa.gml - qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.kappa),     # intervalo esquerdo IC p/ kappa
      ic.d.kappa = kappa.gml + qnorm(1 - alfa.ic/2, 0, 1)*sqrt(inv.fisher.kappa))     # intervalo direito IC p/ kappa           
})

# Número de estações c/ kappa pequeno
df.ic.param[df.ic.param$kappa.gml > -0.001 & df.ic.param$kappa.gml < 0.001,] %>% nrow
df.ic.quantis[df.ic.quantis$ic.e < 0,] %>% nrow

# INTERVALO DE CONFIANÇA QUANTIS -----------------------------------------------------------

library(purrr)

alfa.ic <- 0.05
df.ic.quantis <- ({
  df.quantis_long %>% 
    select(-ends_with(c(".gu", ".lmom", ".ml")), -starts_with("erro"), -c(lmom.gev, ml.gev)) %>% 
    mutate(
      
      # Criar identificador
      id = paste0(estacao, ".", tempo.retorno),
      
      # Transforma de 1x6 p/ 3x3 e inverte a matriz
      hess.gl = lapply(hess.gl, function(vetor) matrix(vetor, ncol = 3, nrow = 3, byrow = TRUE)),
      cov.gl = lapply(hess.gl, solve), # no final essa já é a inversa da hessiana
      
      # Atualiza "p" conforme a p.0
      p.ajustado = ((1 - 1/tempo.retorno) - p.0)/(1 - p.0), # p = 1 - 1/Tr
      
      # Derivadas
      deriv.csi = 1,
      deriv.alpha = fun.deriv.alpha(p = p.ajustado, kappa = kappa.gml),
      deriv.kappa = fun.deriv.kappa(p = p.ajustado, alpha = alpha.gml.gev, kappa = kappa.gml),
      # Derivadas (linha a linha)
      # deriv.alpha = mapply(fun.deriv.alpha, p = p.ajustado, kappa = kappa.gml),
      # deriv.kappa = mapply(fun.deriv.kappa, p = p.ajustado, alpha = alpha.gml.gev, kappa = kappa.gml),
      
      # Matriz gradiente ∇yp
      grad.q = purrr::pmap(list(deriv.csi, deriv.alpha, deriv.kappa), c), # usar pmap() pq somente c() estava juntando todas as 'deriv' p/ mesma estação
      
      # Variância do quantil Var[yp]
      sd.q = purrr::pmap(list(grad.q, cov.gl), function(grad_q, cov_gl){t(grad_q) %*% cov_gl %*% grad_q}) %>% as.numeric %>% sqrt,
      
      # Intervalos de confiança usando nível de confiança de 95%
      ic.e = gml.gev - qnorm(1 - alfa.ic/2, 0, 1)*sd.q,
      ic.d = gml.gev + qnorm(1 - alfa.ic/2, 0, 1)*sd.q
      )
})

# GRÁFICOS ----------------------------------------------------------------

# Desvio padrão kappa vs tamanho das séries
# plot.sd.kappa <- ({
#   
#   ggplot(data = df.ic.param) +
#     geom_point(aes(x = n.serie, y = sqrt(inv.fisher.kappa)),
#                alpha = 0.5, color = "orange")
#     # geom_smooth(stat_smooth(x = n.serie, y = sqrt(inv.fisher.kappa)))
#   
# }); plot.sd.kappa

# Boxplot desvio padrão kappa vs tamanho das séries
plot.boxplot.sd.kappa <- ({
  ggplot(data = df.ic.param) +
    stat_boxplot(aes(x = intervalo, y = sqrt(inv.fisher.kappa)), geom ='errorbar', width = 0.6) +
    geom_jitter(aes(x = intervalo, y = sqrt(inv.fisher.kappa)), alpha = 0.2, size = 0.3, color = "grey70") +
    geom_boxplot(aes(x = intervalo, y = sqrt(inv.fisher.kappa)), fill = "#a086d2", linewidth = 0.5, alpha = 0.6, outlier.shape = 4) +
    labs(x = "N [anos]", y = expression(sigma[kappa])) +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif",
                              size = 10,
                              color = "black"))
}); plot.boxplot.sd.kappa

ggsave(filename = "Plotagens/Incerteza/Boxplot Desvio Padrão Kappa.png",
       plot = plot.boxplot.sd.kappa,
       width = 12, height = 7, units = "cm",
       dpi = 300)

# Curvas de intervalo de confiança de kappa p/ cada 
plot.ic.kappa <- ({
  df.ic.param %>%
    group_by(n.serie) %>%
    summarise(media.kappa = mean(kappa.gml),
              media.ic.d = mean(ic.d.kappa),
              media.ic.e = mean(ic.e.kappa)) %>%
    ggplot() +
    geom_smooth(aes(x = n.serie, y = media.ic.d), linewidth = 0.6, se = FALSE, color = "darkorange") +
    geom_smooth(aes(x = n.serie, y = media.ic.e), linewidth = 0.6, se = FALSE, color = "darkorange") +
    geom_smooth(aes(x = n.serie, y = media.kappa), linewidth = 0.6, se = FALSE, color = "red") +
    labs(x = "N [anos]", y = expression(kappa)) +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 10, color = "black"))
  }); plot.ic.kappa

df.ic.param %>%
  ggplot(aes(x = n.serie)) +
  geom_ribbon(aes(ymin = predict(loess(ic.e.kappa~n.serie)), ymax = predict(loess(ic.d.kappa~n.serie))),
              alpha = 0.5, color = "black", linetype = "dashed", linewidth = 0.4, fill = "lavender") +
  geom_smooth(aes(y = kappa.gml), se = FALSE, color = "black", linewidth = 1) +
  geom_point(data = df.ic.param %>% group_by(n.serie) %>% summarise(media.kappa = mean(kappa.gml)),
             aes(y = media.kappa), shape = 8, color = "red") +
  labs(x = "N [anos]", y = expression(kappa)) +
  theme_minimal() +
  theme(panel.background = NULL,
        plot.background = element_rect(fill = "white", color = "white"),
        axis.text = element_text(color = "black"),
        text = element_text(family = "serif", size = 10, color = "black"))

# IC GREGOS -------------------------------------------

# Definição do tamanho dos intervalos
m <- 5
intervalos <- seq(from = min(df.ic.param$n.serie), to = max(df.ic.param$n.serie) + m, by = m)

# Intervalo de confiança a partir dos percentis
df.ic_sintetico <- ({
  df.ic.param %>% 
    select(estacao, n.serie, kappa.gml, ic.d.kappa, ic.e.kappa, inv.fisher.kappa) %>% 
    mutate(sd.kappa = sqrt(inv.fisher.kappa),
           intervalo = cut(x = n.serie,
                           breaks = intervalos,
                           right = FALSE,
                           labels = intervalos[-1])) %>% 
    group_by(intervalo) %>% 
    summarise(media.kappa = mean(kappa.gml),
              media.sd = mean(sd.kappa),
              perc.kappa.pos = 100*mean(kappa.gml > 0),  # soma dos kappas positivos
              q.e = quantile(kappa.gml, 0.025, na.rm = TRUE),
              q.d = quantile(kappa.gml, 0.975, na.rm = TRUE)) %>% 
    mutate(intervalo = intervalo %>% as.integer,
           ic1.96.e = media.kappa - qnorm(1 - alfa.ic/2, 0, 1)*media.sd,
           ic1.96.d = media.kappa + qnorm(1 - alfa.ic/2, 0, 1)*media.sd)
})

# Plot gregos
plot.kappa.gregos <- ({
  
  # a)
    (ggplot(data = df.ic_sintetico, aes(x = intervalo*m)) +
       # Função ajustada (linha preta)
       geom_smooth(aes(y = media.kappa, color = "Ajuste κ", linetype = "Ajuste κ"), linewidth = 0.7, se = FALSE, method = "gam") +
       geom_smooth(aes(y = ic1.96.e, color = "Ajuste IC 95%", linetype = "Ajuste IC 95%"), linewidth = 0.45, se = FALSE, method = "gam") +
       geom_smooth(aes(y = ic1.96.d, color = "Ajuste IC 95%", linetype = "Ajuste IC 95%"), linewidth = 0.45, se = FALSE, method = "gam") +
       # Pontos κ (vermelho) e IC 95% (roxo)
       geom_point(aes(y = media.kappa, color = "κ", shape = "κ"), size = 1.8) +
       geom_point(aes(y = ic1.96.e, color = "IC 95%", shape = "IC 95%"), size = 1.5) +
       geom_point(aes(y = ic1.96.d, color = "IC 95%", shape = "IC 95%"), size = 1.5) +
       # Escalas manuais para cores, linetypes e shapes
       scale_color_manual(values = c("Ajuste κ" = "black", "Ajuste IC 95%" = "black", "κ" = "red", "IC 95%" = "purple")) +
       scale_linetype_manual(values = c("Ajuste κ" = "solid", "Ajuste IC 95%" = "dashed", "κ" = NA, "IC 95%" = NA)) +
       scale_shape_manual(values = c("Ajuste κ" = NA, "Ajuste IC 95%" = NA, "κ" = 8, "IC 95%" = 20)) +
       # Customização do tema
       labs(x = "N [anos]", y = "κ", color = "", shape = "", linetype = "", title = "a)") +
       guides(shape = "none", linetype = "none",
              color = guide_legend(override.aes = list(shape = c(NA, NA, 20, 8), linetype = c("dashed", "solid", "blank", "blank")))) +
       theme_minimal() +
       theme(panel.background = NULL,
             legend.position = "bottom",
             legend.margin = margin(t = -0.3, unit = "cm"),
             plot.background = element_rect(fill = "white", color = "white"),
             axis.text = element_text(color = "black"),
             text = element_text(family = "serif", size = 10, color = "black"))
     ) + plot_spacer() +
    # b)
    (ggplot(data = df.ic_sintetico, aes(x = intervalo*m)) +
       # Ajuste do desvio padrão de kappa
       geom_smooth(aes(y = media.sd, color = "Ajuste", linetype = "Ajuste"), linewidth = 0.7, se = FALSE, method = "loess") +
       # Pontos média desvio padrão de kappa
       geom_point(aes(y = media.sd, color = "Desvio-padrão κ", shape = "Desvio-padrão κ"), size = 1.8) +
       labs(x = "N [anos]", y = expression(sigma[kappa]), color = "", shape = "", linetype = "", title = "b)") +
       ylim(c(0, 0.14)) +
       scale_color_manual(values = c("Ajuste" = "black", "Desvio-padrão κ" = "red")) +
       scale_shape_manual(values = c("Ajuste" = NA, "Desvio-padrão κ" = 8)) +
       guides(shape = "none", linetype = "none",
              color = guide_legend(override.aes = list(shape = c(NA, 8), linetype = c("solid", "blank")))) +
       theme_minimal() +
       theme(panel.background = NULL,
             legend.position = "bottom",
             legend.margin = margin(t = -0.3, unit = "cm"),
             plot.background = element_rect(fill = "white", color = "white"),
             axis.text = element_text(color = "black"),
             text = element_text(family = "serif", size = 10, color = "black"))
     ) + plot_spacer() +
    # c)
    (ggplot(data = df.ic_sintetico, aes(x = intervalo*m)) +
       # Ajuste do número de estações com kappa positivo
       geom_smooth(aes(y = perc.kappa.pos, color = "Ajuste", linetype = "Ajuste"), linewidth = 0.7, se = FALSE, method = "gam") +
       # Pontos número de estações com kappa positivo
       geom_point(aes(y = perc.kappa.pos, color = "% κ > 0", shape = "% κ > 0"), size = 1.8) +
       labs(x = "N [anos]", y = "Estações com κ > 0 [%]", title = "c)", color = "") +
       scale_color_manual(
         values = c("Ajuste" = "black", "% κ > 0" = "red"),
         breaks = c("Ajuste", "% κ > 0")) +  # Especifica a ordem na legenda
       scale_shape_manual(
         values = c("% κ > 0" = 8, "Ajuste" = NA),
         breaks = c("Ajuste", "% κ > 0"))+  # Especifica a ordem na legenda
       guides(shape = "none", linetype = "none",
              color = guide_legend(override.aes = list(shape = c(NA, 8), linetype = c("solid", "blank")))) +
       theme_minimal() +
       theme_minimal() +
       theme(panel.background = NULL,
             legend.position = "bottom",
             legend.margin = margin(t = -0.3, unit = "cm"),
             plot.background = element_rect(fill = "white", color = "white"),
             axis.text = element_text(color = "black"),
             text = element_text(family = "serif", size = 10, color = "black"))
     ) +
    
    plot_layout(ncol = 1, heights = c(1, -0.1, 1, -0.1, 1))
    
}); plot.kappa.gregos

ggsave(filename = "Plotagens/Incerteza/IC e Kappa Médios Gregos.png",
       plot = plot.kappa.gregos,
       width = 9.5, height = 22, units = "cm",
       dpi = 300)

# GRÁFICO IC QUANTIS ------------------------------------------------------

# Largura relativa dos intervalos de confiança
df.ic.quantis <- df.ic.quantis %>% 
  mutate(l.r = (ic.d - ic.e)/gml.gev)

# Boxplot largura relativa dos intervalos de confiança
plot.boxplot.lr <- ({
  df.ic.quantis %>% 
    # filter(!between(kappa.gml, -0.001, 0.001)) %>%
    mutate(intervalo = cut(n.serie,
                           breaks = c(15, 25, 35, 45, 55, Inf),
                           right = FALSE,
                           labels = c("15-25", "25-35", "35-45", "45-55", "≥55"))) %>% 
    ggplot() +
    facet_wrap(~tempo.retorno, labeller = labeller(tempo.retorno = function(x) paste(x, "anos"))) +
    stat_boxplot(aes(x = intervalo, y = l.r), geom ='errorbar', width = 0.6) +
    geom_jitter(aes(x = intervalo, y = l.r), alpha = 0.2, size = 0.3, color = "grey70") +
    geom_boxplot(aes(x = intervalo, y = l.r), alpha = 0.6, linewidth = 0.4, fill = "#a086d2",
                 outlier.shape = 4) +
    ylim(c(0, 1.6)) +
    labs(x = "N [anos]", y = "Largura relativa do IC 95%") +
    theme_minimal() +
    theme(panel.background = NULL,
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
}); plot.boxplot.lr

ggsave(filename = "Plotagens/Incerteza/Boxplot Largura Relativa IC por TR.png",
       plot = plot.boxplot.lr,
       width = 16, height = 16, units = "cm",
       dpi = 300)

# IC Estações 2549006 e 49009
df.ic.49009.2549006 <- ({
  rbind(
    # Pontos observados 2549006
    df.precip.max.anual %>% 
      filter(Estacao_codigo == 2549006) %>% 
      select(Estacao_codigo, Pdmax) %>%
      rename("estacao" = "Estacao_codigo", "pd.max" = "Pdmax") %>% 
      arrange(desc(pd.max)) %>% 
      mutate(n = length(estacao), i = row_number(), p.weibull = i/(n +1), Tr = 1/p.weibull),
    # Pontos observados 49009
    df.precip.max.anual %>% 
      filter(Estacao_codigo == 49009) %>% 
      select(Estacao_codigo, Pdmax) %>%
      rename("estacao" = "Estacao_codigo", "pd.max" = "Pdmax") %>% 
      arrange(desc(pd.max)) %>% 
      mutate(n = length(estacao), i = row_number(), p.weibull = i/(n +1), Tr = 1/p.weibull)
    )})

# Quantis e IC ajustados
df.q.ic.ajustados.409009.2549006 <- 
  # Quantis e IC calculados 2549006
  rbind(
    data.frame(estacao = rep(2549006, 11)) %>% 
    mutate(csi = rep(df.quantis$csi.gml.gev[df.quantis$estacao == 2549006], 11),     # pega csi de df.quantis
           alpha = rep(df.quantis$alpha.gml.gev[df.quantis$estacao == 2549006], 11), # pega alpha de df.quantis
           kappa = rep(df.quantis$kappa.gml[df.quantis$estacao == 2549006], 11),     # pega kappa de df.quantis
           hess.gl = rep(df.quantis$hess.gl[df.quantis$estacao == 2549006], 11),
           tempo.retorno = c(1.5, 2, 3, 5, 7, 10, 20, 50, 100, 200, 500),
           p.0 = rep(df.quantis$p.0[df.quantis$estacao == 2549006], 11),
           p = 1 - 1/tempo.retorno,
           p.ajustado = (p - p.0)/(1 - p.0)) %>%
    rowwise() %>%
    mutate(q = fun.q.gev(p = p.ajustado, param = c(csi, alpha, kappa))) %>%
    ungroup() %>%
    mutate(hess.gl = lapply(hess.gl, function(vetor) matrix(vetor, ncol = 3, nrow = 3, byrow = TRUE)),
           cov.gl = lapply(hess.gl, solve),
           deriv.csi = 1,
           deriv.alpha = fun.deriv.alpha(p = p.ajustado, kappa = kappa),
           deriv.kappa = fun.deriv.kappa(p = p.ajustado, alpha = alpha, kappa = kappa),
           grad.q = purrr::pmap(list(deriv.csi, deriv.alpha, deriv.kappa), c),
           sd.q = purrr::pmap(list(grad.q, cov.gl), function(grad_q, cov_gl){t(grad_q) %*% cov_gl %*% grad_q}) %>% as.numeric %>% sqrt,
           ic.e = q - qnorm(1 - alfa.ic/2, 0, 1)*sd.q,
           ic.d = q + qnorm(1 - alfa.ic/2, 0, 1)*sd.q),
    data.frame(estacao = rep(49009, 11)) %>% 
      mutate(csi = rep(df.quantis$csi.gml.gev[df.quantis$estacao == 49009], 11),     # pega csi de df.quantis
             alpha = rep(df.quantis$alpha.gml.gev[df.quantis$estacao == 49009], 11), # pega alpha de df.quantis
             kappa = rep(df.quantis$kappa.gml[df.quantis$estacao == 49009], 11),     # pega kappa de df.quantis
             hess.gl = rep(df.quantis$hess.gl[df.quantis$estacao == 49009], 11),
             tempo.retorno = c(1.5, 2, 3, 5, 7, 10, 20, 50, 100, 200, 500),
             p.0 = rep(df.quantis$p.0[df.quantis$estacao == 49009], 11),
             p = 1 - 1/tempo.retorno,
             p.ajustado = (p - p.0)/(1 - p.0)) %>%
      rowwise() %>%
      mutate(q = fun.q.gev(p = p.ajustado, param = c(csi, alpha, kappa))) %>%
      ungroup() %>%
      mutate(hess.gl = lapply(hess.gl, function(vetor) matrix(vetor, ncol = 3, nrow = 3, byrow = TRUE)),
             cov.gl = lapply(hess.gl, solve),
             deriv.csi = 1,
             deriv.alpha = fun.deriv.alpha(p = p.ajustado, kappa = kappa),
             deriv.kappa = fun.deriv.kappa(p = p.ajustado, alpha = alpha, kappa = kappa),
             grad.q = purrr::pmap(list(deriv.csi, deriv.alpha, deriv.kappa), c),
             sd.q = purrr::pmap(list(grad.q, cov.gl), function(grad_q, cov_gl){t(grad_q) %*% cov_gl %*% grad_q}) %>% as.numeric %>% sqrt,
             ic.e = q - qnorm(1 - alfa.ic/2, 0, 1)*sd.q,
             ic.d = q + qnorm(1 - alfa.ic/2, 0, 1)*sd.q))

plot.ic.quantil <- ({
  df.ic.quantis %>% 
    filter(estacao %in% c(2549006, 49009)) %>% # filtrar somente estações 49009 (18 anos) e 2549006 (112 anos)
    ggplot(aes(x = tempo.retorno)) +
    facet_wrap(~estacao, labeller = labeller(estacao = function(x) paste("Estação", x))) +
    # Valor do quantil
    geom_line(aes(y = gml.gev, linetype = "Quantil"), linewidth = 0.6) +
    # Precipitações observadas
    geom_point(data = df.ic.49009.2549006, aes(x = Tr, y = pd.max, color = "Observado", shape = "Observado"), size = 1.5) +
    # IC
    geom_line(aes(y = ic.e, linetype = "IC 95%"), linewidth = 0.4) +
    geom_line(aes(y = ic.d, linetype = "IC 95%"), linewidth = 0.4) +
    scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500),
                  labels = c("10", "20", "50", "100", "200", "500")) +
    scale_linetype_manual(values = c("Quantil" = "solid", "IC 95%" = "dashed")) +
    scale_color_manual(values = c("Observado" = "red")) +
    scale_shape_manual(values = c("Observado" = 8)) +
    labs(x = "TR [anos]", y = "Precipitação Diária Máxima Anual [mm]", linetype = "", color = "", shape = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
}); plot.ic.quantil

plot.ic.quantil2 <- ({
  df.q.ic.ajustados.409009.2549006 %>% 
    ggplot(aes(x = tempo.retorno)) +
    facet_wrap(~estacao, labeller = labeller(estacao = function(x) paste("Estação", x))) +
    # Valor do quantil
    geom_line(aes(y = q, linetype = "Quantil"), linewidth = 0.6) +
    # Precipitações observadas
    geom_point(data = df.ic.49009.2549006, aes(x = Tr, y = pd.max, color = "Observado", shape = "Observado"), size = 1.5) +
    # IC
    geom_line(aes(y = ic.e, linetype = "IC 95%"), linewidth = 0.4) +
    geom_line(aes(y = ic.d, linetype = "IC 95%"), linewidth = 0.4) +
    scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500),
                  labels = c("10", "20", "50", "100", "200", "500")) +
    scale_linetype_manual(values = c("Quantil" = "solid", "IC 95%" = "dashed")) +
    scale_color_manual(values = c("Observado" = "red")) +
    scale_shape_manual(values = c("Observado" = 8)) +
    labs(x = "TR [anos]", y = "Precipitação Diária Máxima Anual [mm]", linetype = "", color = "", shape = "") +
    theme_minimal() +
    theme(panel.background = NULL,
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 11, color = "black"),
          strip.text = element_text(size = 11, hjust = 0, family = "serif", color = "black"))
}); plot.ic.quantil2

ggsave(filename = "Plotagens/Incerteza/Estações 49009 e 2549006 2.png",
       plot = plot.ic.quantil2,
       width = 18, height = 12, units = "cm",
       dpi = 300)