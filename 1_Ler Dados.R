
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)


# IMPORTAR DADOS ----------------------------------------------------------

# Ler dados de precipitação máxima anual (dados brutos ANA organizados pelo código Quarto do Dirceu)
df.precip.completa <- ({
  read.table(file = "Fonte de Dados/df.precip.diaria.txt",
             header = TRUE,
             sep = "\t",
             dec = ".",
             fileEncoding = "UTF-8") %>% 
    select(c(Date, Pdmax, Estacao_codigo, BaciaCodigo, SubBaciaCodigo, Altitude, AreaDrenagem, NomedoRio, Latitude, Longitude))
})

# Salvar como RDS
saveRDS(object = df.precip.completa, file = "Dados Gerados/df.precip.completa.rds")


# OUTLIERS ----------------------------------------------------------------

# Identificar estações que apresentem Pdmax > 700 mm
outlier.limiar <- 700

# Criar um data.frame somente c/ dados acima do limiar definido
df.outlier <- ({
  df.precip.completa %>% 
    select(Estacao_codigo, Date, Pdmax) %>% 
    rename(estacao = Estacao_codigo, data = Date, pd.max = Pdmax) %>% 
    filter(estacao %in% unique(df.precip.completa[which(df.precip.completa$Pdmax > outlier.limiar), 3])) %>% # filtrar estações em que Pdmax > outlier.limiar
    mutate(data = as_date(data, format = "%d/%m/%Y")) %>%
    arrange(estacao, data) %>% 
    group_by(estacao) %>%
    mutate(dif.anos = as.numeric(difftime(data, lag(data), units = "days")) / 365.25) %>%                    # calcular diferença entre anos consecutivos
    ungroup()
})

# Verificações
df.outlier$estacao %>% unique %>% length                           # número de observações acima de 'outlier.limiar'
df.outlier$estacao[df.outlier$dif.anos >= 5] %>% unique %>% length # número de estações com mais de 5 anos consecutivos de falha

# Remover outliers e estações que atendam os critérios acima
df.precip.max.anual <- ({
  df.precip.completa %>% 
    filter(Pdmax < outlier.limiar,                                            # remove observações c/ Pdmax > 700 mm
           !Estacao_codigo %in% df.outlier$estacao[df.outlier$dif.anos >= 5]) # remove as estações c/ Pdmax > 700 e +5 anos sem dados
})