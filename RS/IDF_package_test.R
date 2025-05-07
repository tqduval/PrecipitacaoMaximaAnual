
# PACOTES -----------------------------------------------------------------

install.packages("IDF")
pacman::p_load(pacman, tidyverse, IDF)


# TESTE DA FUNÇÃO ---------------------------------------------------------

# 0. Criar dados amostrais sintéticos de precipitação horária
set.seed(999)

dates <- seq(as.POSIXct("2000-01-01 00:00:00"), # série de datas por hora
             as.POSIXct("2019-12-31 23:00:00"),
             by = "hour")

sample_rainfall <- rgamma(n = length(dates),    # série sintética com a distribuição gamma
                          shape = 0.05,
                          rate = 0.4)

df_precip <- data.frame(date = dates,             # dataframe c/ série sintética
                        rain_h = sample_rainfall) # por padrão o pacote escolhe o nome RR

# 1. Acumular séries e extrair máximos
durations <- 2^(0:6) # durações acumuladas [h]
annual_max <- IDF::IDF.agg(data = list(df_precip),      # lista séries temporais para cada estação (nesse caso 1)
                           ds = durations,              # durações p/ acumular
                           na.accept = 0.1,             # porcentagem de falhas
                           names = c("date", "rain_h")) # nome das colunas

plot(annual_max$ds, annual_max$xdat, log = "xy", xlab = "Duration [h]", ylab = "Intensity [mm/h]")

# 2. Ajustar a GEV às intensidades anuais máximas
fit <- IDF::gev.d.fit(xdat = annual_max$xdat,
                      ds = annual_max$ds,
                      sigma0link = make.link("log"),
                      eta2_zero = FALSE,
                      tau_zero = FALSE)

gev.d.diag(fit = fit, pch = 1, ci = TRUE) # avaliar o ajuste

gev_par <- gev.d.params(fit) # obter os parâmetros

# 3. Plotar as curvas
plot(annual_max$ds,
     annual_max$xdat,
     log = "xy",
     xlab = "Duration [h]", ylab = "Intensity = [mm/h]")

IDF.plot(durations, gev_par, add = TRUE)


# TESTAR AJUSTES DIFERENTES -----------------------------------------------


