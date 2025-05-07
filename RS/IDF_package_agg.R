
# FUNÇÃO DE AGREGAÇÃO DE DURAÇÕES -----------------------------------------

# Função conforme https://gitlab.met.fu-berlin.de/Rpackages/idf_package/-/blob/public/R/IDF.R?ref_type=heads

IDF.agg <- function(data,
                    ds,
                    na.accept = 0,
                    which.stations = NULL,
                    which.mon = list(0:11),
                    names = c("date", "RR"),
                    cl = 1){
  
  # Confere se o argumento 'data' é uma lista
  if(!inherits(data, "list")){
    stop("Argument 'data' must be a list, insted it is a: ", class(data))
    }
  
  # Função 2: análises preliminares de consistência e encontrar a resolução temporal da estação 'dtime'
  # Essa função calcula qual a resolução temporal que a estação apresenta baseado na diferença entre a data das duas primeiras observações
  # depois confere se existe alguma observação nessa estação que a resolução sejea diferente e para o progresso da função caso encontre alguma
  # que a diferença seja maior que a tolerância estabelecida e roda ao final a função de agregação
  agg.station <- function(station){
    
    # Extrai o data.frame para uma dada estação da lista
    data.s <- data[[station]]
    
    # Confere se a lista contém data.frames
    if(!is.data.frame(data.s)){
      stop("Elements of 'data' must be data.frames. But element ", station, " contains: ", class(data.s))
    }
    
    # Usa is.element() parar conferir se os nomes do argumento 'names' são os mesmos nomes das colunas dos data.frames
    if(sum(is.element(names[1:2], names(data.s))) != 2){
      stop("Data.frame of station ", station, " does not contain $", names[1], " or $", names[2], ".")
    }
    
    # Calcula a resolução temporal da estação baseado na diferença de tempo entre as duas primeiras observações do data.frame
    dtime <- as.numeric(x = (data.s[, names[1]][2] - data.s[, names[1]][1]),
                        units = "hours")
    
    # Confere as durações ewm 'ds' são múltiplas da resolução temporal da estação
    if(any((ds/dtime) %% 1 > 1e-8)){
      stop("At least one of the given aggregation durations is not multiple of the time resolution = ", dtime, " hours at station ", station, ".")
    }
    
    # Função 1: agregar por durações e encontrar o máximo anual de cada uma
    # Usa o pacote RcppRoll para implementar funções do C++ de janela móvel (roll_sum e roll_mean) mais rápido que alternativas do R
    agg.ts <- function(ds){
      
      # Criar um vetor de somas móveis da coluna de precipitação, usando uma janela de ds/dtime
      runsum <- RcppRoll::roll_sum(x = data.s[, names[2]], # vetor com dados da coluna 2 (precipitação)
                                   n = round(ds/dtime),    # tamanho da janela móvel (ds é o vetor c/ as durações e dtime a resolução da estação)
                                   fill = NA,              # preenche as pontas com NA
                                   align = "right")        # alinha a janela à direita (no índice i)
      
      # Converte a soma acumulada em intensidade por hora
      runsum <- runsum/ds
      
      # Analisa cada data.frame dentro da lista 'data' (por isso o lapply) 
      max.subset <- 
        lapply(1:length(which.mon), 
               function(m.i){
                 
                 # Retornar um vetor lógico que marca como TRUE quais os meses que estão dentro do argumento 'which.mon'
                 subset <- is.element(as.POSIXlt(data.s[, names[1]])$mon, which.mon[[m.i]])
                 
                 # Agrupa as observações por ano e compara o percentual de falhas com o na.accept
                 # ctapply() é uma versão otimizada da função tapply() que faz comparações baseadas em uma categoria (nesse caso o ano)
                 max <- fastmatch::ctapply(X = runsum[subset],                                           # somas móveis de cada mês analisado
                                           INDEX = (as.POSIXlt(data.s[, names[1]])$year + 1900)[subset], # ano de cada mês com subset = TRUE
                                           FUN = function(vec){                                          # máximo anual se n.falhas <= na.accept
                                             
                                             # Número de NAs no ano agrupado
                                             n.na <- sum(is.na(vec))
                                             
                                             # Retorna o máximo anual caso a função tenha número de falhas aceitável
                                             max <- ifelse(n.na <= na.accept*length(vec),
                                                           max(vec, na.rm = TRUE),
                                                           NA)
                                             
                                             return(max)
                                             
                                           })
                 
                 # Data.frame resultante de máximos
                 df <- data.frame(xdat = max,                      # intensidades máximas anuais
                                  ds = ds,                         # durações (cte em todo o conjunto)
                                  year = as.numeric((names(max))), # converte os nomes de max (contém strings dos anos) em numérico
                                  mon = deparse(which.mon[[m.i]]), # sequencia de meses analisados
                                  station = station,               # nome da estação que está sendo analisada
                                  stringsAsFactors = FALSE)
                 
                 # Investigar: guardar a data em que ocorreu o máximo
                 
                 return(df)
                 
               })
      
      df <- do.call(rbind, max.subset)
      
      return(df) # máximos para uma duração
      
    }
    
    # Chamar a função agg.ts dentro do parLapply (processamento paralelo) para calcular os máximos para todas as durações
    if(cl > 1){ # Aplicar processamento paralelo
      
      clust <- parallel::makeCluster(cl, type = "PSOCK")
      data.agg <- parallel::parLapply(cl = clust,
                                      X = ds,
                                      fun = agg.ts)
      parallel::stopCluster(clust)
      
    } else{     # Não aplicar processamento paralelo
      data.agg <- lapply(X = ds, agg.ts)
    }
    
    # Data.frame com máximos anuais para todas as durações em 'ds'
    df <- do.call(rbind, data.agg)
    
    return(df)
    
  }
  
  # Define quais estações utilizar
  if(is.null(which.stations))which.stations <-
    if(is.null(names(data))){
      1:length(data)
    } else{
        names(data)
    }
  
  # Chama a função agg.station com lapply p/ agregar as durações para todas as estações
  station.list <- pbapply::pblapply(which.stations, agg.station)
  
  return(do.call('rbind', station.list))
  
}


# TESTAR FUNÇÃO -----------------------------------------------------------

library(tidyverse)

# Uma estação qualquer
subset.table.data1 <- rs_data_table[rs_data_table$gauge_code == "87399000",]
subset.table.data1 <- list("87399000" = data.frame(datetime = subset.table.data1$datetime,
                                                   rain_mm = subset.table.data1$rain_mm))

# Aplicar função
subset.table.data1.agg <- IDF.agg(data = subset.table.data1,
                                  ds = 2^(0:4),
                                  na.accept = 0,
                                  which.stations = NULL,
                                  which.mon = list(0:11),
                                  names = c("datetime", "rain_mm"),
                                  cl = 2)
