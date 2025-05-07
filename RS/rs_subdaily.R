
# PACKAGES ----------------------------------------------------------------

BiocManager::install("rhdf5")
pacman::p_load(pacman, dplyr,
               hdf5r,       # pacote para manipulação de arquivos HDF5
               BiocManager, # test with Bioconductor rhdf5 package
               rhdf5)       # pacote p/ manipulação de arquio HDF5


# READ DATA ---------------------------------------------------------------

path <- "C:/Users/tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/Fonte de Dados/Dados_RS/Dados subdiarios consolidados"

# Read files with ".h5" file extension
files <- list.files(path = path, pattern = "*.h5")

# Isolate a single file (2001)
files.2001 <- files[1]
complete.path <- paste0(path, "/", files.2001)

# Test with rhdf5 package, functions are not able to read VLEN files and turns them into NA
df.info <- rhdf5::h5read(paste0(path, "/", files.2001), "table_info")
df.data <- rhdf5::h5read(paste0(path, "/", files.2001), "table_data")

# Test with hdf5r package
h5file <- H5File$new(complete.path, mode = "r")
df.info <- h5file[["table_info"]]
df.data <- h5file[["table_data"]]
