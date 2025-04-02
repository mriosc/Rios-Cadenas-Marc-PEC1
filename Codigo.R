###
### Código para crear el objeto SummarizedExperiment
###



# Cargar librerias
library(SummarizedExperiment)
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Ruta del archivo
file_path <- "Gastric_NMR.xlsx"

# Leer hojas del archivo Excel
data_sheet <- read_excel(file_path, sheet = "data")
peak_sheet <- read_excel(file_path, sheet = "peak")

# Mostrar primeras filas
head(data_sheet)
head(peak_sheet)

# Asegurar que la matriz de expresión tenga metabolitos en filas y muestras en columnas
expr_data <- data_sheet %>%
  select(starts_with("M")) %>%
  as.matrix()

rownames(expr_data) <- data_sheet$Sample_id  # Etiquetar filas con IDs de muestra
expr_data <- t(expr_data)  # Transponer para que metabolitos sean filas

# Asegurar que meta_samples tenga filas con nombres de muestra
meta_samples <- data_sheet %>%
  select(Sample_id, Sample_Type, Batch) 

rownames(meta_samples) <- meta_samples$Sample_id  # Asignar Sample_id como rownames
meta_samples$Sample_id <- NULL  # Eliminar columna redundante


# Asegurar que meta_metabolites tenga filas con nombres de metabolitos
meta_metabolites <- peak_sheet %>%
  rename(Metabolite = Label) %>%
  select(Name, Metabolite)

rownames(meta_metabolites) <- meta_metabolites$Name  # Asignar Name como rownames
meta_metabolites$Name <- NULL  # Eliminar la columna redundante


# Crear objeto SummarizedExperiment con dimensiones correctas
se <- SummarizedExperiment(
  assays = list(counts = expr_data),
  colData = meta_samples,
  rowData = meta_metabolites
)

# Mostrar estructura del objeto
se

# Guardar el objeto SummarizedExperiment en formato .Rda
save(se, file = "SummarizedExperiment_Gastric_NMR.Rda")



###
### Código para el análisis exploratorio
###



library(tidyr)

df_long <- data_sheet %>%
  select(Sample_Type, starts_with("M")) %>%
  pivot_longer(cols = starts_with("M"), names_to = "Metabolite", values_to = "Intensity")

## Boxplot de distribución de intensidad
ggplot(df_long, aes(x = Sample_Type, y = Intensity, fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Distribución de Intensidad de Metabolitos por Tipo de Muestra",
       x = "Tipo de Muestra", y = "Intensidad (log10)")

# Verificación de valores no finitos
cat("Número de valores NA en la matriz de expresión:", sum(is.na(assay(se))), "\n")
cat("Número de valores Inf o NaN en la matriz de expresión:", sum(!is.finite(assay(se))), "\n")


##
## Prueba de normalidad
##

set.seed(123)  

# Obtener el número total de observaciones en cada grupo
n_qc <- sum(df_long$Sample_Type == "QC")
n_sample <- sum(df_long$Sample_Type == "Sample")

# Seleccionar la cantidad mínima entre 5000 y el tamaño real
sample_qc <- sample(df_long$Intensity[df_long$Sample_Type == "QC"], min(n_qc, 5000))
sample_sample <- sample(df_long$Intensity[df_long$Sample_Type == "Sample"], min(n_sample, 5000))


##
##Prueba de Shapiro-Wilk
##

shapiro.test(sample_qc)

shapiro.test(sample_sample)

ggplot(df_long, aes(sample = Intensity, color = Sample_Type)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Sample_Type) +
  theme_minimal()


##
## Prueba de Mann-Whitney U (para datos no normales)
##

wilcox.test(Intensity ~ Sample_Type, data = df_long)

ggplot(df_long, aes(x = Sample_Type, y = Intensity, fill = Sample_Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "Comparación de Intensidades entre Grupos",
       subtitle = paste("Wilcoxon p-value:", signif(0.00285, 3))) +
  theme_minimal()
r_wilcoxon <- abs(qnorm(0.00285) / sqrt(nrow(df_long)))
r_wilcoxon


##
## Pheatmap
##

library(pheatmap)

# Tratamiento de valores faltantes
expr_matrix <- assay(se)
expr_matrix <- apply(expr_matrix, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x[!is.finite(x)] <- median(x, na.rm = TRUE)
  return(x)
})
expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"

# Confirmaciones
cat("Número de valores NA después de limpieza:", sum(is.na(expr_matrix)), "\n")

cat("Número de valores Inf o NaN después de limpieza:", sum(!is.finite(expr_matrix)), "\n")

pheatmap(
  t(expr_matrix), 
  annotation_col = as.data.frame(colData(se)),  # Convertir colData a data.frame
  main = "Mapa de Calor de Metabolitos",
  cluster_cols = TRUE,
  cluster_rows = TRUE
)