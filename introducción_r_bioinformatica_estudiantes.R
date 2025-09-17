# Curso introducción al análisis de datos con R - BIOINFORMÁTICA
# Sesión 1: Fundamentos de R aplicados a datos biológicos

# Autor: dialvas-anahuac Córdoba
# Versión: 2025

#### INSTRUCCIONES: 
#### Este script tiene ejercicios básicos. Para que recuerden que hacer cada función o apartado recuerden agregar sus propios comentarios con la tecla # 
#### O pueden hacer sus notas en la hoja de trucos que ha entregado la docente

#### Empecemos: 

# Determinar directorio de trabajo
getwd()
# 


##=========================##
# 1. calculadora: 
##=========================##

# Cálculos básicos
50 + 100  

# divisiones
2/10000  

# Funciones matemáticas
log(1)   
exp(2)   
sqrt(16) 

##============##
# 2. Comparaciones: 
##============##

pH_muestra1 <- 7.4
pH_muestra2 <- 6.8
pH_muestra1 == pH_muestra2  # ¿Son iguales los pH?
pH_muestra1 != pH_muestra2  # ¿Son diferentes?
pH_muestra1 > pH_muestra2   # ¿Cuál es más básico?
pH_muestra1 >= 7.0          # ¿Es neutro o básico?


##=========================##
# 3. Asignar variables 
##=========================##

concentracion_DNA <- 250.5  # ng/µL
numero_muestras <- 48
temperatura_incubacion <- 37.0
especie_estudiada <- "Escherichia coli"

# Mostrar las variables
concentracion_DNA
numero_muestras
especie_estudiada

# Cálculos con variables
volumen_total <- concentracion_DNA * 0.002  
volumen_total


##================##
# 4. Administrar entorno
##================##

ls()  # Listar todas las variables
rm(concentracion_DNA)  # Eliminar una variable específica
# rm(list = ls())  # Borrar todo (¡cuidado!)


##========================##
# 5. Exploración de paquetes para bioinformática
##========================##

# Paquetes comunes en bioinformática (no ejecutar install si ya están instalados)
# install.packages("seqinr")      # Análisis de secuencias
# install.packages("stringr")     # Manipulación de texto/secuencias
# install.packages("dplyr")       # Manipulación de datos
# install.packages("ggplot2")     # Gráficos

library(stringr)  # Para trabajar con secuencias como texto

##==========##
# 6. Ayuda
##==========##

help(paste)
?paste

##==========##
# EJERCICIO 1 -
##==========##

# ¿Cuál será el valor de cada variable?
peso_molecular <- 150.5   # kDa
tiempo_experimento <- 120  # minutos
peso_molecular <- peso_molecular * 2.3  # Conversión de unidades


# Escriba un comando para comparar peso_molecular con tiempo_experimento


#======================================#
# EJERCICIO: Limpia tu entorno de trabajo. 
#=====================================#
rm(list = ls())

##==================##
# EJERCICIO 2 - Secuencias de ADN
##==================##

# Ejemplos con secuencias de ADN
secuencia1 <- c("ATG", "GCT")
secuencia2 <- "TAA"

paste(secuencia1, secuencia2)                    # Unir con espacio
paste(secuencia1, secuencia2, sep = "")          # Unir sin separador
paste(secuencia1, secuencia2, collapse = "_")     # Colapsar en una sola cadena


# ¿Cuál es la diferencia entre sep y collapse?



#### pero yo deseo un resultado así "ATGGCTTAA" ¿cómo puedo hacerlo?

paste(c(secuencia1, secuencia2), collapse = "")


##=============##
# 7. Tipos de datos 
##=============##


# Vectores - Tienen una dimensión y un tipo de dato. 
expresion_genes <- c(2.5, 1.8, 0.9, 3.2, 1.1, 2.7)  # Niveles de expresión
nombres_genes <- c("ACTB", "GAPDH", "TP53", "BRCA1", "MYC", "EGFR")
expresado <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)   # ¿Está expresado?

# Secuencias numéricas útiles
posiciones <- 1:10           # Posiciones en una secuencia
potencias_2 <- 2^(1:5)       # Diluciones seriadas

expresion_genes
nombres_genes
str(expresion_genes)


# Matrices - Datos de expresión génica
matriz_expresion <- matrix(1:12, nrow = 3, ncol = 4)
matriz_expresion

colnames(matriz_expresion) <- c("Muestra1", "Muestra2", "Muestra3", "Muestra4")
rownames(matriz_expresion) <- c("Gen_A", "Gen_B", "Gen_C")
matriz_expresion

# Explorar la matriz
class(matriz_expresion)
str(matriz_expresion)
dim(matriz_expresion)    # dimensiones
nrow(matriz_expresion)   # número de filas (genes)
ncol(matriz_expresion)   # número de columnas (muestras)

# ¿Cuántos valores total tiene la matriz?
length(matriz_expresion)

# Listas - Datos heterogéneos de un experimento
secuencia_fibonacci <- c(0,1,1,2,3,5,8,13)  # Patrón en la naturaleza
especies_plantas <- c("Rosa", "Tulipán", "Girasol")
datos_crecimiento <- matrix(1:12, ncol = 3)

experimento_completo <- list(secuencia_fibonacci, especies_plantas, datos_crecimiento)
experimento_completo

##===================##
# 8. Data.frame - Tabla de datos biológicos
##===================##

#### primero limpiar espacio de trabajo
rm(list = ls())

# Usar dataset iris (clásico en biología)
iris  # Dataset de flores incorporado en R

# Crear nuestro propio dataframe biológico
datos_bacterias <- data.frame(
  cepa = c("E.coli_K12", "E.coli_DH5α", "B.subtilis_168", "S.aureus_NCTC"),
  crecimiento_hr = c(0.8, 0.9, 1.2, 1.5),              # Tiempo de duplicación en horas
  temperatura_optima = c(37, 37, 30, 37),               # Temperatura óptima °C
  resistente_ampicilina = c(FALSE, TRUE, FALSE, TRUE)   # Resistencia a antibiótico
)

datos_bacterias

# Guardar el dataframe en archivo CSV
write.csv(datos_bacterias, "Data/datos_bacterias.csv", row.names = FALSE)


# Limpiar entorno y volver a cargar
rm(datos_bacterias)

# Buscar función para cargar archivos CSV (común en bioinformática)
?read.csv
datos_bacterias <- read.csv(" ", stringsAsFactors = TRUE)

# Explorar columnas con el operador $
datos_bacterias$cepa
datos_bacterias$crecimiento_hr
datos_bacterias$temperatura_optima

# ¿Qué sucede si sumamos datos incompatibles?
datos_bacterias$crecimiento_hr + datos_bacterias$cepa 

# Verificar tipos de datos
class(datos_bacterias$cepa)
class(datos_bacterias$crecimiento_hr)
class(datos_bacterias$resistente_ampicilina)

str(datos_bacterias)

##=============##
# Coerción de datos
##=============##

# Mezclar tipos en vectores
vector_mixto <- c('control', TRUE)         # Fuerza todo a carácter
str(vector_mixto)

vector_numerico_logico <- c(0, TRUE)       # TRUE se convierte a 1
str(vector_numerico_logico)

# Convertir concentraciones de texto a número
concentraciones_texto <- c('2.5', '1.8', '3.2')
str(concentraciones_texto)

concentraciones_num <- as.numeric(concentraciones_texto)
concentraciones_num
str(concentraciones_num)

# Convertir a valores lógicos (útil para expresión ON/OFF)
expresion_logica <- as.logical(c(1, 0, 1, 0))
expresion_logica

# Ejemplo con nuestros datos de bacterias
datos_bacterias$resistente_ampicilina
class(datos_bacterias$resistente_ampicilina)

##============##
# Factores - Categorías biológicas
##============##

# Los factores son perfectos para categorizar datos biológicos
tipos_celula <- c('neurona', 'hepatocito', 'miocito', 'fibroblasto', 'linfocito')
str(tipos_celula)

# Convertir a factor (categoría)
categorias_celula <- factor(tipos_celula)
class(categorias_celula)
str(categorias_celula)

##==================##
# Filas y columnas de data.frame
##==================##

# Acceder a columnas de diferentes formas
datos_bacterias$cepa      # Con $
datos_bacterias[,1]       # Con índice de columna
datos_bacterias[,"cepa"]  # Con nombre de columna

# Acceder a filas
datos_bacterias[1,]       # Primera fila (primera bacteria)
datos_bacterias[1:2,]     # Primeras dos filas

##========================================##
# EXPRESIONES REGULARES EN BIOINFORMÁTICA
##========================================##

library(stringr)
library(dplyr)

# Texto con información de secuencias y genes
info_genomica <- "El gen TP53 (tumor protein p53) está localizado en el cromosoma 17p13.1 y codifica una proteína de 393 aminoácidos. La secuencia de referencia es NM_000546.5 con una longitud de 1,182 pb. Otros genes importantes incluyen BRCA1 (NM_007294.3, 5,592 pb), EGFR (NM_005228.4, 4,293 pb) y MYC (NM_002467.4, 2,185 pb). La expresión de estos genes se midió mediante qPCR con eficiencias del 98.5%, 97.2%, 99.1% y 96% respectivamente."

info_genomica

# Reemplazar "aminoácidos" por "aa"
str_replace(info_genomica, "aminoácidos", "aa")

# Extraer números (longitudes, eficiencias, etc.)
numeros_genomicos <- str_extract_all(info_genomica, "\\d+,?\\d*,?\\d*\\.?\\d*")
str(numeros_genomicos)


# Extraer porcentajes (eficiencias de qPCR)
eficiencias <- str_extract_all(info_genomica, "\\d+.\\d*%")
eficiencias

# Extraer números de acceso (formato NM_######.#)
numeros_acceso <- str_extract_all(info_genomica, "NM_\\d{6}\\.\\d+")
numeros_acceso

# Extraer nombres de genes (palabras en mayúsculas de 3-5 letras)
nombres_genes_patron <- str_extract_all(info_genomica, "\\b[A-Z]{2,5}\\d*\\b")
nombres_genes_patron


##===============================================##
# EXPRESIONES REGULARES CON DATAFRAMES BIOLÓGICOS
##===============================================##

library(tibble)

# Crear dataframe con información de publicaciones en bioinformática
df_papers_bio <- tibble(
  autor = c(
    "geneticista_molecular",
    "bioinformatico_1",
    "especialista_genomica",
    "analista_proteinas",
    "biologo_computacional",
    "investigador_crispr",
    "experto_rnaseq",
    "cientifico_datos",
    "bioestadistico",
    "genomica_funcional"
  ),
  abstract = c(
    "Nuestro estudio identifica 127 genes diferencialmente expresados en cáncer de mama.",
    "¿Cuáles son las mejores prácticas para el análisis de datos de RNA-seq?",
    "La tecnología CRISPR-Cas9 ha revolucionado la edición génica en los últimos años.",
    "El análisis proteómico reveló 89 proteínas upreguladas y 156 downreguladas.",
    "Desarrollamos un pipeline bioinformático para analizar variantes genéticas raras.",
    "La expresión de genes relacionados con apoptosis aumentó significativamente.",
    "Cómo optimizar los parámetros de alineamiento para secuencias cortas?",
    "Nuestro algoritmo de machine learning predice la función de proteínas desconocidas.",
    "El estudio genome-wide association identificó 23 SNPs asociados con la enfermedad.",
    "La secuenciación de célula única reveló heterogeneidad en la expresión génica."
  )
)

print(df_papers_bio)

# Filtrar abstracts que son preguntas (investigación exploratoria)

papers_con_preguntas <- df_papers_bio %>%
  filter(str_detect(abstract, "¿?.*\\?"))

print(papers_con_preguntas)

# filtro papers sobre expresión génica
temas_expresion <- df_papers_bio %>%
  filter(str_detect(abstract, "[Ee]xpresión|[Ee]xpresados|[Ee]xpresado"))

print(temas_expresion)


##=========================##
# EJERCICIOS PRÁCTICOS FINALES
##=========================##
 ## 1. Limpia tu espacio de trabajo. 

rm(list = ls())

#################################################
##################################################
# EJERCICIO A: Crear tu propio dataset de genes
##################################################
##################################################

# Crea un dataframe con la función data.frame con los siguientes genes  BRCA1, TOX3,CASP8
###El dataframe debe incluir la siguiente información que se puede encontrar en https://www.ensembl.org/index.html

# - Nombre del gen : 
# - Transcripts:
# - Cromosoma: 
# - Descripción: 


### Nota: recuerda que cada variable se observa como columna. 



#################################################
##################################################
# EJERCICIO B: Texto que contiene secuencias
#################################################
##################################################

# Usando expresiones regulares, extrae información específica del siguiente texto:
texto_secuencias <- "La secuencia ATGCGTACGTTTAAA tiene 15 nucleótidos. El contenido GC es del 46.7%. Se encontraron 3 sitios de restricción para EcoRI y 2 para BamHI."
texto_secuencias

# Extrae:

# 1. La secuencia de ADN
### respuesta



# 2. El número de nucleótidos



# 3. El porcentaje GC



##################################################
##################################################
# EJERCICIO C: Filtra tu dataset
##################################################
##################################################

# Usando tu dataframe de genes del ejercicio A, filtra Genes con mas de 10 transcritos

