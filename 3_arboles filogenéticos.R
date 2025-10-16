##### instalar paquetes de bioconductor
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Biostrings")
##BiocManager::install("pwalign")
#BiocManager::install("msa")
#BiocManager::install ("seqinr")
#library(devtools)
#devtools::install_github("YuLab-SMU/ggmsa")

#install.packages("ape")

#### arbol filogenetico
library(rentrez)
library(Biostrings)
library(msa)
library(ggmsa)
library(seqinr)
library(ape)

#### Vamos a estudiar la familia de genes Shroom 
###(este ejemplo tiene 14 genes)
shroom_family <- c("CAA78718" , "X. laevis Apx" ,         "Rana-xShroom1",
                  "NP_597713" , "H. sapiens APXL2" ,     "Humano-hShroom1",
                  "CAA58534" , "H. sapiens APXL",        "Humano-hShroom2",
                  "ABD19518" , "M. musculus Apxl" ,      "Raton-mShroom2",
                  "AAF13269" , "M. musculus ShroomL" ,   "Raton-mShroom3a",
                  "AAF13270" , "M. musculus ShroomS" ,   "Raton-mShroom3b",
                  "NP_065910", "H. sapiens Shroom" ,     "Humano-hShroom3",
                  "ABD59319" , "X. laevis Shroom-like",  "Rana-xShroom3",
                  "NP_065768", "H. sapiens KIAA1202" ,   "Humano-hShroom4a",
                  "AAK95579" , "H. sapiens SHAP-A" ,     "Humano-hShroom4b",
                  "ABA81834" , "D. melanogaster Shroom", "Mosca-dmShroom",
                  "EAA12598" , "A. gambiae Shroom",      "Mosquito-agShroom",
                  "XP_392427" , "A. mellifera Shroom" ,  "Abeja-amShroom",
                  "XP_783573" , "S. purpuratus Shroom" , "Erizo_de_mar-spShroom") #sea urchin
shroom_family


##### ahora convertimos este vector a matriz

shroom_family_matrix <- matrix(shroom_family,
                              byrow = T,
                              nrow = 14)
shroom_family_matrix

#### y ahora convertimos la matriz a dataframe
shroom_family_df <- data.frame(shroom_family_matrix, 
                           stringsAsFactors = F)

### limpio espacio de trabajo
rm(shroom_family_matrix)
rm(shroom_family)

### asignamos nombres a las columnas 

names(shroom_family_df) <- c("accession", "name.orig","name.new")
shroom_family_df


#### ahora vamos a descargar multiples secuencias de las proteinas previas
shroom_family_df$accession

### vamos a crear una funcion nueva para 
###descargar los datos de cada referencia de acceso
entrez_fetch_list <- function(db, id, rettype){
  
  #setup list for storing output
  n.seq <- length(id)
  list.output <- as.list(rep(NA, n.seq))
  names(list.output) <- id
  
  # get output
  for(i in 1:length(id)){
    list.output[[i]] <- rentrez::entrez_fetch(db = db,
                                              id = id[i],
                                              rettype = rettype)
  }
  
  
  return(list.output)
}

entrez_fetch_list


####### Ahora corremos la funcion para que haga 
####### la descarga de codigos fasta de manera ciclica y los guarde en una lista
####### Este paso puede tardar un poco.

shrooms_lista <- entrez_fetch_list(db = "protein", 
                                  id = shroom_family_df$accession, 
                                  rettype = "fasta")

##### ¿como observo algun elemento de la lista de manera individual (secuencia FASTA)?
shrooms_lista[[14]]

### Ahora debemos limpiar las secuencias porque tienen "\n"
##### para esto creamos la funcion fasta_cleaner
fasta_cleaner <- function(fasta_object, parse = TRUE){
  
  fasta_object <- sub("^(>)(.*?)(\\n)(.*)(\\n\\n)","\\4",fasta_object)
  fasta_object <- gsub("\n", "", fasta_object)
  
  if(parse == TRUE){
    fasta_object <- stringr::str_split(fasta_object,
                                       pattern = "",
                                       simplify = FALSE)
  }
  
  return(fasta_object[[1]])
}

fasta_cleaner
#####Debo hacer la misma limpieza de la sencuencia fasta pero como son 14 elementos FASTA 
#### No es muy util correr 14 veces el mismo codigo y para esto haremos un ciclo for 
##y corremos la funcion fasta_cleaner que creamos previamente


shrooms_lista[[14]]


###### Para extraer las secuencia de la 
###lista tenemos que crear un vector vacio
shrooms_vector <- rep(NA, length(shrooms_lista))
shrooms_vector


###


shrooms_vector

### Necesitamos poner nombre a los elementos del vector 
###(esto para saber a que proteina corresponde cada secuencia)

names(shrooms_vector) <- shroom_family_df$name.new
### limpio espacio de trabajo
rm(shroom_family_df)

shrooms_vector

###### ¡¡¡Muy bien!!!
##### Ya tenemos multiples secuencias. 
#### ahora debemos convertirlas a un objeto StringSet
### que es un tipo de objeto en R 
## para manejar multiples secuencias de forma eficiente. 
#### para esto corremos esta linea de comando

shrooms_vector_ss <- Biostrings::AAStringSet(shrooms_vector)
shrooms_vector_ss
##### si ya tenemos el elemento, 
####  ahora vamos a hacer el alineamiento
###   y construir el arbol filogenetico.

#### La funcion que utlizaremos sera el 
#### MSA (Multiple Sequence Alignment)
### Esta funcion se basa en el algoritmo ClustalW con R

shrooms_align <- msa(shrooms_vector_ss,
                     method = "ClustalW")

shrooms_align

### analicemos el tipo de objeto que se genera
class(shrooms_align)

#### ahora Se utilizara la funcion msaConvert() 
#### para ajustar el objeto msa y hacerlo compatible con
#### funciones del paquete seqinr. Ademas, 
#### se cambiara el nombre del objeto de shrooms_align 
####   shrooms_align_seqinr.

shrooms_align_seqinr <- msaConvert(shrooms_align, 
                                   type = "seqinr::alignment")

shrooms_align_seqinr

#### noten que a pesar de esto la visualizacion 
###  sigue siendo compleja

###########################################################
##### visualizacon de los alineamientos###################
###########################################################
### primero debemos tenen un objeto del tipo 
### AAMultipleAlignment
### para esto corremos la siguinte funcion 

class(shrooms_align) <- "AAMultipleAlignment"

### a continuacion graficamos el alineamiento y lo exportamos a pdf




###Si bien un MSA es una buena manera de examinar una secuencia, 
###es dificil evaluar visualmente toda la informacion. 
##Un arbol filogenetico permite resumir patrones en un MSA. 
##La forma mas rapida de crear arboles filogeneticos es 
##resumir primero un MSA utilizando una matriz de distancia genetica. 
##Cuantos mas aminoacidos sean identicos entre si, 
##menor sera la distancia genetica entre ellos y menor sera la evolucion.
##Normalmente trabajamos en terminos de diferencia o 
##distancia genetica (tambien conocida como distancia evolutiva), 
##aunque a menudo tambien hablamos de similitud 
#o identidad.


#### vamos a calcular distancia genetica

shrooms_dist <- seqinr::dist.alignment(shrooms_align_seqinr, 
                                              matrix = "identity")
class(shrooms_dist)

##### por fin el arbol filogenetico
### metodos NJ
tree <- nj(shrooms_dist)

#### plot tree


#source: https://brouwern.github.io/lbrb/worked-example-building-a-phylogeny-in-r.html
