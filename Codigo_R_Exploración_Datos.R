
##########################################
## Obtención y preparación de los datos ##
##########################################

#Instalar y cargar los paquetes necesarios
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

if(!require(SummarizedExperiment)) {
  BiocManager::install("SummarizedExperiment")
}
require(SummarizedExperiment)

library(readr)

#Cargamos los archivos
DataInfo_S013 <- read_csv("Datasets/2018-MetabotypingPaper/DataInfo_S013.csv")
DataValues_S013 <- read_csv("Datasets/2018-MetabotypingPaper/DataValues_S013.csv")
head(DataInfo_S013)
head(DataValues_S013)

# Para trabajar con las trablas las transformamos a Dataframes
DataInfo_S013 <- as.data.frame(DataInfo_S013)
DataValues_S013 <- as.data.frame(DataValues_S013)

# Asignamos la primera columna como rownames
rownames(DataInfo_S013) <- DataInfo_S013[,1]
DataInfo_S013[,1] <- NULL

rownames(DataValues_S013) <- DataValues_S013[,1]
DataValues_S013[,1] <- NULL

# Extraemos la información de la matriz de valores
values <- DataValues_S013[,6:695]
values <- t(values)

# Extraemos la información de las columnas o samples
samples <- DataValues_S013[,1:5]

# Extraemos la información de las filas o features
features <- DataInfo_S013[6:695,]

# Cargamos la información del dataset como metadata con la función readLines()
metadata <- readLines("Datasets/2018-MetabotypingPaper/description.md")




#####################################################
## Creación del Contenedor "Summarized Experiment" ##
#####################################################

#contruimos el contenedor Summarized Experiment:
se <- SummarizedExperiment(
  assays = list(counts = values),
  rowData = features,
  colData = samples,
  metadata = list(title = "Metabotypes of response to bariatric surgery independent 
                    of the magnitude of weight loss",
                  repository = "https://github.com/nutrimetabolomics/Metabotyping2018",
                  info = metadata)
)

#Comprobarmos que el contenedor "Summarized Experiment" este correctamente generado:
se
str(se)

# Visualizamos parte de los datos y metadatos
head(assays(se)$counts) #Matriz de valores
rowData(se) #filas (características)
colData(se) #columnas (muestras)
metadata(se) #metadatos




##############################
## Exploración de los datos ##
##############################

# Resumenes estadísticos (de las variables peso y bmi)
row_weight <- which(grepl(paste("PESO", collapse = "|"), rownames(se)))
summary(t(assays(se)$counts [row_weight,]))

row_bmi <- which(grepl(paste("bmi", collapse = "|"), rownames(se)))
summary(t(assays(se)$counts [row_bmi,]))

# Representación en histogramas (del las variables peso y bmi)
opt <- par(mfrow=c(2,4))
hist(assays(se)$counts [10,], xlab = "Peso", main = "PESO_T0", xlim = c(50,210), ylim = c(0,20))
hist(assays(se)$counts [182,], xlab = "Peso", main = "PESO_T2", xlim = c(50,210), ylim = c(0,20))
hist(assays(se)$counts [356,], xlab = "Peso", main = "PESO_T4", xlim = c(50,210), ylim = c(0,20))
hist(assays(se)$counts [528,], xlab = "Peso", main = "PESO_T5", xlim = c(50,210), ylim = c(0,20))
hist(assays(se)$counts [11,], xlab = "BMI", main = "bmi_T0", xlim = c(20,70), ylim = c(0,20))
hist(assays(se)$counts [183,], xlab = "BMI", main = "bmi_T2", xlim = c(20,70), ylim = c(0,20))
hist(assays(se)$counts [357,], xlab = "BMI", main = "bmi_T4", xlim = c(20,70), ylim = c(0,20))
hist(assays(se)$counts [529,], xlab = "BMI", main = "bmi_T5", xlim = c(20,70), ylim = c(0,20))
par(opt)

# Determinamos las filas que se corresponden a las variables del colesterol
Metabol <- c("COL", "LDL", "HDL", "VLDL")
filas <- which(grepl(paste(Metabol, collapse = "|"), rownames(se)))
filas

# Representación en diagramas de cajas y bigotes (para las variables del colesterol)
groupColors <- c("green","green2","green3","green4")
opt <- par(mfrow=c(1,2))
boxplot(t(assays(se)$counts [c(18,190,364,536),]), col=groupColors, main="Total Cholesterol",
        xlab="Times",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)
boxplot(t(assays(se)$counts [c(19,191,365,537),]), col=groupColors, main="LDL",
        xlab="Times",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)
boxplot(t(assays(se)$counts [c(20,192,366,538),]), col=groupColors, main="HDL",
        xlab="Times",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)
boxplot(t(assays(se)$counts [c(21,193,367,539),]), col=groupColors, main="VLDL",
        xlab="Times",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)
par(opt)

# Representación en diagramas de cajas y bigotes de los valores de las muestras
boxplot(assays(se)$counts, col=groupColors,
        xlab="Samples",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)

# Representación en diagramas de cajas y bigotes del logaritmo de los valores de las muestras
boxplot(log(assays(se)$counts), col=groupColors, main=,
        xlab="Samples",
        ylab="Level", las=2, cex.axis=0.7, cex.main=0.7)

#Comprobación y eliminación de variables con NAs
all(is.finite( assays(se)$counts )) #Comprobamos si hay valores no finitos
counts_filtered <- na.omit(assays(se)$counts) #Eliminamos variables con NAs
all(is.finite(counts_filtered)) #Comprobamos si hay valores no finitos

#Analisis de componentes principales (PCA)
PC_se <- prcomp(t(counts_filtered), scale= FALSE)

#Representación del PCA
loads <- round(PC_se$sdev^2/sum(PC_se$sdev^2)*100,1)
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(PC_se$x[,1:2],xlab=xlab,ylab=ylab,
     main ="Principal components (PCA)")
names2plot<-names(as.data.frame(counts_filtered))
text(PC_se$x[,1],PC_se$x[,2],names2plot, pos=3, cex=.6)

#Coeficientes de las componentes principales 1 y 2
PC1 <- PC_se$rotation[,1]
pC1_sorted <- sort(PC1, decreasing = TRUE)
pC1_sorted[1:20]

PC2 <- PC_se$rotation[,2]
pC2_sorted <- sort(PC2, decreasing = TRUE)
pC2_sorted[1:20]

#Representación del PCA en base a la cirugía
palette <- c("red", "blue")
surgery <- as.factor(colData(se)$SURGERY)
colors_surgery <- palette[as.numeric(surgery)]

plot(PC_se$x[,1:2],xlab=xlab,ylab=ylab, col = colors_surgery,
     main ="Principal components (PCA) - Surgery")
text(PC_se$x[,1],PC_se$x[,2],names2plot, pos=3, cex=.6)
legend( "bottomright"  , inset = c(0.01,0.01), cex =1, bty = "y", legend = c("by pass", "tubular"), 
        text.col = c("red", "blue"), col = c("red", "blue"), pt.bg = c("red","blue"), pch = c(1,2))

#Representación del PCA en base al genero
gender <- as.factor(colData(se)$GENDER)
colors_gender <- palette[as.numeric(gender)]

plot(PC_se$x[,1:2],xlab=xlab,ylab=ylab, col = colors_gender,
     main ="Principal components (PCA) - Gender")
text(PC_se$x[,1],PC_se$x[,2],names2plot, pos=3, cex=.6)
legend( "bottomright"  , inset = c(0.01,0.01), cex =1, bty = "y", legend = c("Female", "Male"), 
        text.col = c("red", "blue"), col = c("red", "blue"), pt.bg = c("red","blue"), pch = c(1,2))

#Representación del PCA en base al grupo
group <- as.factor(colData(se)$Group)
colors_group <- palette[as.numeric(group)]

plot(PC_se$x[,1:2],xlab=xlab,ylab=ylab, col = colors_group,
     main ="Principal components (PCA) - Group")
text(PC_se$x[,1],PC_se$x[,2],names2plot, pos=3, cex=.6)
legend( "bottomright"  , inset = c(0.01,0.01), cex =1, bty = "y", legend = c("1", "2"), 
        text.col = c("red", "blue"), col = c("red", "blue"), pt.bg = c("red","blue"), pch = c(1,2))

#Agrupación jerárquica con el método average y representación
clust.euclid.average <- hclust(dist(t(counts_filtered)),method="average")
plot(clust.euclid.average, hang=-1)

#Agrupación jerárquica con el método complete y representación
clust.euclid.complete <- hclust(dist(t(counts_filtered)),method="complete")
plot(clust.euclid.complete, hang=-1)


#######################################
## Reposición de los datos en Github ##
#######################################

# Guardamos el objeto contenedor se en un archivo .Rda
save(se, file = "Contenedor_SummarizedExperiment.Rda")

#guardamos la matriz de valores del objeto contenedor en formato .csv
write.csv(assays(se)$counts, file = "Datos.csv", row.names = TRUE)

# Convertir los metadatos del dataset, de las filas y colunmas a formato markdown
metadata_dataset <- knitr::kable(metadata(se), format = "markdown")
metadatos_samples_md <- knitr::kable(as.data.frame(colData(se)), format = "markdown")
metadatos_features_md <- knitr::kable(as.data.frame(rowData(se)), format = "markdown")
# Crear el archivo markdown
writeLines(c("# Metadatos del Dataset", metadata_dataset, "# Metadatos de columnas", 
             metadatos_samples_md, "# Metadatos de filas", metadatos_features_md), "metadatos.md")

