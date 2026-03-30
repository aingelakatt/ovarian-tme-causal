#-------Dataset Acquisition and Exploration-------

library(Seurat) 
library(tidyverse)
library(pcalg)
library(bnlearn)



#Load Dataset 
tme_data <- readRDS("~/Downloads/sc data in ovarian cancer zemin zhang/raw_object.rds")

#Metadata 
print(head(tme_data@meta.data))
