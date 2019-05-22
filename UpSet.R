setwd("/data/ownCloud/Werk/Doctoraat/Pharmacogenetics/CancerPanel/Paper/PanCancerTools")
library(UpSetR)

# database resources
sourceTable <- read.csv('restults/sourceTable.csv')
colnames(sourceTable) = c('X','CGI','CIViC','DEPO','JAX-CKB','OncoKB')

upset(sourceTable,
      point.size = 5,
      line.size = 1.5,
      mainbar.y.label = "Pharmacogenomic interactions originating from each knowledge base",
      sets.x.label = "Pharmacogenomic interactions per knowledge base",
      text.scale = c(1.2,1,1.2,1.5,2,2)
      )

variants
sourceTable <- read.csv('results/UpSet_panelVariants.csv')

upset(sourceTable,
      point.size = 5,
      line.size = 1.5,
      mainbar.y.label = "Variants targeted by each panel",
      sets.x.label = "Number of Variants targeted by each panel",
      text.scale = c(1.2,1,1.2,1.5,2,2)
)
