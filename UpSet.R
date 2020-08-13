library(UpSetR)

# variants
sourceTable <- read.csv('results/UpSet_panelVariants.csv')
colnames(sourceTable) = c('X',
                          'xGen Exome Research Panel v2',
                          'xGen Pan-Cancer Panel v2.4',
                          'TrueSight Oncology 500',
                          'Ion AmpliSeq Comprehensive Cancer Panel',
                          'QIAseq Targeted Human Comprehensive Cancer Panel',
                          'AVENIO ctDNA Targeted Kit',
                          'AVENIO ctDNA Expanded Kit',
                          'AVENIO ctDNA Surveillance Kit',
                          'FoundationOne CDx',
                          'FoundationOne Liquid ')

upset(sourceTable,
      nsets = 11,
      nintersects = NA,
      point.size = 5,
      line.size = 1.5,
      mainbar.y.label = "Variants in the meta-knowledgebase targeted by each panel",
      sets.x.label = "Number of variants in the meta-knowledgebase targeted by each panel",
      text.scale = c(1.5,1.5,1.1,1.5,2,2)
)
