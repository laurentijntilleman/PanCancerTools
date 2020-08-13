# PanCancerTools

This GitHub directory contains the script used in manuscript 'Pan-cancer pharmacogenetics: targeted sequencing panels or exome sequencing?'
Please cite:
- Laurentijn Tilleman , Björn Heindryckx, Dieter Deforce & Filip Van Nieuwerburgh *Pan-cancer pharmacogenetics: targeted sequencing panels or exome sequencing?* Pharmacogenomics 2020, (doi: [10.2217/pgs-2020-0035](http://dx.doi.org/10.2217/pgs-2020-0035))

## Data preprocessing

* Harvest the databases:
  * Run the harvester tool in https://github.com/laurentijntilleman/g2p-aggregator for CGI, CIViC, OncoKB, DEPO, and JAX-CKB
  * Copy cgi.json, civic.json, oncokb.json, depo.json, and jax.json from g2p-aggregator/harvester/ to the data folder
  * Copy cosmic_lookup_table.tsv from g2p-aggregator/harvester/ to the data folder

* Pan-Cancer data:
  * Download Mutation.CTAT.3D.Scores.txt from https://gdc.cancer.gov/about-data/publications/pancan-driver in the data folder (Bailey *et al.*, 2018)
  * Download 1-s2.0-S2211124718303954-mmc2.xlsx from supplemental information from Gao *et al.* (2018)

* Panels:
  * Data received from the supplier can be found in the data folder
  * Sources of the data can be found in the manuscript (Tilleman *et al.*, 2020)

* GRCh38 to GRCh37 conversion:
  * Download the chain file from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz into the data folder

* UCSC exone coordinates:
  * Download UCSC_exons.bed from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables):
    * track: UCSC Genes
    * table: knownGene
    * output format: BED - browser extensible data
    * press get output
    * select option CDS, intron, 3' exons
    * press get BED
  * Modify the file to separate the transcript name of the rest of information:

  `awk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"$6}}' UCSC_exon.bed > UCSC_exons_modif.bed`

  `awk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"$6}}' UCSC_intron.bed > UCSC_intron_modif.bed`

  `awk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"$6}}' UCSC_3_exon.bed > UCSC_3_exons_modif.bed`

  * Download USCS_canonical.bed from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
    * track: UCSC Genes
    * table: knownCanonical
    * output format: select fields from primary and related tables
    * press get output
    * select fields from hg19.knownCanonical: chrom, chromStart, chromEnd, transcript
    * select fields from hg19.kgXref: geneSymbol
    * press get output
  * Join the sorted files based on the transcript identificator:

  `join -1 4 -2 4 <(sort -k4 UCSC_exons_modif.bed ) <(sort -k4 UCSC_canonical.bed) | awk '{print $2"\t"$3"\t"$4"\t"$10"\t"$5"\t"$6}' | bedtools sort -i "-" > UCSC_exons_modif_canonical.bed`

  `join -1 4 -2 4 <(sort -k4 UCSC_intron_modif.bed ) <(sort -k4 UCSC_canonical.bed) | awk '{print $2"\t"$3"\t"$4"\t"$10"\t"$5"\t"$6}' | bedtools sort -i "-" > UCSC_introns_modif_canonical.bed`

  `join -1 4 -2 4 <(sort -k4 UCSC_3_exons_modif.bed ) <(sort -k4 UCSC_canonical.bed) | awk '{print $2"\t"$3"\t"$4"\t"$10"\t"$5"\t"$6}' | bedtools sort -i "-" > UCSC_3_exons_modif_canonical.bed`

## Files in the data folder

* CCP.20131001.submitted.bed : The targeted regions of the Ion AmpliSeq Comprehensive Cancer Panel
* List of drugs that are corrected comma separtated
* cosmic_lookup_table.tsv : COSMIC lookup table from the harvester tool described in Wagner et al. (2018)
* expanded.txt : The genes of the AVENIO ctDNA Expanded Kit
* Expanded_Panel_Regions_Information_Tumor_08374325001.xlsx : The targeted regions of the AVENIO ctDNA Expanded Kit
* FoundationOne.txt : The genes of the FoundationOne Liquid panel from FoundationMedicine
* FoundationOneCDx.txt : The genes of the FoundationOne CDx panel from FoundationMedicine
* FoundationOneCDxIntron.txt : The intronic regions of the genes of the FoudnationOne CDx panel from FoundationMedicine
* FoundationOneIntron.txt : The intronic regions of the genes of the FoudnationOne Liquid panel from FoundationMedicine
* hg38ToHg19.over.chain : The chain file for converting the GRCh38 coordinates to the GRCh37 coordinates
* illuminaTSO500_list.txt : The genes of the TrueSight Oncology 500 panel
* Mutation.CTAT.3D.Scores.txt : The Pan-Cancer mutations (Bailey et al. 2018)
* PanCancer V2 gene list.csv : The genes of the xGen Pan-Cancer v2.4 Panel
* qiagen_compCancer_list.txt : The genes of the QIAseq Targeted Human Comprehensive Cancer Panel
* surveilance.txt : The genes of the AVENIO ctDNA Surveillance Kit
* Surveillance_Panel_Regions_Information_Tumor_08374333001.xlsx : The targeted regions of the AVENIO ctDNA Surveillance Kit
* target.txt : The genes of the AVENIO ctDNA Targeted Kit
* Targeted_Panel_Regions_Information_Tumor_08374317001.xlsx : The targeted regions of the AVENIO ctDNA Targeted Kit

## Scripts

* Python2:
 * PanCancer_mutations.py: lookup of the genomic coordinates of the Pan-Cancer driver mutations

* Python3:
 * loadMongoDB.py: load output of the harvester tool in the mongoDB
 * queryPGX.py: compared the targeted sequencing panels based on genes
 * distinctCA.py: removes duplicated pharmacogenomic interactions and build the main tables

* R-script:
 * UpSet.R: make Figure 2

# Bibliography
M. H. Bailey *et al.*, *Comprehensive Characterization of Cancer Driver Genes and Mutations*, Cell, vol. 173, no. 2, pp. 371-385.e18, Apr. 2018.

A. H. Wagner *et al.*, *A harmonized meta-knowledgebase of clinical interpretations of cancer genomic variants*, bioRxiv, p. 366856, Oct. 2018.

Q. Gao *et al.*, *Diver fusions and their implications in the development and treatment of human cancers*, Cell, vol. 23, no. 1, pp. 227–238, 2018.

L. Tilleman *et al.*, *Pan-cancer pharmacogenetics: targeted sequencing panels or exome sequencing?* Pharmacogenomics 2020, (doi: [10.2217/pgs-2020-0035](http://dx.doi.org/10.2217/pgs-2020-0035))
