#Pre-requisites:
```r
#Dependencies pacakage
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
############
BiocManager::install(c("Biostrings","liftOver", "rtracklayer", "BSgenome","GenomeInfoDb","GenomicRanges","IRanges"))

##########
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmusculus.UCSC.mm10"))

#############
install.packages(c("stringr","devtools"))

# development version from GitHub:
require("devtools")
devtools::install_github("XuLabSJTU/COL")
