
#' Dependencies pacakage
#' @import BiocManager
#' @importFrom BiocManager install
#' @import rtracklayer
#' @import Biostrings
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @import liftOver
#' @import stringr
#' @importFrom utils write.table
#' @importFrom IRanges IRanges
#' @importFrom stats lm
#'@importFrom stats cooks.distance
#'@importFrom stats residuals
#'@import RobustRankAggreg
#'@import parallel
#'@import doParallel
#'@import foreach
#'@importFrom data.table data.table

##############################################################################################

##########################################################
##################################### Detect outliers based on COOKs distance
#'@title outlier
#'@name outlier
#'@description
#'outlier function used to detect outliers based on COOKs distance
#'@param cirRNA cirularRNA data
#'@param STAR_output output file from STAR
#'@param k the Cook's distance cutoff for outliers
#example outlier(cirRNA=cirRNA,STAR_output=STAR_output)
############
#'@export
outlier<-function(cirRNA,STAR_output,k=4){
  BCPAP.rep.1<-cirRNA
  BCPAP.rep.1.cleanedSJ.out <-STAR_output
  BCPAP.rep.1 = BCPAP.rep.1[which(BCPAP.rep.1$X.junction_reads>1),]
  rownames(BCPAP.rep.1) <- 1:nrow(BCPAP.rep.1)

  jilunum = 0
  for (i in BCPAP.rep.1$chr) {
    jilunum = jilunum + 1
    BCPAP.rep.1[jilunum,"hu_bsjl"] = paste(BCPAP.rep.1[jilunum,"chr"],
                                           BCPAP.rep.1[jilunum,"circRNA_start"],
                                           BCPAP.rep.1[jilunum,"strand"],sep="@")
    BCPAP.rep.1[jilunum,"hu_bsjr"] = paste(BCPAP.rep.1[jilunum,"chr"],
                                           BCPAP.rep.1[jilunum,"circRNA_end"],
                                           BCPAP.rep.1[jilunum,"strand"],sep="@")
    BCPAP.rep.1[jilunum,"Info"] = paste(BCPAP.rep.1[jilunum,"chr"],
                                        BCPAP.rep.1[jilunum,"circRNA_start"],
                                        BCPAP.rep.1[jilunum,"circRNA_end"],
                                        BCPAP.rep.1[jilunum,"strand"],sep="@")
  }

  BCPAP.rep.1["BsjRds"] = BCPAP.rep.1["X.junction_reads"]

  # left site total circular reads
  BCPAP.rep.1["hu_bsjlRds"] = 0
  jilunum = 0
  for (i in BCPAP.rep.1$hu_bsjl){
    jilunum = jilunum + 1
    indexs1 = which(BCPAP.rep.1$hu_bsjl==i)
    for (index1 in indexs1) {
      BCPAP.rep.1[jilunum,"hu_bsjlRds"] = BCPAP.rep.1[jilunum,"hu_bsjlRds"]+BCPAP.rep.1[index1,"X.junction_reads"]
    }
  }

  # right site total circular reads
  BCPAP.rep.1["hu_bsjrRds"] = 0
  jilunum = 0
  for (i in BCPAP.rep.1$hu_bsjr){
    jilunum = jilunum + 1
    indexs1 = which(BCPAP.rep.1$hu_bsjr==i)
    for (index1 in indexs1) {
      BCPAP.rep.1[jilunum,"hu_bsjrRds"] = BCPAP.rep.1[jilunum,"hu_bsjrRds"]+BCPAP.rep.1[index1,"X.junction_reads"]
    }
  }


  # clean STAR output information
  BCPAP.rep.1.cleanedSJ.out = BCPAP.rep.1.cleanedSJ.out[which(BCPAP.rep.1.cleanedSJ.out$V7>1),]
  rownames(BCPAP.rep.1.cleanedSJ.out) <- 1:nrow(BCPAP.rep.1.cleanedSJ.out)

  jilunum = 0
  for (i in BCPAP.rep.1.cleanedSJ.out$V1) {
    jilunum = jilunum + 1
    if(BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"]==1){
      BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"] = "+"
    }
    if(BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"]==2){
      BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"] = "-"
    }
    if(BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"]==0){
      BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"] = "none"
    }

    BCPAP.rep.1.cleanedSJ.out[jilunum,"hu_bsjl"] = paste(BCPAP.rep.1.cleanedSJ.out[jilunum,"V1"],
                                                         as.character(as.numeric(BCPAP.rep.1.cleanedSJ.out[jilunum,"V2"])-1),
                                                         BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"],sep="@")

    BCPAP.rep.1.cleanedSJ.out[jilunum,"hu_bsjr"] = paste(BCPAP.rep.1.cleanedSJ.out[jilunum,"V1"],
                                                         as.character(as.numeric(BCPAP.rep.1.cleanedSJ.out[jilunum,"V3"])+1),
                                                         BCPAP.rep.1.cleanedSJ.out[jilunum,"V4"],sep="@")
  }


  # left site total linear reads
  BCPAP.rep.1["hu_bsjlLsjRds"] = 0
  jilunum = 0
  for (i in BCPAP.rep.1$hu_bsjl){
    jilunum = jilunum + 1
    indexs1 = which(BCPAP.rep.1.cleanedSJ.out$hu_bsjr==i)
    for (index1 in indexs1) {
      BCPAP.rep.1[jilunum,"hu_bsjlLsjRds"] = BCPAP.rep.1[jilunum,"hu_bsjlLsjRds"]+BCPAP.rep.1.cleanedSJ.out[index1,"V7"]
    }
  }

  # right site total linear reads
  BCPAP.rep.1["hu_bsjrLsjRds"] = 0
  jilunum = 0
  for (i in BCPAP.rep.1$hu_bsjr){
    jilunum = jilunum + 1
    indexs1 = which(BCPAP.rep.1.cleanedSJ.out$hu_bsjl==i)
    for (index1 in indexs1) {
      BCPAP.rep.1[jilunum,"hu_bsjrLsjRds"] = BCPAP.rep.1[jilunum,"hu_bsjrLsjRds"]+BCPAP.rep.1.cleanedSJ.out[index1,"V7"]
    }
  }

  BCPAP.rep.1["TotalRds"] = BCPAP.rep.1["hu_bsjlRds"]+BCPAP.rep.1["hu_bsjrRds"]+BCPAP.rep.1["hu_bsjlLsjRds"]+BCPAP.rep.1["hu_bsjrLsjRds"]
  BCPAP.rep.1["BsjRt"] = 2*BCPAP.rep.1["X.junction_reads"]/BCPAP.rep.1["TotalRds"]


  ### calculate cook's distance and chose outliers
  huRdsWO0 = log10(BCPAP.rep.1[,"TotalRds"])
  huRtWO0 = log10(BCPAP.rep.1[,"BsjRt"])
  huMod = lm(huRtWO0~huRdsWO0)
  # calculate cook's distance
  huCooksd <- cooks.distance(huMod)
  BCPAP.rep.1[,"huCooksd"] = 0

  jilunum=0
  for (i in BCPAP.rep.1$Info) {
    jilunum = jilunum+1
    BCPAP.rep.1[jilunum,"huCooksd"] = huCooksd[jilunum]
  }

  ###确认site的outliers状态，overLm是1，underLm是-1，nonOutliers是0
  BCPAP.rep.1["OutlStatus"] = 0
  huCsIndex1 = huCooksd > k*mean(huCooksd, na.rm=T)
  huCsIndex2 = (huCooksd > k*mean(huCooksd, na.rm=T))& (residuals(huMod)>0)
  huCsIndex2 = which(huCsIndex2==TRUE)
  huCsIndex3 = (huCooksd > k*mean(huCooksd, na.rm=T))& (residuals(huMod)<0)
  huCsIndex3 = which(huCsIndex3==TRUE)


  for (i in huCsIndex2) {
    BCPAP.rep.1[i,"OutlStatus"] = 1
  }
  for (i in huCsIndex3) {
    BCPAP.rep.1[i,"OutlStatus"] = -1
  }
  testDt = BCPAP.rep.1[15:24]
  jilunum = 0
  for (i in testDt$Info) {
    jilunum = jilunum + 1
    infos = strsplit(i,"@")[[1]]
    testDt[jilunum,"chr"] = paste("chr",infos[1],sep = "")
    testDt[jilunum,"start"] = infos[2]
    testDt[jilunum,"end"] = infos[3]
    testDt[jilunum,"strand"] = infos[4]
  }
  return(testDt)
}

##############################################################################
############################### function of detecting motif
##############################################################################

#'@title cmotif
#'@name cmotif
#'@description
#'cmotif function used to find downstream and upstream motif (2nt) from start and end of the chromosome,respectively.
#'@param genome.UCSC genome data, by default human genome file
#'@param genomeOther other DNA sequence (FASTA) instead of Homo sapiens
#'@param cRNA circular RNA data
# example cmotif(genome.UCSC="Hgenome",cRNA=cRNA)
#'@export
cmotif <- function(genome.UCSC = "Hgenome", genomeOther, cRNA) {

  # Read FASTA file
  gen <- genome(genome.UCSC)
  chr <- names(GenomeInfoDb::seqlengths(gen))

  # Filter cRNA to only include chromosomes in gen
  cRNA <- cRNA[cRNA$chr %in% chr, ]

  # Set up parallel cluster
  cl <- makeCluster(detectCores() - 1)  # Use all but one core
  registerDoParallel(cl)
  #
  chrF<-NULL
  # Parallel processing for each chromosome
  output <- foreach(chrF = unique(cRNA$chr), .combine = rbind, .packages = c("Biostrings","GenomeInfoDb")) %dopar% {
    genI <- as.character(gen[[chrF]])  # Extract genome sequence for the current chromosome

    ######### Forward Strand
    cRNAI <- cRNA[cRNA$chr == chrF & cRNA$strand == "+", ]
    outtF<-data.table()
    if (nrow(cRNAI) > 0) {
      motifs_start <- mapply(substr,genI, as.numeric(cRNAI$start) - 2, as.numeric(cRNAI$start) - 1)
      motifs_end <-  mapply(substr,genI, as.numeric(cRNAI$end) + 1, as.numeric(cRNAI$end) + 2)
      cRNAI$motif <- paste(motifs_start, motifs_end, sep = "")
      outtF<-cRNAI
    }

    ######### Reverse Strand
    cRNAR <- cRNA[cRNA$chr == chrF & cRNA$strand == "-", ]
    outtR <- data.table()
    if (nrow(cRNAR) > 0) {
      motifs_start <- mapply(substr,genI, as.numeric(cRNAR$start) - 2, as.numeric(cRNAR$start) - 1)
      motifs_end <-  mapply(substr,genI, as.numeric(cRNAR$end) + 1, as.numeric(cRNAR$end) + 2)
      motifs <- paste(motifs_start, motifs_end, sep = "")
      cRNAR$motif <- as.character(reverseComplement(DNAStringSet(motifs)))
      outtR<-cRNAR
    }

    # Combine forward and reverse strand results
    outchr <- rbind(outtF, outtR)
    return(outchr)
  }

  # Stop parallel cluster
  stopCluster(cl)

  # Return the output
  return(as.data.frame(output))
}
###########################################################################################
##########################################################################################
########################### Run liftover
#'@title liftOVER
#'@name liftOVER
#'@description
#'liftOVER function is built based on the UCSC liftover tools, which are designed for lifting features from specific region rather than
#'an interval region from one genome to other.
#'@param output circular RNA data, usually data have chromosome,start,end and so on,
#'liftOVER function accept the data format like chr,start,end,strand and so on
#'@param chainObject chainObject usually imported using \code{\link[rtracklayer]{import.chain}} and downloaded data from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
#'@return lifting result with chromosome, start, end, strand and so on
#'@seealso \code{\link[rtracklayer]{liftOver}} for more details
# example
# lifData<-liftOVER(outlier,hg38ToMm10.over.chain) #human to mouse:hg38ToMm10.over.chain ; mouse to human:mm10ToHg38.over.chain
#'@references https://genome.ucsc.edu/cgi-bin/hgLiftOver
#'@export
liftOVER<-function(output,chainObject){
  chainObject<-rtracklayer::import.chain(chainObject)
  #########################
  liftover.data.F<-data.frame()
  for (ii in 1:length(unique(output$chr))){
    cRNA.F<-output[output$chr %in% unique(output$chr)[ii],]

    ################
    liftover.data.FF<-data.frame()
    for (lf in 1:dim(cRNA.F)[1]){
      # specify coordinates to liftover
      grObject.S <- GenomicRanges::GRanges(seqnames=unique(cRNA.F$chr), ranges=IRanges(start=as.numeric(cRNA.F$start[lf])+1,end=as.numeric(cRNA.F$start[lf]))+1,strand =as.character(cRNA.F$strand[lf]))
      result.S<- as.data.frame(rtracklayer::liftOver(grObject.S, chainObject))
      if(dim(result.S)[1]>0){
        end<-result.S$start
      }else{
        end<-0
      }
      ####################
      grObject.E <- GenomicRanges::GRanges(seqnames=unique(cRNA.F$chr), ranges=IRanges(start=as.numeric(cRNA.F$end[lf])+1,end=as.numeric(cRNA.F$end[lf]))+1,strand =as.character(cRNA.F$strand[lf]))
      result.E<- as.data.frame(rtracklayer::liftOver(grObject.E, chainObject))
      if(dim(result.E)[1]>0){
        start<-result.E$start
      }else{
        start<-0
      }
      strand=ifelse(dim(result.E)[1]>0,as.character(result.E$strand),as.character(result.S$strand))
      if(dim(result.S)[1] & dim(result.E)[1]>0){
        if(as.character(result.E$seqnames)[1]==as.character(result.S$seqnames)[1]){
          seqnames<-ifelse(dim(result.E)[1]>0,as.character(result.E$seqnames),as.character(result.S$seqnames))
        }else{
          seqnames<-"Multi"
        }
      }else{
        seqnames<-seqnames<-ifelse(dim(result.E)[1]>0,as.character(result.E$seqnames),as.character(result.S$seqnames))
      }
      lift<-cbind(seqnames=seqnames,strand=strand,start,end,seqnamesB=unique(cRNA.F$chr),strandB=cRNA.F$strand[lf],startB=cRNA.F$start[lf],endB=cRNA.F$end[lf],motifB=cRNA.F$motif[lf])
      liftover.data.FF<-rbind(liftover.data.FF,data.frame(lift))
    }

    liftover.data.F<-rbind(liftover.data.F,liftover.data.FF)
  }
  ################# extension
  liftoverData <- liftover.data.F[liftover.data.F$start != 0 & liftover.data.F$end != 0 & liftover.data.F$seqnames != "Multi",]
  liftfinal <- data.frame()

  for (e in 1:dim(liftoverData)[1]) {
    liftoverData11 <- liftoverData[e,]
    liftoverData1<-liftoverData11
    # Check for missing values in start and end columns
    if (!is.na(liftoverData1$start) && !is.na(liftoverData1$end)) {
      if (as.numeric(liftoverData1$start) < as.numeric(liftoverData1$end)) {
        liftoverData1$end <- as.numeric(liftoverData1$end) + 1
      } else {
        liftoverData1$start <- as.numeric(liftoverData1$end)- 1
        liftoverData1$end <- as.numeric(liftoverData11$start)
      }
      liftfinal <- rbind(liftfinal, liftoverData1)
    } else {
      # Handle missing values (if needed)
      # For example, you could skip the current iteration or handle it differently based on your requirements
    }
  }
  ###save the result
  return(liftfinal)
}

####################################################################################
############################## comparison of human to mouse motif or others
#'@title conserve
#'@name conserve
#'@description
#'conserve function is built to identify putatively functional back-splicing based on the molecular error hypothesis.COL combined (1) substantially higher back-splicing rates than expected from the total splicing amounts; (2) conserved splicing motifs between different species; (3) high back-splicing levels to identify functional back-splicing.
#'@param cirRNA cirularRNA data
#'@param STAR_output output file from STAR
#'@param genome1 DNA sequence (FASTA), by default Homo sapiens
#'@param genomeOther1 other DNA sequence (FASTA) instead of Homo sapiens
#'@param genome2 DNA sequence (FASTA), by default Mus musculus
#'@param genomeOther2 other DNA sequence (FASTA) instead of Mus musculus
#'@param chainObject chainObject usually imported using \code{\link[rtracklayer]{import.chain}} and downloaded data from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
#'@param outpath output path
#'@return Conserved motif result
#'@seealso position specific liftover function: \code{\link[COL]{liftOVER}} , outlier identification function: \code{\link[COL]{outlier}} , Moftif identification function: \code{\link[COL]{cmotif}}
# example
# conserve(cirRNA=cirRNA,STAR_output=cirRNA,genome1="Hgenome",genome2="Mgenome",chainObject=hg38ToMm10.over.chain,"/outpath")
#'@export
conserve<-function(cirRNA,STAR_output,genome1="Hgenome",genomeOther1,genome2="Mgenome",genomeOther2,chainObject,outpath){
  message("COL::INFO:Computation takes times... ...")
  #Run outlier
  cRNA<-outlier(cirRNA,STAR_output)
  #Run comotif
  motifDataP<-cmotif(genome.UCSC=genome1,genomeOther=genomeOther1,cRNA=cRNA)
  #Run liftOVER
  liftData<-liftOVER(output =motifDataP,chainObject)
  #compare
  colnames(liftData)<-c("chr","strand","start","end","chrB","strandB","startB","endB","motifB")
  liftData$start<-as.numeric(liftData$start)+1
  convert<-cmotif(genome.UCSC=genome2,genomeOther=genomeOther2,cRNA=liftData)
  convert$conserve<-ifelse(convert$motifB==convert$motif,1,0)
  #combined outputs of outlier and comotif
  jilunum = 0
  for (i in convert$chr) {
    jilunum = jilunum + 1
    convert[jilunum,"Info"] = paste(substring(convert[jilunum,"chrB"],4),
                                    convert[jilunum,"startB"],
                                    convert[jilunum,"endB"],
                                    convert[jilunum,"strandB"],sep="@")
  }
  jilunum = 0
  for (i in cRNA$Info) {
    jilunum = jilunum + 1
    index1 = which(convert$Info==i)
    if(length(index1)>0){
      cRNA[jilunum,"huMoConMotif"] = convert[index1[1],"conserve"]
    }
  }

  #RRA rank in COL
  huBsjRdsInd=order(cRNA[,"BsjRds"],decreasing = TRUE)
  cRNA[,"BSJRank"] = NA
  jilunum = 0
  for (i in huBsjRdsInd) {
    jilunum = jilunum + 1
    cRNA[i,"BSJRank"] = jilunum
  }

  num1 = length(which((cRNA$huMoConMotif==1)&(cRNA$OutlStatus==1)))
  cRNA["COL_Candidates"] = 0
  cRNA[which((cRNA$huMoConMotif==1)&(cRNA$OutlStatus==1)&(cRNA$BSJRank<=num1)),"COL_Candidates"] = 1

  testDt = cRNA[which(cRNA$COL_Candidates==1),]
  huBsjRtInd1 = order(testDt[,"BsjRds"],decreasing = TRUE)
  huBsjRtInd2 = order(testDt[,"huCooksd"],decreasing = TRUE)

  # RRA
  glist <- list(gene_set1=testDt[huBsjRtInd1,"Info"],
                gene_set2=testDt[huBsjRtInd2,"Info"])

  #应用RRA算法，对基因进行整合排序
  ag=RobustRankAggreg::aggregateRanks(glist,N=length(cRNA$Info))
  #添加基因出现的次数
  cRNA[,"RRAScore"]=NA
  jilunum = 0
  for (i in ag$Name) {
    jilunum=jilunum+1
    index1=which(cRNA[,"Info"]==i)
    cRNA[index1,"RRAScore"]= ag$Score[jilunum]
  }

  huBsjCombinedInd=order(cRNA[,"RRAScore"],decreasing = FALSE)
  cRNA[,"RRARank"] = NA
  jilunum = 0
  for (i in huBsjCombinedInd) {
    jilunum = jilunum + 1
    cRNA[i,"RRARank"] = jilunum
  }

  FunctionalCirc<-cRNA[which(cRNA$RRAScore<0.05),]
  ###save the result
  outpathO<-stringr::str_c(outpath,"/OutliersData.txt")
  write.table(cRNA,outpathO,row.names =F)
  outpathF<-stringr::str_c(outpath,"/Motifdata.txt")
  write.table(convert,outpathF,row.names =F)
  outpathAll<-stringr::str_c(outpath,"/FunctionalCirc.txt")
  write.table(FunctionalCirc,outpathAll,row.names =F)
  list(OutliersData=cRNA,FunctionalCirc=FunctionalCirc,finalData=convert)
}
###############
################### Genome sequences
#'@title genome
#'@name genome
#'@description
#'Genome sequences for Homo sapiens (see more \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}) and Mus musculus (see more \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}})
#'
#'@param genome.UCSC genome.UCSC DNA sequence (FASTA format); for example Hgenome for Homo sapiens and Mgenome for Mus musculus
#'@param genomeOther other DNA sequence (FASTA) instead of Homo sapiens; if user use independent genome FASTA where parameter #genome.UCSC="others".

#'@export
genome <- function(genome.UCSC, genomeOther) {
  gen <- NULL
  if (genome.UCSC == "Hgenome") {
    gen <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (genome.UCSC == "Mgenome") {
    gen <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else {
    gen <- genomeOther
  }
  return(gen)
}

#################################################################################################
#################################### END #######################################################
################################################################################################


