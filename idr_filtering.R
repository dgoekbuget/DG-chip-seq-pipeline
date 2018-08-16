library(ChIPpeakAnno)
library(GenomicAlignments)

idr_filtering <- function(path){
  #This function takes a file path for one experiment (with /bam, /bw and /macs2 subfolders containing the respective files)
  
  #Define IDRfilter function for 3 replicates
  
  IDRfilterN3 <- function(peaksA, peaksB, peaksC, bamfileA, bamfileB, bamfileC,
                          maxgap=-1L, minoverlap=0L, singleEnd=TRUE,
                          IDRcutoff=0.01){
    stopifnot(class(peaksA)=="GRanges")
    stopifnot(class(peaksB)=="GRanges")
    stopifnot(class(peaksC)=="GRanges")
    stopifnot(file.exists(bamfileA))
    stopifnot(file.exists(bamfileB))
    stopifnot(file.exists(bamfileC))
    ol <- findOverlapsOfPeaks(peaksA, peaksB, peaksC, 
                              maxgap=maxgap, 
                              minoverlap=minoverlap)
    ol <- ol$peaklist[grepl("\\/\\/\\/", names(ol$peaklist))][[1]]
    if(length(ol)<1) return(GRanges())
    names(ol) <- paste0("olp", 
                        formatC(1:length(ol), 
                                width=nchar(as.character(length(ol))),
                                flag="0"))
    coverage <- GenomicAlignments::summarizeOverlaps(features = ol, 
                                                     reads =  Rsamtools::BamFileList(c(bamfileA, bamfileB, bamfileC)),
                                                     mode=Union, 
                                                     ignore.strand = FALSE, 
                                                     singleEnd=singleEnd)
    idr <- idr::est.IDR(assay(coverage)/width(rowRanges(coverage)),
                        mu=2.07, sigma=1.34, rho=0.89, p=0.84)
    ol[idr$IDR<IDRcutoff]
  }
  
  #Define path to MACS2 output, Bam files and BigWig files
  ##Convert .xls MACS2 output to Grange
  peaks <- dir(file.path(path,"/macs2"),"peaks.xls")
  for (i in 1:length(peaks)){
    assign(paste0("peaks_N",i),toGRanges(file.path(path,"/macs2",peaks)[i],format = "MACS2"))
  }
  
  ##File path to experimental .bam files
  bams <- dir(file.path(path,"/bam"),".bam")
  for (i in 1:length(bams)){
    assign(paste0("bam_N",i),file.path(path,"/bam",bams)[i])
  }
  
  #IDR filtering
  
  if (length(peaks) == 2){
    idr_peaks <- IDRfilter(peaks_N1,peaks_N2,bam_N1,bam_N2)
  } else if (length(peaks) == 3) {
    idr_peaks <- IDRfilterN3(peaks_N1, peaks_N2, peaks_N3, bam_N1, bam_N2, bam_N3)
  } else if (length(peaks) == 1) {
    print("no IDR filtering done. n=1.")
    idr_peaks <- peaks_N1
  } else {
    print("Sample size unknown. IDR failed!")
  }
  
  #Save data
  saveRDS(idr_peaks, paste0(path,"/",gsub(".sorted.dedup.bam","",bams[1]),".",Sys.Date(),".idr.rds"))
}
