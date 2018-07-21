peak_overlap <- function(path1,path2,heatmap_width,peak_dist){
  require(RColorBrewer)
  require(ChIPpeakAnno)
  require(rtracklayer)
  
  peaks1.idr <- readRDS(path1)
  peaks2.idr <- readRDS(path2)
  
  if (peak_dist != 0) {
    #Define peaks1 center
    features <- peaks1.idr
    wid <- width(features)
    feature.recentered <- feature.center <- features
    start(feature.center) <- start(features) + floor(wid/2)
    width(feature.center) <- 1
    start(feature.recentered) <- start(feature.center) - peak_dist
    end(feature.recentered) <- end(feature.center) + peak_dist
    
    peaks1.idr.recentered <- feature.recentered
    
    #Define peaks2 center
    features <- peaks2.idr
    wid <- width(features)
    feature.recentered <- feature.center <- features
    start(feature.center) <- start(features) + floor(wid/2)
    width(feature.center) <- 1
    start(feature.recentered) <- start(feature.center) - peak_dist
    end(feature.recentered) <- end(feature.center) + peak_dist
    
    peaks2.idr.recentered <- feature.recentered
    #Overlap of recentered peaks
    ol <- findOverlapsOfPeaks(peaks1.idr.recentered,peaks2.idr.recentered) 
  } else if (peak_dist == 0){
    #Overlap Sall4/Foxd3
    ol <- findOverlapsOfPeaks(peaks1.idr,peaks2.idr) 
  } else {
    print("Error. Wrong print_dist value.")
  }
  #Make Venn diagram of overlap
  pdf(file.path(dirname(path1),"venn.pdf"))
  makeVennDiagram(ol,totalTest = 1.87e+9 * 0.03 / (2 * 150),NameOfPeaks = c(basename(dirname(path1)),basename(dirname(path2))))
  dev.off()
  
  #Heatmap and density plot of overlapping peak regions
  
  features <- ol$peaklist[[length(ol$peaklist)]]
  
  wid <- width(features)
  feature.recentered <- feature.center <- features
  start(feature.center) <- start(features) + floor(wid/2)
  width(feature.center) <- 1
  start(feature.recentered) <- start(feature.center) - heatmap_width/2
  end(feature.recentered) <- end(feature.center) + heatmap_width/2
  
  ## Consider importData function in bioconductor trackViewer package 
  ## to import the coverage.
  ## compared to rtracklayer, it saves time when handling huge datasets.
  
  #Define .bw file locations
  ##File path to experimental .bw files
  bws1 <- dir(file.path(dirname(path1),"/bw"),"bw")
  bws2 <- dir(file.path(dirname(path2),"/bw"),"bw")
  
  ##Import .bw files
  cvglists1 <- sapply(file.path(dirname(path1),"/bw", bws1), import, 
                        format="BigWig", 
                        which=feature.recentered, 
                        as="RleList")
  
  cvglists2 <- sapply(file.path(dirname(path2),"/bw",bws2), import, 
                        format="BigWig", 
                        which=feature.recentered, 
                        as="RleList")

  names(cvglists1) <- gsub(".dedup.sorted.bw", "", bws1)
  names(cvglists2) <- gsub(".dedup.sorted.bw", "", bws2)
  
  sig1 <- featureAlignedSignal(cvglists1, feature.center, 
                               upstream=heatmap_width/2, downstream=heatmap_width/2) 
  
  sig2 <- featureAlignedSignal(cvglists2, feature.center, 
                               upstream=heatmap_width/2, downstream=heatmap_width/2)

  pdf(file.path(dirname(path1),"peaks1_heatmap.pdf"))
  heatmap1 <- featureAlignedHeatmap(sig1, feature.center, 
                                    upstream=heatmap_width/2, downstream=heatmap_width/2,
                                    upper.extreme=c(2,2,2,2,2), color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
  dev.off()
  
  pdf(file.path(dirname(path2),"peaks2_heatmap.pdf"))
  heatmap2 <- featureAlignedHeatmap(sig2, feature.center, 
                                    upstream=heatmap_width/2, downstream=heatmap_width/2,
                                    upper.extreme=c(2,2,2,2,2), color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
  dev.off()
  
  pdf(file.path(dirname(path2),"peaks1_distribution.pdf"))
  featureAlignedDistribution(sig1, feature.center, 
                             upstream=heatmap_width/2, downstream=heatmap_width/2,
                             type="l")
  dev.off()
  
  pdf(file.path(dirname(path2),"peaks2_distribution.pdf"))
  featureAlignedDistribution(sig2, feature.center, 
                             upstream=heatmap_width/2, downstream=heatmap_width/2,
                             type="l")
  dev.off()
}
