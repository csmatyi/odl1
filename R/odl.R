#' Analyzes Organelle DNA Lineages
#' @export
odl <- function() {
  # Import necessary packages
  install.packages("BiocManager")
  install.packages("msa")
  install.packages("Biostrings")
  install.packages("ape")
  install.packages("cluster")
  library(msa)
  library(Biostrings)
  library(ape)
  library(cluster)

    # read in species and accession list
  acclist <- NULL
  tryCatch(
    acclist <- as.data.frame(read.table(choose.files(filters = Filters[c("zip", "All"),]),header=T,row.names=1,sep="\t",check.names=FALSE)),
            error = function(e)
              print("Could not read species & accession list!"))

  if ((dim(acclist) < 2) | (dim(acclist)[2] != 2)) {
    print("Improper format of accession file!")
    exit()
  }

  colnames(acclist) <- c("species","accession")

  # Exit if there are species with multiple accessions
  if (length(as.list(names(which(table(acclist$species)>1)))) > 0) {
    print("Multiple accessions detected for at least one species! Please edit accession file!")
    exit()
  }

  # number of species
  n <- dim(acclist)[1]
  # height of gene order for each species
  offset <- 300
  # length of stretch of species name
  species_name_offset <- 5000
  # total height of figure
  a <- n * offset

  # list of all products: genes, tRNA, rRNA
  prod <- c()
  # maximum length of DNA
  maxlen <- 0
  # maximum number of genes
  max_n_genes <- 0

  # add sequence column to acclist
  acclist$sequence <- NA
  acclist$sequence_length <- 0

  # go through list, download GenBank annotation, package into variables
  for (acc in acclist$accession) {
    print(paste("Downloading ",acc,"...",sep=""))
    tryCatch( annot <- getAnnotationsGenBank(acc),
              error = function(e)
                print("getAnnotationsGenBank unsuccessful!"))
    # put gene/RNA products into prod
    prod <- unique(c(sort(prod),sort(annot$product)))
    # get end of last element
    alen <- annot$end[length(annot$end)]
    # if longer than longest position, overwrite it
    if (alen > maxlen) maxlen <- alen
    # get number of genes/RNAs
    glen <- length(annot$product)
    # if more than largest number, overwrite
    if (glen > max_n_genes) max_n_genes <- glen
    # sleep so there aren't too many requests all at once
    Sys.sleep(1)

    # get sequence
    GBi <- read.GenBank(acc,as.character=T)
    Sys.sleep(1)
    # the name of the column in GBi is some accession umber, so we rename it
    names(GBi) <- "id"
    seqlen <- length(GBi$id)
    seq <- paste(GBi$id,sep="",collapse="")
    acclist[acclist$accession==acc,]$sequence <- seq
    acclist[acclist$accession==acc,]$sequence_length <- seqlen
  }

  # width of figure: species name offset for species names,
  # plus longest genome plus 250 padding for sequence length
  b <- maxlen + species_name_offset + 250

  genes <- prod
  ngenes <- length(genes)
  # select ngenes number of colors from the rainbow palette
  # and put them into gene_colors
  colors <- rainbow(ngenes)
  gene_colors <- structure(names=genes,colors)

  # genemx is a matrix containing info for n species and all genes
  genemx <- matrix(nrow=n,ncol=max_n_genes)
  rownames(genemx) <- acclist$species

  # create gene order plot
  c <- 0
  jpeg("genomes_gene_order.jpg",height=a,width=b)
  par(mfrow=c(1,2))

  plot(c(1, b), c(1, a), type= "n", xlab = "", ylab = "", axes=T)
  par(mar=c(0,0,0,0))
  # go through all accessions
  for (acc in acclist$accession) {
    c <- c + 1
    species <- acclist$species[c]
    # get species' annotation
    tryCatch( annot <- getAnnotationsGenBank(acc),
              error = function(e)
                print("getAnnotationsGenBank unsuccessful!"))
    Sys.sleep(1)

    # go through all of species' genes/RNAs (products)
    # place species name at left
    text(1,offset*(c-1)+(offset/2),species,cex=11,adj=0)
    for (i in 1:length(annot$product)) {
      if (!is.na(annot$product[i])) {
        # get start and end point
        s <- annot$start[i]
        e <- annot$end[i]
        # draw product
        rect(s+species_name_offset,offset*(c-1) + 1,e+species_name_offset,offset*(c-1) + (offset-30),
             col=gene_colors[annot$product[i]])

        # add product elements to genemx
        genemx[species,i] = annot$product[i]
      }
    }
    # place sequence length at right
    seqlen <- acclist[acclist$accession==acc,]$sequence_length
    text(25+e+species_name_offset,offset*(c-1)+(offset/2),seqlen,cex=11,adj=0)
    # list of genes/RNAs for species
    sp_genes_rnas <- paste(species," |",sep="")
    # go through all products and add them to sp_genes_rnas and write them into a file
    for (gr in annot$product) {sp_genes_rnas <- paste(sp_genes_rnas,gr,sep="\t")}
    write(sp_genes_rnas, file=paste("genes_rnas.txt"), sep="\t", append = T)
  }

  # draw legend with no plot, just legend
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  par(mar=c(0,0,0,0))
  legend('topleft',legend=names(gene_colors),fill=gene_colors,cex=15,ncol=3,bty='o')
  dev.off()

  ### MATRIXES and HEATMAPS ###
  # calculate distance matrix and sequence similarity array

  # distance matrix
  dmx <- matrix(nrow = n, ncol = n)
  # sequence similarity matrix
  simmx <- matrix(nrow = n, ncol = n)
  # column and row names for both are species, since both are square matrixes
  colnames(dmx) <- acclist$species
  rownames(dmx) <- acclist$species
  colnames(simmx) <- acclist$species
  rownames(simmx) <- acclist$species

  # go through all species
  for (i in seq(1,n)) {
    # get first species (i)
    speciesi <- acclist[i,1]
    print(paste("Analyzinging ",speciesi,"...",sep=""))
    # get its accession number
    acci <- acclist[i,2]
    # get its sequence
    seqi <- acclist[acclist$accession == acci,]$sequence

    for (j in seq(i,n)) {
      # get data and sequence for speciesj
      speciesj <- acclist[j,1]
      accj <- acclist[j,2]
      # get its sequence
      seqj <- acclist[acclist$accession == accj,]$sequence

      # % similarity with percent identity = pid function
      pcnt_sim <- pid(pairwiseAlignment(seqi,seqj,type="global"),"PID2")/100
      # for pixles in diagonal
      if (i == j) {
        # distance is 0, and similarity is 100%
        dmx[speciesi,speciesj] = 0
        simmx[speciesi,speciesj] = 1
      }
      # for all other species pairs
      if (i < j) {
        simmx[speciesi,speciesj] = pcnt_sim
        simmx[speciesj,speciesi] = pcnt_sim
        # initialize distance to 0
        d <- 0
        # get gene products for speciesi, and exlude NA elements
        x <- genemx[speciesi,]
        x <- x[!is.na(x)]
        # same for speciesj
        y <- genemx[speciesj,]
        y <- y[!is.na(y)]
        # get common elements to species i and j
        xy <- intersect(x,y)
        # go through all common elements (products/genes/RNAs)
        for (k in seq(1,length(xy))) {
          # an element might be present more than once (duplicate),
          # thus take the mean index for speciesi and speciesj
          d <- d + abs(mean(which(x==xy[k])) - mean(which(y==xy[k])))
        }
        # Add mean distance of common elements and also add number of elements unique to speciesi
        # and elements unique to speciesj
        # by subtracting number of common elements from number of elements in speciesi
        # and in speciesj
        d <- d/length(xy) + (length(x) - length(xy)) + (length(y) - length(xy))
        dmx[speciesi,speciesj] <- d
        dmx[speciesj,speciesi] <- d
      }
    }
  }

  # Create similarity heatmap and write matrix
  # gene order matrix is 1 minus normalized distance matrix
  # normalize so you get a scale from 0 to 1,
  # but subtract norm.distance from 1 to get gene order similarity
  mx <- 1 - dmx/max(dmx)
  # write matrix output
  write.table(dmx, file="distance_matrix.txt", sep="\t", quote = F)
  write.table(mx, file="similarity_matrix.txt", sep="\t", quote = F)

  # elements for heatmap
  species <- row.names(mx)
  myBreaks <- c(seq(0,1,by=0.01))
  cexx <- length(species)/(4*length(species))
  ceyy <- cexx
  # normalize heatmap values to 0,1
  mx_hm <- mx
  mx2 <- (mx_hm - min(mx_hm))/(max(mx_hm) - min(mx_hm))
  mx_hm <- mx2
  clr <- colorRampPalette(c("white","yellow","orange","red"))(100)

  clusmeth="ward.D2" # ward.D ward.D2 single median average mcquitty complete centroid
  heatmap_name="gene_order_similarity.jpg"
  jpeg(filename = heatmap_name, height = 1000, width =1000, units = "px", res=300) # topo.colors(100) 5500, 5000
  h <- heatmap(mx_hm, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = clr, breaks = myBreaks, border_color=NA, na.color="white", margin = c(10,10), # gray.colors(100)
               cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth), # dendrogram="none",
               labCol=as.expression(lapply(colnames(mx), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(mx), function(a) bquote(italic(.(a))))))
  invisible(dev.off())

  # Create sequence similarity heatmap
  # write sequence similarity matrix
  write.table(simmx, file="seq_sim_matrix.txt", sep="\t", quote = F)
  species <- row.names(simmx)
  myBreaks <- c(seq(0,1,by=0.01))
  cexx <- length(species)/(4*length(species))
  ceyy <- cexx

  # normalize
  sim_mx_hm <- simmx
  simmx2 <- (sim_mx_hm - min(sim_mx_hm))/(max(sim_mx_hm) - min(sim_mx_hm))
  sim_mx_hm <- simmx2
  clr <- colorRampPalette(c("white","yellow","orange","red"))(100)

  clusmeth="ward.D2" # ward.D ward.D2 single median average mcquitty complete centroid
  ssim_heatmap_name="sequence_similarity.jpg"
  jpeg(filename = ssim_heatmap_name, height = 1000, width =1000, units = "px", res=300) # topo.colors(100) 5500, 5000
  h <- heatmap(sim_mx_hm, symkey =F, symbreaks=F, scale="none", dendrogram = F, Rowv=F, Colv=F,col = clr, breaks = myBreaks, border_color=NA, na.color="white", margin = c(10,10), # gray.colors(100)
               cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth), # dendrogram="none",
               labCol=as.expression(lapply(colnames(mx), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(mx), function(a) bquote(italic(.(a))))))
  invisible(dev.off())

  ### WRITE RESULTS and Silhouette plots ###

  # Hopkins clustering stat for gene order
  res_mx <- get_clust_tendency(mx, n = nrow(mx)-1, graph = FALSE)
  mx_hop <- res_mx$hopkins_stat
  # Hopkins clustering stat for sequence similarity
  res_simmx <- get_clust_tendency(simmx, n = nrow(simmx)-1, graph = FALSE)
  simmx_hop <- res_simmx$hopkins_stat

  # Silhouette plots
  jpeg("Silhouette_geneorder.jpg")
  mx_sil <- fviz_nbclust(mx, pam, method = "silhouette") + theme_classic()
  dev.off()

  jpeg("Silhouette_sequence.jpg")
  simmx_sil <- fviz_nbclust(simmx, pam, method = "silhouette") + theme_classic()
  dev.off()

  # max clusters sequence similarity
  k_simmx <- which(simmx_sil$data$y == max(simmx_sil$data$y))

  # get clusters
  clusmeth = "ward.D2"
  row.clusters = hclust(dist(simmx),method=clusmeth)
  ctk <- cutree(row.clusters,k=k_simmx)
  filename=paste("clusters.txt")
  write.table(ctk, file=filename, col.names=F, quote=F, sep="\t")

  # get cluster stats
  header = "cluster\tno. species\tDNA length?sd\tmin\tmean\tmax\tstdev\tp-value\tneglog"
  write(header, file=paste("stats.txt"), sep="\t", append=T)

  cluster_sizes <- table(ctk)
  non_group_cc = c()

  for (n_cluster in 1:k_simmx) {
    csize = cluster_sizes[n_cluster]
    if (csize >= 2) {
      m1 = as.matrix(mx[ctk == n_cluster,ctk == n_cluster])

      x = m1[upper.tri(m1)]
      xx = as.numeric(as.list(x))
      nn = dim(m1)[1]

      m2 = as.matrix(cbind(mx[ctk != n_cluster,ctk == n_cluster],t(mx[ctk == n_cluster,ctk != n_cluster])))
      m2b = m2[!duplicated(colnames(m2))]
      non_group_cc = c(non_group_cc, m2b)

      mean_dna_length <- sprintf("%.3f",mean(acclist[ctk == n_cluster,]$sequence_length))
      sd_dna_length <- sprintf("%.3f",sd(acclist[ctk == n_cluster,]$sequence_length))

      if (csize>= 3) {
        t = t.test(x,m2b)
      } else if (csize == 2) {
        t = t.test(m2b,mu=x)
      }
      pval = t$p.value
      nglog = -log10(pval)
      min = min(x)
      max = max(x)

      mean2 = sprintf("%.3f", mean(x))
      sd2 = sprintf("%.3f", sd(x))
      min2 = sprintf("%.3f", min)
      max2 = sprintf("%.3f", max)
      pval2 = sprintf("%.3f", pval)
      nglog2 = sprintf("%.3f", nglog)

      stats = paste(n_cluster, nn, paste(mean_dna_length,"?",sd_dna_length,sep=""), min2, mean2, max2, sd2, pval, nglog2, sep="\t")
      stats2 = gsub("\n","\t",stats)
      write(stats, file="stats.txt", sep="\t", append = T)
    }
  }

  write(date(), file="results.txt", sep="\t", append=T)
  write(paste("Number of species: ",n,sep=""), file="results.txt", append=T)
  write(paste("Hopkins clustering stat: ",round(simmx_hop,4),sep=""), file="results.txt", append=T)
  write(paste("Number of optimal clusters: ",k_simmx,sep=""), file="results.txt", append=T)
}
