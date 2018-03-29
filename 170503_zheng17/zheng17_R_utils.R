#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

# ----------------------
# load (pure types)
# ----------------------
load_purified_pbmc_types<-function(pure_select_file,pure_gene_ens) {
  pure_select_id <- pure_select_file$pure_id   # from pure_select_file
  pure_select_avg <- pure_select_file$pure_avg # from pure_select_file
  pure_use_genes <- pure_select_file$pure_use_genes # from pure_select_file
  pure_use_genes_ens<-pure_gene_ens[pure_use_genes]
  avg<-data.frame(t(pure_select_avg))
  rownames(avg)<-pure_use_genes_ens
  names(avg)<-pure_select_id
  return(avg)
}
# -----------------------------------
# train a multinomial (averging)
# -----------------------------------
.train_multinomial <- function(x) {
  (1+colSums(x)) / sum(1+colSums(x))
}
# -----------------------------------
# train a multinomial (averging)
# -----------------------------------
.compare_by_cor<-function(m_filt,use_gene_ids,dmap_data) {
  
  sig_genes <- intersect(use_gene_ids, rownames(dmap_data))
  m_forsig <- as.matrix(m_filt[,which(use_gene_ids %in% sig_genes)])
  sig_data_filt <- dmap_data[match(use_gene_ids[which(use_gene_ids %in% sig_genes)], rownames(dmap_data)),]
  
  z <- lapply(1:ncol(sig_data_filt), function(j) sapply(1:nrow(m_forsig), function(i) cor(m_forsig[i,], sig_data_filt[,j], method='spearman')))
  z <- do.call(cbind, z)
  colnames(z) <- colnames(sig_data_filt)
  z
}
# -------------------------------------
# normalize expression of one gene
# -------------------------------------
# v - expression vector of one gene
.normalize_by_gene<-function(v) {
  m <- log(1+v)
  m<-m-mean(m)
  m<-m/sd(m)
}
# -------------------------------
# get cluster specific genes
# -------------------------------
.get_cluster_specific_genes<-function(km_avg,use_genes,clus_id,n=10) {
  gene_scores <- do.call(rbind,lapply(1:nrow(km_avg), function(i) {
    data.frame(subclu=i,
               gene=use_genes,
               abs_score=colMeans(abs(sweep(log2(km_avg[-i,]), 2, log2(km_avg[i,]), '-'))),
               pos_score=colMeans(sweep(log2(km_avg[-i,]), 2, log2(km_avg[i,]), '-')),
               stringsAsFactors=F)
  } ) )
  subclu_specific_genes_pos_score <- unique((gene_scores %>% group_by(subclu) %>% arrange(pos_score) %>% top_n(n, -pos_score) %>% ungroup())$gene)
  
  subclu_specific_gene_data <- km_avg[,match(subclu_specific_genes_pos_score, use_genes)]
  colnames(subclu_specific_gene_data) <- subclu_specific_genes_pos_score
  rownames(subclu_specific_gene_data)<-clus_id
  
  sub_nor<-scale(subclu_specific_gene_data)
  list(sub_nor=sub_nor,gene_scores=gene_scores)
}
# ------------------------------------------------
# save figures to file with specified format
# ------------------------------------------------
.save_figure <- function(ggp,dir,id,marker,width,height) {
  # create file name
  x <- id
  x <- gsub('+','',x,fixed = TRUE); x <- gsub('-','',x,fixed = TRUE)
  x <- gsub(' ','_',x,fixed = TRUE); x <- gsub('/','_',x,fixed = TRUE)  
  filepath <- file.path(dir,paste(x,'_',marker,'.png',sep=""))
  ggsave(file=filepath,ggp,width=width,height=height)
  cat('Saved figure to: ',filepath,'\n')
}
# ------------------------------------------------
# plot gene expression based on gene markers
# ------------------------------------------------
.plot_gene_expression <- function(markers, genes, mat, tsne_df, style, title, dir=NULL) {
  for (marker in markers) {
    idx<-match(marker,genes)
    goi<-.normalize_by_gene(mat[,idx])
    p <- ggplot(tsne_df,aes(X1,X2,col=goi))+geom_point(size=0)+ggtitle(paste(title,'(',marker,')')) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    midpnt <- 0.5
    if (style == 'type1') {
      min_lim <- min(goi); max_lim <- max(goi)
      if (marker == 'CD8A') {midpnt <- 0; min_lim <- -1}      
      p <- p + scale_colour_gradient2(low="royalblue3",high="red",mid="white",midpoint=midpnt,limits=c(min_lim,max_lim),
                                      guide=guide_colorbar(barwidth = 0.8, barheight = 5)) 
    } else if (style == 'type2') {
      if (marker == 'FCER1A') {midpnt <- 1}
      if (marker == 'S100A8') {midpnt <- 0.3}
      p <- p + scale_colour_gradient2(low="olivedrab1",high="red",mid="white",midpoint=midpnt, 
                                      guide=guide_colorbar(barwidth = 0.8, barheight = 5))
    }
    if (is.null(dir)) {
      print(p)
    } else {
      .save_figure(p,dir,title,marker,width=4.0,height=3.5)
    }
  }
}
# ---------------------------------------------------------
# plot cluster heatmap and save cluster specific genes
# ---------------------------------------------------------
.plot_heatmap_and_write_genes <- function(n_clust,tsne_df,mat,genes,dir) {
  clu<-c(1:n_clust)
  km_centers <- do.call(cbind, lapply(1:n_clust,function(i) .train_multinomial(mat[which(tsne_df$k == clu[i]),])))
  km_avg<-t(km_centers)
  sub_nor<-.get_cluster_specific_genes(km_avg,genes,clu,10)
  write.table(sub_nor$gene_scores,file = dir,sep="\t",quote=F,row.names = F)
  pheatmap(sub_nor$sub_nor,cluster_rows = T,cluster_cols = T,fontsize_col=8,color=grey.colors(8))
}
# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}
# -----------------------------------------------------
# compute top n principle components with propack 
# -----------------------------------------------------
.do_propack <- function(x,n) {
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  m <- sweep(m, 2, apply(m, 2, sd), '/')
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(ppk=ppk,pca=pca, m=m,use_genes=use_genes)
}
# ---------------------------------------------
# normalize the gene barcode matrix by umi
# ---------------------------------------------
.normalize_by_umi <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=x_use_genes)
}
# --------------------------------------
# down-sample a gene-barcode matrix
# --------------------------------------
.downsample_gene_bc_mtx <- function(json, orig_matrix_data, mol_info, tgt_rpc, target_type='raw_reads', transcriptome='hg19') {
  if (!(target_type %in% c('raw_reads', 'conf_mapped_reads'))) {
    stop(sprintf('Unsupported target_type: %s', target_type))
  }
  if (length(orig_matrix_data) > 2) {
    warning('Multiple transcriptomes not yet implemented')
  }
  cat("Filtering mol info\n")
  mol_info <- mol_info[reads > 0]
  cat("Sorting mol info\n")
  setkey(mol_info, barcode, gene, umi)
  cat("Aggregating mol info\n")
  bc_gene_umi <- mol_info[, j=list(reads=sum(reads)), by=c('barcode', 'gene', 'umi')]
  n_cells <- json[[sprintf("%s_filtered_bcs", transcriptome)]]
  tot_reads <- json$total_reads
  candidate_reads <- sum(bc_gene_umi$reads)
  candidate_read_frac <- candidate_reads / json$total_reads
  orig_mat_barcodes <- sub('-.*$', '', orig_matrix_data[[transcriptome]]$barcodes)
  subsampled_mats <- lapply(tgt_rpc, function(tgt_rpc_i) {
    cat('.')
    if (target_type == 'raw_reads') {
      tgt_candidate_reads <- tgt_rpc_i * n_cells * candidate_read_frac
      candidate_rpc <- tgt_candidate_reads / n_cells
      raw_reads_per_cell <- tgt_rpc_i
    } else if (target_type == 'conf_mapped_reads') {
      tgt_candidate_reads <- tgt_rpc_i * n_cells
      candidate_rpc <- tgt_candidate_reads / n_cells
      raw_reads_per_cell <- candidate_rpc / candidate_read_frac
    }
    if (tgt_candidate_reads > candidate_reads) {
      return(NA)
    }
    subsample_rate <- tgt_candidate_reads / candidate_reads
    cat("Subsampling\n")
    bc_gene_umi_subsampled <- bc_gene_umi %>% mutate(reads=rbinom(length(reads), reads, subsample_rate))
    cat("Sorting\n")
    setkey(bc_gene_umi_subsampled, barcode, gene)
    cat("Re-aggregating\n")
    bc_gene_counts <- bc_gene_umi_subsampled[barcode %in% orig_mat_barcodes, j=list(count=sum(reads > 0)), by=c('barcode', 'gene')]
    cat("Building matrix\n")
    with(bc_gene_counts, sparseMatrix(i = match(barcode, orig_mat_barcodes),
                                      j = 1 + gene, x = count, dims=dim(orig_matrix_data[[transcriptome]]$mat)))
  } )
  return(subsampled_mats)
}
# ------------------------------------------
# reassign purified ids, due to overlap
# ------------------------------------------
.reassign_pbmc_11<-function(z) {
  unlist(lapply(1:nrow(z),function(i) {
    best<-which.max(z[i,])
    x<-z[i,]
    nextbest<-which.max(x[x!=max(x)])
    # if best is CD4+ T helper, and the next best is cd4+/cd25+, or cd4+/cd45ro+ or cd4+/cd45ra+/cd25-, use the next best assignment
    if (best==9 & (nextbest==3 || nextbest==4 || nextbest==6)) {
      best=nextbest
    }
    # if best is CD8+, and the next best is CD8+/CD45RA+, use next best assignment
    if (best==7 & nextbest==5) {
      best=5
    }
    best
  }))
}
# ----------------------------------------
# set color 68k colors for aesthetics
# ----------------------------------------
.set_pbmc_color_11<-function() {
  myColors <- c( "dodgerblue2",
                      "green4", 
                      "#6A3D9A", # purple
                       "grey",
                       "tan4",
                       "yellow", 
                      "#FF7F00", # orange
                      "black",
                      "#FB9A99", # pink
                      "orchid",
                      "red")
  id<-c("CD19+ B","CD14+ Monocyte","Dendritic","CD56+ NK","CD34+","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T","CD4+/CD45RO+ Memory","CD4+ T Helper2","CD8+/CD45RA+ Naive Cytotoxic","CD8+ Cytotoxic T")
  names(myColors)<-id
  scale_colour_manual(name = "cls",values = myColors)
}
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2",
         "#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
