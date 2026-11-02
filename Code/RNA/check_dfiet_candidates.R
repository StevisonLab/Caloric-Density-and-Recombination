# Standalone companion to Rna_diet (9)(1).R. Run from your project root.
# Uses the ORIGINAL genotype/day DEG tables, never the filtered *_sig objects.
# No packages needed. This is an exploratory literature-guided screen.
# IMPORTANT: if *_deg.xls files contain only DEGs, missing genes are UNKNOWN,
# not unchanged. Prefer complete DESeq2 results with every tested gene.
setwd("~/GitHub/Caloric-Density-and-Recombination")
deg_path <- "Raw Data/RNA/DEG"
out_path <- "Summaries and Outputs/RNA/Diet_candidates"
p_cutoff <- 0.05
lfc_cutoff <- 1 # strict >1, matching the supplied script
# Optional annotation CSV with columns gene_id,gene_name; needed if names absent.
annotation_file <- "" 
# Optional harmonized Novak 2026 results CSV; see README for required format.
novak_file <- ""
dir.create(out_path, recursive=TRUE, showWarnings=FALSE)

sets <- list(
  Novak2026_FigS12_verified_subset=c("armi","aub","cuff","egg","zuc","mael","piwi","spn-E","AGO3","Hsp83","Su(var)205"),
  Brookheart2017_redox=c("Sod1","Sod2","GstD1"),
  Brookheart2017_mitochondria=c("srl","delg","park","mtTFB1"),
  AbulBasar2019_nutrient_oogenesis="srl",
  Wang2026_ovarian_IIS_CrebA=c("CrebA","InR","foxo"),
  Kang2017_amino_acid_response=c("Thor","Gcn2","crc"),
  Zid2009_translation_mitochondria="Thor",
  Tettweiler2005_starvation_redox=c("Thor","foxo"),
  PritchettMcCall2012_IIS_TOR=c("InR","chico","Tor","S6k","Dcp-1"),
  Meiotic_context_not_diet_DEG_claim=c("mei-W68","mei-P22","mei-217","mei-218","rec","mei-9","Blm")
)
# Only explicit aliases: no substring matching. Verify against your annotation release.
aliases <- c("spargel"="srl", "dpgc-1"="srl", "sod"="sod1",
             "gst-d1"="gstd1", "mt-tfb1"="mttfb1", "parkin"="park",
             "4e-bp"="thor", "4ebp"="thor", "d4e-bp"="thor",
             "atf4"="crc", "datf4"="crc", "dfoxo"="foxo")
key <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[is.na(x)] <- NA_character_
  m <- match(y,names(aliases)); use <- !is.na(m)
  y[use] <- unname(aliases[m[use]])
  y
}
catalog <- unique(do.call(rbind,lapply(names(sets),function(s)
  data.frame(candidate=sets[[s]],evidence_group=s,stringsAsFactors=FALSE))))
catalog$key <- key(catalog$candidate)
write.csv(catalog,file.path(out_path,"candidate_catalog.csv"),row.names=FALSE)
annotation <- NULL
if(nzchar(annotation_file)) {
  annotation <- read.csv(annotation_file,stringsAsFactors=FALSE)
  stopifnot(all(c("gene_id","gene_name") %in% names(annotation)))
  if(anyDuplicated(annotation$gene_id)) stop("Annotation must have one row per gene_id.")
}
standardize <- function(d,label) {
  req <- c("gene_id","log2FoldChange","padj")
  if(!all(req %in% names(d))) stop(label,": missing ",paste(setdiff(req,names(d)),collapse=", "))
  if(anyNA(d$gene_id) || any(!nzchar(trimws(d$gene_id)))) stop(label,": blank gene_id")
  if(anyDuplicated(d$gene_id)) stop(label,": duplicate gene_id; use gene-level tables.")
  d$gene_id <- as.character(d$gene_id)
  if(!"gene_name" %in% names(d)) d$gene_name <- NA_character_
  d$gene_name <- as.character(d$gene_name)
  blank <- is.na(d$gene_name) | !nzchar(trimws(d$gene_name))
  if(!is.null(annotation)) {
    d$gene_name[blank] <- annotation$gene_name[match(d$gene_id[blank],annotation$gene_id)]
  }
  blank <- is.na(d$gene_name) | !nzchar(trimws(d$gene_name))
  # Some files put symbols directly in gene_id. Do not treat FBgn IDs as symbols.
  use <- blank & !grepl("^FBgn[0-9]+$",d$gene_id)
  d$gene_name[use] <- d$gene_id[use]
  for(v in c("log2FoldChange","padj")) {
    old <- d[[v]]; d[[v]] <- suppressWarnings(as.numeric(as.character(old)))
    if(any(!is.na(old) & nzchar(as.character(old)) & is.na(d[[v]])))
      stop(label,": nonnumeric ",v)
  }
  if(any(!is.na(d$padj) & (d$padj<0 | d$padj>1))) stop(label,": invalid padj")
  d$key <- key(d$gene_name)
  d$FDR_significant <- !is.na(d$padj) & d$padj<p_cutoff
  d$passes_original_filter <- d$FDR_significant & is.finite(d$log2FoldChange) & abs(d$log2FoldChange)>lfc_cutoff
  d$status <- ifelse(is.na(d$padj),"padj unavailable",
                     ifelse(!is.finite(d$log2FoldChange),"effect unavailable",
                            ifelse(d$passes_original_filter,"FDR and effect cutoff",
                                   ifelse(d$FDR_significant,"FDR only","not FDR significant"))))
  d
}
all_results <- list(); audits <- list(); i <- 0
for(g in c("S42","S399")) for(day in c("D2","D6")) for(pair in list(c("HC","C"),c("LC","C"),c("LC","HC"))) {
  i <- i+1
  file <- paste0(g,"_",pair[1],"_",day,"vs",g,"_",pair[2],"_",day,"_deg.xls")
  path <- file.path(deg_path,file)
  if(!file.exists(path)) stop("Missing input: ",path)
  d <- standardize(read.delim(path,check.names=FALSE,quote="",stringsAsFactors=FALSE),file)
  d$genotype <- g; d$day <- day
  d$contrast <- paste(pair,collapse="_vs_")
  # Positive log2FC means higher in the FIRST diet, if file names encode numerator.
  d$direction <- ifelse(!is.finite(d$log2FoldChange),NA_character_,
                        ifelse(d$log2FoldChange>0,paste0("higher in ",pair[1]),
                               ifelse(d$log2FoldChange<0,paste0("higher in ",pair[2]),"zero estimate")))
  all_results[[i]] <- d[,c("gene_id","gene_name","key","genotype","day","contrast",
                           "log2FoldChange","padj","FDR_significant","passes_original_filter","status","direction")]
  audits[[i]] <- data.frame(file=file,n_rows=nrow(d),n_FDR=sum(d$FDR_significant),
                            n_original_filter=sum(d$passes_original_filter),n_without_symbol=sum(is.na(d$key)|!nzchar(d$key)))
}
all <- do.call(rbind,all_results)
audit <- do.call(rbind,audits)
# Full supplied results enable selection beyond the literature candidate list.
write.csv(all,file.path(out_path,"all_gene_results.csv"),row.names=FALSE)
write.csv(audit,file.path(out_path,"input_audit.csv"),row.names=FALSE)
if(all(is.na(all$key)|!nzchar(all$key))) stop("No gene symbols: supply annotation_file with gene_id,gene_name.")

# An explicit grid ensures absent candidates are exported rather than silently lost.
conditions <- unique(all[,c("genotype","day","contrast")])
grid <- merge(catalog,conditions,by=NULL)
hits <- merge(grid,all,by=c("key","genotype","day","contrast"),all.x=TRUE,sort=FALSE)
missing <- is.na(hits$gene_id)
hits$status[missing] <- "not matched in supplied table; expression/testing unknown"
hits <- hits[order(hits$evidence_group,hits$candidate,hits$genotype,hits$day,hits$contrast),]
write.csv(hits,file.path(out_path,"all_candidate_results.csv"),row.names=FALSE)
write.csv(hits[!is.na(hits$passes_original_filter)&hits$passes_original_filter,],
          file.path(out_path,"candidate_hits_original_filter.csv"),row.names=FALSE)
# Primary comparison to Novak: LC vs HC, separated by genotype and day.
write.csv(hits[hits$contrast=="LC_vs_HC",],file.path(out_path,"candidates_LC_vs_HC.csv"),row.names=FALSE)
# Optional BASE-R heatmap: actual log2FC (not expression or pathway activity), no clustering.
u <- unique(hits[,c("candidate","genotype","day","contrast","log2FoldChange","padj")])
labels <- paste(conditions$genotype,conditions$day,conditions$contrast,sep=" | ")
genes <- unique(catalog$candidate)
mat <- matrix(NA_real_,length(genes),length(labels),dimnames=list(genes,labels))
stars <- mat
for(j in seq_len(nrow(u))) {
  rr <- match(u$candidate[j],genes)
  cc <- match(paste(u$genotype[j],u$day[j],u$contrast[j],sep=" | "),labels)
  if(!is.na(rr)&&!is.na(cc)) {
    mat[rr,cc] <- u$log2FoldChange[j]
    stars[rr,cc] <- as.numeric(!is.na(u$padj[j]) && u$padj[j]<p_cutoff)
  }
}
finite <- mat[is.finite(mat)]
if(length(finite)) {
  lim <- max(abs(finite),1); plotmat <- mat; plotmat[!is.finite(plotmat)] <- NA_real_
  pdf(file.path(out_path,"candidate_log2FC_heatmap.pdf"),width=13,height=12)
  par(mar=c(11,10,3,2))
  image(seq_along(labels),seq_along(genes),t(plotmat),
        col=colorRampPalette(c("#2166AC","white","#B2182B"))(101),
        zlim=c(-lim,lim),axes=FALSE,xlab="",ylab="",
        main=paste0("Candidate log2FC; blue negative, red positive; scale +/-",round(lim,2)))
  axis(1,seq_along(labels),labels,las=2,cex.axis=.65)
  axis(2,seq_along(genes),genes,las=2,cex.axis=.7)
  for(r in seq_along(genes)) for(c in seq_along(labels)) {
    if(is.na(plotmat[r,c])) text(c,r,"NA",cex=.4,col="grey50")
    else if(!is.na(stars[r,c])&&stars[r,c]==1) text(c,r,"*",cex=.8)
  }
  mtext("* padj < 0.05; NA = missing/unavailable; colors are fold changes, not absolute expression",side=3,cex=.65)
  dev.off()
}

# Optional direct gene-by-gene comparison, requiring FULL harmonized Novak data.
# Required CSV: gene_id,gene_name,genotype,contrast,log2FoldChange,padj
# genotype S42 or S217; contrast LC_vs_HC; log2FC MUST mean LC / HC.
# Only S42 is shared. Do not label S399 and S217 as the same genetic background.
if(nzchar(novak_file)) {
  n <- read.csv(novak_file,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(all(c("genotype","contrast") %in% names(n)))
  n <- n[n$genotype=="S42" & n$contrast=="LC_vs_HC",]
  if(!nrow(n)) stop("Novak input needs S42 LC_vs_HC rows.")
  n <- standardize(n,"Novak S42 LC_vs_HC")
  n <- n[,c("gene_id","gene_name","log2FoldChange","padj")]
  names(n)[-1] <- paste0("Novak_",names(n)[-1])
  our <- all[all$genotype=="S42" & all$contrast=="LC_vs_HC",]
  both <- merge(our,n,by="gene_id",all.x=TRUE)
  both$shared_FDR <- !is.na(both$padj)&!is.na(both$Novak_padj)&both$padj<p_cutoff&both$Novak_padj<p_cutoff
  both$same_direction <- ifelse(is.finite(both$log2FoldChange)&is.finite(both$Novak_log2FoldChange),
                                sign(both$log2FoldChange)==sign(both$Novak_log2FoldChange),NA)
  write.csv(both,file.path(out_path,"direct_Novak2026_S42_comparison.csv"),row.names=FALSE)
} else message("Direct Novak DEG overlap not run: no harmonized Novak results supplied.")
message("Done: ",normalizePath(out_path))
message("Check input_audit.csv and verify contrast directions before interpretation.")
