##### set global variables #####################################################

#define the path where to find the xlsx file with the DESeq2 results
#one sheet should be one comparison
xlsx_path <- "host_data/results/01_DESeq/data/DEG_tables/xlsx/DEG_shrink.xlsx"

#define p-value thresholds
p_thresholds <- c(0.05,0.01)

#define a single fc threshold for volcano plots
fc_cut <- 0.5

#define maximal number of go-clusters
max_cluster <- 10

#path for the resuls
res_path <- "host_data/results/02_DEG_analysis/"
go_res_path <- paste0(res_path,"go_tables/")
go_plots_path <- paste0(res_path,"go_plots/")

##### create folders ###########################################################

p_paths <- as.character(p_thresholds)
p_paths <- gsub(pattern = "0.",replacement = "",x = p_paths,fixed = T)
p_paths <- paste0(res_path,"cutoff_",p_paths,"/")
names(p_paths) <- p_thresholds
for(p in p_paths){
  dir.create(path = p)
}
dir.create(go_res_path)
dir.create(go_plots_path)

##### load libraries ###########################################################

library("readxl")
library("writexl")
library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
library("topGO")
library("biomaRt")
library("rrvgo")
source("host_data/scripts/analysis/lib/helper_functions.R")

##### make go table if necessary ###############################################

ontologies <- c("BP","MF","CC")
go_path <- "host_data/data/gene_ontology/go_table.txt"
if(!file.exists(go_path)){
  ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  at = c("ensembl_gene_id", "go_id", "namespace_1003")
  go = getBM(attributes = at, mart = ensembl)
  colnames(go)[3] <- "ontology"
  go <- go[go$go_id!="",]
  go[,3] <- c("biological_process"="BP",
              "molecular_function"="MF",
              "cellular_component"="CC")[go[,3]]
  write.table(x = go,file = go_path,quote = F,sep = "\t",col.names = T)
  
  ontology_list <- lapply(ontologies,FUN=function(ont){
    mat <- go[go[,"ontology"]==ont,,drop=F]
    uniqueGenes <- unique(mat[,"ensembl_gene_id"])
    goByGene <- lapply(uniqueGenes,FUN=function(g){
      mat[mat[,"ensembl_gene_id"]==g,"go_id"]
    })
    names(goByGene) <- uniqueGenes
    return(goByGene)
  })
  names(ontology_list) <- ontologies
  saveRDS(ontology_list,
          file = "host_data/data/gene_ontology/ontology_list.rds")
  
  rm(list = c("go","ontology_list"))
}
go_table <- read.table("host_data/data/gene_ontology/go_table.txt",
                       sep="\t",
                       header = T)
ontology_list <- readRDS(file = 
                           "host_data/data/gene_ontology/ontology_list.rds")
sim_file <- "host_data/data/gene_ontology/revigo.rds"
if(!file.exists(sim_file)){
  sim_matrix_list <- lapply(ontologies,FUN=function(ont){
    calculateSimMatrix(unique(go_table$go_id),
                       orgdb="org.Mm.eg.db",
                       ont=ont,
                       method="Rel")
  })
  names(sim_matrix_list) <- ontologies
  saveRDS(object = sim_matrix_list,
          file = sim_file)
}
sim_matrix_list <- readRDS(file = sim_file)

##### load data ################################################################

deseq_results <- xlsx_all_sheets(xlsx_path)
contrast_names <- names(deseq_results)

##### volcano plots ############################################################

for(p in p_thresholds){
  volcano_path <- paste0(p_paths[as.character(p)],"volcano/")
  dir.create(volcano_path)
  for(contrast in contrast_names){
    res <- deseq_results[[contrast]]
    EnhancedVolcano(toptable = res,
                    lab = res$Symbol,
                    x = "log2FoldChange",
                    y="padj",
                    title = contrast,
                    subtitle = "",
                    pCutoff = p,
                    FCcutoff = fc_cut)
    double_ggsave(path = paste0(volcano_path,"volcano_",contrast),
                  height = 7,
                  width = 7)
  }
}

##### MA plots #################################################################

for(p in p_thresholds){
  ma_path <- paste0(p_paths[as.character(p)],"MA/")
  dir.create(ma_path)
  for(contrast in contrast_names){
    res <- deseq_results[[contrast]]
    MAplot_gg(df = res,
              y = "log2FoldChange",
              x = "baseMean",
              p_threshold = p,
              p_column = "padj",
              neutral="darkgrey",
              main=contrast)
    double_ggsave(path = paste0(ma_path,"MA_",contrast),
                  height = 10,
                  width = 10)
  }
}


##### GO enrichment ############################################################

#prepare a df that contains the p values for each gene in each comparison
p_df <- p_df_from_deseq(deseq_list = deseq_results)
gene_names <- rownames(p_df)
go_tables <- lapply(p_thresholds,FUN=function(p){
  go_from_p_df(p_df=p_df,ontology_list=ontology_list,p=p,nodesize = 20)
})
names(go_tables) <- paste0("sig_level_",p_thresholds)

#export go results
go_export <- lapply(names(go_tables),FUN=function(topname){
  out <- unlist(go_tables[[topname]],recursive = F)
  write_xlsx(x = out,path = paste0(go_res_path,topname,".xlsx"))
})

##### REVIGO ###################################################################

#precompute a similarity matrix for terms
sim_matrix_list <- lapply(ontologies,FUN=function(ont){
  calculateSimMatrix(unique(go_tables[[1]][[1]][[ont]][,"GO.ID"]),
                     orgdb="org.Mm.eg.db",
                     ont=ont,
                     method="Rel")
})
names(sim_matrix_list) <- ontologies

#cluster terms and pick the x most significant clusters and the most significant
#term for each
reduced_terms_by_pthr <- lapply(go_tables,FUN=function(res_list_thr){
  reduced_by_set <- lapply(res_list_thr,FUN=function(res_list_set){
    reduced_by_ont <- lapply(ontologies,FUN=function(ont){
      res <- res_list_set[[ont]]
      sim_mat <- sim_matrix_list[[ont]]
      scores <- setNames(-log10(res$pvalue), res$GO.ID)
      reducedTerms <- suppressMessages(reduceSimMatrix(sim_mat,
                                                       scores,
                                                       threshold=0.5,
                                                       orgdb="org.Mm.eg.db"))
      min_p_by_cluster <- sapply(1:max(reducedTerms$cluster),FUN=function(i){
        terms <- reducedTerms[reducedTerms$cluster==i,"go"]
        minp <- min(res[terms,"pvalue"])
        return(minp)
      })
      parent_per_cluster <- sapply(1:max(reducedTerms$cluster),FUN=function(i){
        parent <- reducedTerms[which(reducedTerms$cluster==i)[1],"parent"]
        return(parent)
      })
      pruned_terms <- parent_per_cluster[head(order(min_p_by_cluster),
                                              n=max_cluster)]
      pruned_res <- res[pruned_terms,]
      return(pruned_res)
    })
    names(reduced_by_ont) <- ontologies
    return(reduced_by_ont)
  })
  names(reduced_by_set) <- colnames(p_df)
  return(reduced_by_set)
})
names(reduced_terms_by_pthr) <- names(go_tables)

#save the the tables with the reduced terms
go_export <- lapply(names(go_tables),FUN=function(topname){
  out <- unlist(reduced_terms_by_pthr[[topname]],recursive = F)
  write_xlsx(x = out,path = paste0(go_res_path,topname,"_reduced.xlsx"))
})

##### enrichment heatmaps ######################################################

#extract the go terms that are among the most significant anywhere
res <- reduced_terms_by_pthr
relevant_terms <- lapply(names(go_tables),FUN=function(topname){
  terms_by_ont <- lapply(ontologies,FUN=function(ont){
    all_terms <- lapply(res[[topname]],FUN=function(go_by_contrast){
      go_by_contrast[[ont]][,"GO.ID"]
    })
    all_terms <- unique(unlist(all_terms))
    return(all_terms)
  })
  names(terms_by_ont) <- ontologies
  return(terms_by_ont)
})
names(relevant_terms) <- names(go_tables)

#compile a df with the enrichment scores
enrichment_df <- lapply(names(go_tables),FUN=function(topname){
  df_by_ont <- lapply(ontologies,FUN=function(ont){
    terms <- relevant_terms[[topname]][[ont]]
    df <- sapply(go_tables[[topname]],FUN=function(res_by_contrast){
      res <- res_by_contrast[[ont]]
      enrichments <- pmax(res[,"Significant"],1)/res[,"Expected"]
      enrichments <- log2(enrichments)
      names(enrichments) <- res$GO.ID
      enrichments <- enrichments[terms]
      return(enrichments)
    })
    colnames(df) <- names(deseq_results)
    rownames(df) <- terms
    return(df)
  })
  names(df_by_ont) <- ontologies
  return(df_by_ont)
})
names(enrichment_df) <- names(go_tables)

#make heatmaps with the enrichment scores
go_plots <- lapply(names(enrichment_df),FUN=function(topname){
  lapply(ontologies,FUN=function(ont){
    filename <- paste0(go_plots_path,topname,"_",ont)
    df <- enrichment_df[[topname]][[ont]]
    go_heat(df=df,path = filename,mid_col = "lightgrey")
  })
})


