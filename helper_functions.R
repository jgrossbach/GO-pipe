cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#function that automatically loads all sheets for an xlsx file with the readxl
#package. The resulting object is a list with names corresponding to the sheets
xlsx_all_sheets <- function(path){
  require(readxl)
  sheets <- excel_sheets(path)
  out <- lapply(sheets,FUN=function(sheet){
    df <- read_xlsx(path,sheet = sheet)
    df <- as.data.frame(df)
    return(df)
  })
  names(out) <- sheets
  return(out)
}


#function that makes pdf and png files with similar dimensions for the last plot
#only works for ggplot objects
double_ggsave <- function(path,...){
  ggsave(filename = paste0(path,".png"),...)
  ggsave(filename = paste0(path,".pdf"),...)
}

#function that makes a base R MA plot
MAplot <- function(df,
                   x,
                   y,
                   p_threshold,
                   p_column,
                   col_up="red",
                   col_down="blue",
                   neutral="grey",
                   xlab="log2(mean expression+1)",
                   ylab="log2FC",
                   main="",
                   cex=2){
  set.seed(1337)
  df <- df[sample(1:nrow(df)),]
  fc <- df[,y]
  me <- log2(df[,x]+1)
  p <- df[,p_column]
  sig <- p<=p_threshold
  up <- fc>0
  down <- fc<0
  sig_up <- which(sig&up)
  sig_down <- which(sig&down)
  col_vec <- rep(neutral,nrow(df))
  col_vec[sig_up] <- col_up
  col_vec[sig_down] <- col_down
  plot(NULL,xlim=range(me),ylim=range(fc),xlab=xlab,ylab=ylab,main=main,cex=cex)
  plot_df <- data.frame(x=me,y=fc,pch=19,col=col_vec)
  plot_df <- plot_df[rep(1:nrow(plot_df),each=2),]
  even <- 1:nrow(plot_df)%%2==0
  plot_df[even,"pch"] <- 1
  plot_df[even,"col"] <- "black"
  points(x = plot_df$x,y = plot_df$y,pch=plot_df$pch,col=plot_df$col)
}

#function that makes a ggplot2 MA plot
MAplot_gg <- function(df,
                      x,
                      y,
                      p_threshold,
                      p_column,
                      col_up="red",
                      col_down="blue",
                      neutral="grey",
                      xlab="log2(mean expression+1)",
                      ylab="log2FC",
                      main=""){
  set.seed(1337)
  df <- df[sample(1:nrow(df)),]
  fc <- df[,y]
  me <- log2(df[,x]+1)
  p <- df[,p_column]
  sig <- p<=p_threshold
  up <- fc>0
  down <- fc<0
  sig_up <- which(sig&up)
  sig_down <- which(sig&down)
  state <- rep("not sig",nrow(df))
  state[sig_up] <- "up"
  state[sig_down] <- "down"
  plot_df <- data.frame(x=me,y=fc,fill=as.factor(state))
  ggplot() + 
    theme_grey(base_size = 22) +
    geom_point(data = plot_df, shape=21, aes(x = x, y = y,fill=fill),size=4) +
    scale_fill_manual(values = c("up"=col_up,
                                 "down"=col_down,
                                 "not_sig"=neutral)) + 
    labs(title = main) +
    ylab("log2FC") + 
    xlab("log2(mean expression + 1)")
}

#function that takes a vector of pvalues and returns a factorized boolean 
#vector. One can specify a significance threshold and/or a minimum number of 
#hits.
is_sig <- function(pvec,sigThreshold=0.05,minHits=0){
  hits <- pvec<sigThreshold
  if(sum(hits)<minHits){
    hits <- pvec<=sort(pvec)[minHits]
  }
  hits <- as.factor(as.numeric(hits))
  names(hits) <- names(pvec)
  return(hits)
}

#function that takes a list of deseq2 results and extracts a df with adjusted
#pvalues in the same order. NAs are converted to 1s.
p_df_from_deseq <- function(deseq_list){
  geneNames <- sort(deseq_results[[1]][,"geneNames"])
  out <- sapply(deseq_list,FUN=function(x){
    pvec <- x[,"padj"]
    names(pvec) <- x[,"geneNames"]
    pvec[is.na(pvec)] <- 1
    #in case the tables are preorderer we need to reorder
    pvec <- pvec[geneNames] 
    return(pvec)
  })
  return(out)
}

#function that takes a df with p values and performs go enrichment for each
#column using fisher tests. algorithm can be specified
go_from_p_df <- function(p_df,ontology_list,p=0.05,algo="classic",nodesize=10){
  ontologies <- names(ontology_list)
  res_by_exp <- lapply(1:ncol(p_df),FUN=function(i){
    hits <- is_sig(p_df[,i],sigThreshold = p)
    res_by_ont <- lapply(ontologies,FUN=function(ont){
      GOdata <- suppressMessages(new("topGOdata",
                                     ontology = ont,
                                     allGenes = hits,
                                     annot = annFUN.gene2GO,
                                     gene2GO = ontology_list[[ont]],
                                     nodeSize=nodesize))
      res <- suppressMessages(runTest(GOdata,
                                      algorithm = algo,
                                      statistic = "fisher"))
      resTable <- GenTable(object=GOdata,
                           pvalue = res,
                           topNodes=length(res@score))
      resTable$pvalue <- res@score[resTable$GO.ID]
      rownames(resTable) <- resTable[,1]
      return(resTable)
    })
    names(res_by_ont) <- ontologies
    return(res_by_ont)
  })
  names(res_by_exp) <- colnames(p_df)
  return(res_by_exp)
}

#function to make heatmaps based on go enrichment scores
go_heat <- function(df,
                    low_col="blue",
                    high_col="red",
                    mid_col="grey",
                    border_col="black",
                    translate_terms=TRUE,
                    cluster=TRUE,
                    path=".",
                    height=10,
                    width=10){
  #prepare df to plot
  df_plot <- expand.grid(terms=rownames(df),contrast=colnames(df))
  df_plot <- cbind(df_plot,enrichment=df[as.matrix(df_plot)])
  if(cluster){
    #cluster terms
    term_order <- hclust(dist(df))$order
    df_plot$terms <- factor(df_plot$terms, levels = rownames(df)[term_order])
  }
  if(translate_terms){
    require("GO.db")
    term_levels <- Term(levels(df_plot$terms))
    df_plot$terms <- factor(Term(as.character(df_plot$terms)),
                            levels = term_levels)
  }
  
  #make the plot
  ggplot(df_plot, aes(contrast, terms, fill= enrichment)) + 
    geom_tile(color = border_col,linewidth=0.5) +
    scale_fill_gradient2(low=low_col,mid=mid_col, high=high_col)
  #save the plot
  double_ggsave(path = path,
                height = height,
                width = width)
}