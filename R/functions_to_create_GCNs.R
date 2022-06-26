
#' getCorrelatedGenes it returns the top x most correlated genes for each hub gene
#'
#' @param data the pseudo-bulk expression matrix of a specific cell type
#' @param selectedGenes the names of the hub genes selected in the best model created with getGLMNETmodelsPerCov
#' @param tam the number of correlated genes we want to add to each hub gene
#' @return a data frame with as many rows as hub genes*tam
#' @export
#' @examples
#'
getCorrelatedGenes <- function(data,
                               selectedGenes,
                               tam=NULL) {

  M <- cor(as.matrix(data))
  M2 <- M[, which(colnames(M) %in% selectedGenes)]
  df <- data.frame()

  for (i in 1:ncol(M2)) {
    corr <- M2[, i]
    if(!is.null(tam)) {
      ind <- sort(abs(corr), decreasing=T)[1:(tam+1)]
    } else {
      ind <- sort(abs(corr), decreasing=T)
    }

    mydf <- data.frame(main.gene = colnames(M2)[i], genes = names(ind), cor = as.vector(ind))
    df <- rbind(df, mydf)
  }
  return(df)
}



#' getLRmodel it applies a linear/logistic regression approach and returns the R2/accuracy.
#' It is used to complete the covariate-specific GCN modules that are only made up of the hub genes at the beginning
#'
#' @param data the pseudo-bulk expression matrix of a specific cell type with as many columns as genes + one covariate
#' @return if the covariate is categorical, it returns the accuracy, if not, it returns the R2
#' @export
#' @examples
#'
getLRmodel <- function(data) {

  tab <- table(data$covariate)

  if(length(tab)==2) { # it is a categorical variable, so create a logistic regression model and report the accuracy
    model <- glm(covariate ~., data = data, family = binomial)
    probabilities <- model %>% predict(data, type = "response")
    predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
    metric <- mean(predicted.classes == data$covariate)

  } else { # it is a continuous variable, so create a linear regression model and report the R2
    mymodel = lm(covariate~ .,data=data)
    metric <- as.numeric(summary(mymodel)[9])
  }

  return(metric)
}




#' getCovGCN it returns the genes that made up each cov-specific GCN module and the corresponding annotations
#'
#' @param models_list the models list create with getGLMNETmodelsPerCov
#' @param selectedGenes the names of the hub genes selected in the best model created with getGLMNETmodelsPerCov
#' @param gcns_dir the directory where the gcn will be stored
#' @param name the name of the file to be stored
#' @return a list containing both the genes that made up each cov-specific GCN module and the corresponding annotations of each module
#' @export
#' @examples
#'
getCovGCN <- function(models_list,
                      selectedGenes,
                      name,
                      gcns_dir) {


  # Load gcn if already created
  result <- list.files(gcns_dir, full.names=TRUE, pattern=paste0(name, "_gcns"))

  if(length(result)>0) {
    result <- readRDS(result)
  } else {
    # Load the libraries
    require(gprofiler2)

    # Prepare the data
    data <- t(models_list$data)
    covariate <- models_list$cov
    cat("Number of genes selected", length(selectedGenes), "\n")

    # Get the most correlated genes for each hub gene
    df <- getCorrelatedGenes(data, selectedGenes, tam=NULL)
    main.genes <- unique(df$main.gene)
    results <- data.frame()

    # For each hub gene, create a linear/logistic regression model
    for (main.gene in main.genes) {
      mydf <- df[df$main.gene==main.gene, ]
      tam <- 2
      indx <- selectRows(colnames(data), mydf[1:tam,2])
      mydata <- data.frame(covariate = covariate, data[, indx])
      adjusted_r2 <- getLRmodel(mydata)
      r2 <- adjusted_r2
      diff <- 1

      # Add correlated genes gradually while the R2 or accuracy continues improving
      while (diff > 0){
        tam <- tam+1
        indx <- selectRows(colnames(data), df[1:tam,2])
        mydata <- data.frame(covariate = covariate, data[, indx])
        ad_r2 <- getLRmodel(mydata)
        diff <- ad_r2 -adjusted_r2
        adjusted_r2 <- ad_r2
        r2 <- c(r2, adjusted_r2)
      }

      r2=c(0, r2)

      if(diff == 0) {
        result <- cbind(mydf[1:tam,], Adjusted.R2 = r2)
      } else if (diff<0) {
        result <- cbind(mydf[1:(tam-1),], Adjusted.R2 = utils::head(r2,n=(tam-1)))
      }

      results <- rbind(results, result)
    }

    annotations <- getAnnotations(results)

    result <- list(GCN=results,
                   annotations=annotations)

    saveRDS(result, paste0(gcns_dir, "/", name, "_gcn", ".rds"))

  }

  return(result)

}


#' getAnnotations it applies a functinal enrichment analysis for the hub genes
#'
#' @param hubGenes it could be both a data frame with at least two columns, one for the main.genes and another one for the most correlated genes,
#' or it could be a vector with only the main.genes
#' @return a data frame containing as many rows as annotations
#' @export
#' @examples
#'
getAnnotations <- function(hubGenes) {

  # Annotate the network created
  all.genes <- list()

  # If we have a data frame with both the main.genes and the most correlated genes
  if(!is.null(dim(hubGenes))) {
    for(gene in unique(hubGenes$main.gene)) {
      all.genes[[gene]] <- hubGenes$genes[hubGenes$main.gene==gene]
    }
  } else { # instead we have a vector with the names of the main.genes only
    all.genes[["hubGenes"]] <- hubGenes
  }

  GO <- gprofiler2::gost(all.genes,
                         correction_method="fdr",
                         sources = c("GO","KEGG","REAC"),
                         organism = "hsapiens",
                         exclude_iea = F)

  return(GO$result)

}




#' getEnrichmentPlot it plots functional enrichment analysis results for the hub genes
#'
#' @param hubGenes it could be both a data frame with at least two columns, one for the main.genes and another one for the most correlated genes,
#' or it could be a vector with only the main.genes
#' @return a list containing both the plot and the annotations
#' @export
#' @examples
#'
getEnrichmentPlot <- function(hubGenes,
                              sources=c("GO:BP", "REAC")) {

  # Load libraries
  require(ggplot2)

  # Get annotations
  df <- getAnnotations(hubGenes)
  queriesToSelect <- names(sort(table(df$query), decreasing=T))[1:5]
  df <- df[which(df$query %in% queriesToSelect), ]

  # Filter the annotations per source
  df <- df[df$source %in% c("GO:BP", "REAC"), ]

  # Order the annotations based on how significant they are
  df$log <- -log10(df$p_value)
  df <- df[order(df$log, decreasing=T), ]

  # Select the most relevant annotation for each hub gene
  terms <- c()

  for (query in unique(df$query)) {
    terms <- c(terms, df$term_name[df$query==query][1:3])
  }

  df <- df[which(df$term_name %in% terms), ]
  df$term_name <- factor(df$term_name, levels=unique(terms))

  p <- ggplot(df, aes(x=query, y=term_name)) +
    geom_point(aes(color=log), size=2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9),
          axis.text.y = element_text(size=9),
          strip.text.y = element_text(size=11,angle=0),
          strip.text.x = element_text(size=11),
          legend.position="bottom") +
    ylab("") +
    xlab("") +
    scale_colour_gradient(low = "green", high = "red", na.value = NA, name="-log10(p-value)")

  return(list(p=p, output=df))

}



