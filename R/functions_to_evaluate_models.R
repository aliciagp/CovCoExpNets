
#' getPseudoBulkMatrix it creates a pseudo-bulk expression matrix from a single-cell gene expression matrix
#'
#' @param name the name of the file that contains the pseudo-bulk matrix if it has been already created
#' @param pbm_dir the directory where the pseudo-bulk matrix has been stored if it has been already created
#' @param matrix a single-cell gene expression matrix with as many rows as genes and as many columns as cells
#' @param metadata a data frame containing the metadata of each cell
#' @param cellID the name of the variable that represents the ID of each cell
#' @param covariates the covariates we want to include in the pseudo-bulk expression matrix
#' @return a list containing both the pseudo-bulk expression matrix and the covariates of the donors
#' @export
#' @examples
#'
getPseudoBulkMatrix <- function(name,
                                pbm_dir=getwd(),
                                matrix=NULL,
                                metadata=NULL,
                                cellID=NULL,
                                covariates=NULL) {

  pbm <- list.files(pbm_dir, full.names=TRUE, pattern=paste0("pseudoBulkMatrix_", name))

  if(length(pbm)>0) {
    pbm <- readRDS(pbm)
  } else {
    colnames(matrix) <- gsub("\\_.*", "", colnames(matrix)) # now the colnames are the sames of the donors
    ids <- unique(colnames(matrix))

    pseudo_bulk_matrix <- data.frame()

    for(id in ids) {# for each donor, estimate the mean expression of each gene
      df <- as.data.frame(as.matrix(matrix[, grep(id, colnames(matrix))]))
      if(is.null(ncol(df))) {
        mean <- df
      } else {
        mean <- rowMeans(df)
      }
      pseudo_bulk_matrix <- rbind(pseudo_bulk_matrix, mean)
    }

    pseudo_bulk_matrix <- as.data.frame(t(pseudo_bulk_matrix)) # Genes are the rows and columns are the samples
    colnames(pseudo_bulk_matrix) <- ids
    rownames(pseudo_bulk_matrix) <- rownames(matrix)

    # Create a metadata table at sample levels instead of cell level
    sampleCovs <- metadata
    sampleCovs[, "sampleID"] <- gsub("\\_.*", "", sampleCovs[, cellID])
    sampleCovs[, cellID] <- NULL
    sampleCovs <- sampleCovs[!duplicated(sampleCovs), ]

    stopifnot(identical(sampleCovs[, "sampleID"], colnames(pseudo_bulk_matrix )))

    # Add biological covariates
    ids <- colnames(pseudo_bulk_matrix)

    for (cov in covariates) {
      pseudo_bulk_matrix <- rbind(pseudo_bulk_matrix,
                                  as.numeric(as.character(sampleCovs[, which(colnames(sampleCovs) %in% cov)])))
    }

    len_matrix <- nrow(pseudo_bulk_matrix)
    len_covs <- length(covariates)

    rownames(pseudo_bulk_matrix)[(len_matrix-(len_covs-1)):len_matrix] <- covariates
    pseudo_bulk_matrix <- as.data.frame(t(apply(pseudo_bulk_matrix, 1, function(x) as.numeric(as.character(x))))) # all variables must be numeric
    colnames(pseudo_bulk_matrix) <- ids

    pbm <- list(pseudo_bulk_matrix=pseudo_bulk_matrix,
                sampleCovs=sampleCovs)
    saveRDS(pbm, paste0(pbm_dir, "pseudoBulkMatrix_", name, ".rds"))
  }

  return(pbm)
}



#' getBestModels it selects the top x best models of the list
#'
#' @param models_list a list of glmnet models for the same condition
#' @param top the number of models to return
#' @return a data frame containing as many rows as annotations
#' @export
#' @examples
#'
getBestModels <- function(models_list,
                          top=5) {

  summary <- models_list$models_summary
  summary <- summary[summary$nzero>=5, ]
  summary <- summary[summary$testSamples>=5, ]

  if(length(grep("R2", colnames(summary)))>0) {
    summary <- summary[order(summary$test_rmse, decreasing=F), ]
  } else {
    summary <- summary[summary$sensitivity>=0.3, ]
    summary <- summary[summary$specificity>=0.3, ]
    summary <- summary[order(summary$accuracy, decreasing=T), ]
  }

  if(nrow(summary)<top) {
    top <- nrow(summary)
  }

  bestModels <- models_list$models[as.numeric(rownames(summary)[1:top])]
  names(bestModels) <- as.numeric(rownames(summary)[1:top])

  return(bestModels)
}



#' getPredictions it returns the predictions for train, test and replication sets if available
#'
#' @param models_list a list of glmnet models for the same condition obtained with getGLMNETmodelsPerCov
#' @param models_selected the position of the model of interest in the models_list
#' @param name the name of the file that contains the pseudo-bulk matrix if it has been already created
#' @param pbm_dir the directory where the pseudo-bulk matrix has been stored if it has been already created
#' @param covariate the name of the covariate we want to predict
#' @return a data frame containing as many rows as annotations
#' @export
#' @examples
#'
getPredictions <- function(models_list,
                           model_selected,
                           name,
                           covariate=c("Brain_Bank_Path_Dx", "Donor_age_years", "Sex"),
                           pbm_dir) {


  # Define parameters
  if(length(grep("R2", colnames(models_list$models_summary)))>0) {
    type="response"
  } else {
    type="class"
  }

  # Get predictions for train and test sets
  mymodel <- models_list$models[[as.numeric(model_selected)]]
  data <- models_list$data

  index <- which(colnames(data) %in% rownames(mymodel$fit.preval))
  cov <- models_list$cov
  data_train <- t(data[, index])
  cov_train <- cov[index]
  data_test <- t(data[, -index])
  cov_test <- cov[-index]

  pred_train <- predict(mymodel, newx = data_train, type = type, s = "lambda.min")
  pred_test <- predict(mymodel, newx = data_test, type = type, s = "lambda.min")


  # Get replication data if available
  pbm <- getPseudoBulkMatrix(name, pbm_dir=pbm_dir)
  tab <- table(pbm$sampleCovs$dataset)

  if(length(tab)==2) {
    # Get pseudo-bulk matrix for the replication set
    name_repl <- names(which.min(tab))
    data_repl <- t(pbm$pseudo_bulk_matrix[, which(pbm$sampleCovs$dataset %in% name_repl)])
    data_repl <- data_repl[, -which(colnames(data_repl) %in% covariate)]

    # Get covariate for the replication set
    cov_repl <- pbm$sampleCovs[pbm$sampleCovs$dataset==name_repl, ]
    cov_repl <- as.vector(unlist(cov_repl[, covariate]))

    # Get predictions
    pred_repl <- predict(mymodel, newx = data_repl, type = type, s = "lambda.min")

  } else {
    cov_repl <- NULL
    pred_repl <- NULL
  }

  # Create a data frame with the results
  observed <- c(cov_train, cov_test, cov_repl)
  predicted <- c(as.vector(pred_train), as.vector(pred_test), as.vector(pred_repl))
  group <- c(rep("train", length(cov_train)),
             rep("test", length(cov_test)),
             rep("replication", length(cov_repl)))
  ids <- c(rownames(data_train), rownames(data_test), rownames(data_repl))

  results <- as.data.frame(cbind(observed, predicted, group, ids))
  results$group <- factor(results$group, levels=c("train", "test", "replication"))
  results$predicted <- as.numeric(results$predicted)
  results$observed <- as.numeric(results$observed)

  return(results)
}



#' plotPredictions it plots the predictions results for train and test sets (and replication set if available).
#' It uses a scatter plot if the variable is continuous or confusion matrix if the variable is categorical
#' @param models_list a list of glmnet models for the same condition obtained with getGLMNETmodelsPerCov
#' @param models_selected the position of the model of interest in the models_list
#' @param name the name of the file that contains the pseudo-bulk matrix if it has been already created
#' @param pbm_dir the directory where the pseudo-bulk matrix has been stored if it has been already created
#' @param covariate the name of the covariate we want to predict
#' @return if the variable to predict is continuous, it returns a scatter plot, ifnot, it returns a confusion matrix for each set
#' @export
#' @examples
#'
plotPredictions <- function(models_list,
                            models_selected,
                            name,
                            pbm_dir,
                            covariate) {

  # Load the libraries
  require(glmnet)
  require(caret)

  # Define parameters
  if(length(grep("R2", colnames(models_list$models_summary)))>0) {
    type="response"
  } else {
    type="class"
  }

  plots_list <- list()

  # For each model selected
  for(model in models_selected) {

    # Get predictions
    data <- getPredictions(models_list=models_list,
                           model_selected=model,
                           name=name,
                           covariate=covariate,
                           pbm_dir=pbm_dir)

    # Plot a scatter plot if it is a continuous variable
    if(type=="response") {
      p <- ggplot(data,
                   aes(x = predicted,
                       y = observed,
                       col=group)) +
        theme_classic() +
        theme(axis.title=element_text(size=14),
              axis.text=element_text(size=12),
              legend.position="bottom",
              legend.text = element_text(size=12),
              legend.title = element_text(size=12)) +
        guides(col=guide_legend(nrow=3,byrow=TRUE)) +
        geom_point(size=2) +
        ggtitle(paste0(model, " model"))

    } else { #If not, plot a confusion matrix for each set
      cm_train <- confusionMatrix(as.factor(as.vector(data$predicted[data$group=="train"])),
                                  as.factor(data$observed[data$group=="train"]), positive="1")

      cm_test <- confusionMatrix(as.factor(as.vector(data$predicted[data$group=="test"])),
                                  as.factor(data$observed[data$group=="test"]), positive="1")

      cm_repl <- confusionMatrix(as.factor(as.vector(data$predicted[data$group=="replication"])),
                                 as.factor(data$observed[data$group=="replication"]), positive="1")

      df <- as.data.frame(rbind(cm_train$table, cm_test$table, cm_repl$table))
      df$group <- rep(c("train", "test", "replication"), each=4)
      df$group <- factor(df$group, levels=c("train", "test", "replication"))

      p <- ggplot(data =  df, mapping = aes(x = Prediction, y = Reference)) +
        facet_wrap(.~group) +
        geom_tile(aes(fill = Freq), colour = "white") +
        geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1, size=6) +
        scale_fill_gradient(low = "lightblue", high = "blue") +
        theme_bw() + theme(legend.position = "none",
                           axis.text=element_text(size=12))
    }

    plots_list[[model]] <- p
  }

  return(plots_list)
}


