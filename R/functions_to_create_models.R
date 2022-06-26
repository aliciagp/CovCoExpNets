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



#' getPCAplot it plots donors in a pca space to
#' i) check for outliers as part of the data processing
#' ii) check if donors seggregate per covariate using the hub genes selected with glmnet
#' @param name the name of the file that contains the pseudo-bulk matrix if it has been already created
#' @param pbm_dir the directory where the pseudo-bulk matrix has been stored if it has been already created
#' @param matrix a single-cell gene expression matrix with as many rows as genes and as many columns as cells
#' @param metadata a data frame containing the metadata of each cell
#' @param cellID the name of the variable that represents the ID of each cell
#' @param covariates the name of the covariates we want to use to color the donors in the plot
#' @param label if TRUE, the ID of the donors will be shown
#' @param loadings if TRUE, loadings will be shown
#' @param label.repel only show the labels of the outliers
#' @param size the size of the points representing the donors
#' @param label.size if label=TRUE, the size of the labels
#' @param selectedGenes a vector of the names of the genes to be used for the PCA
#' @return a list containing all the plots generated
#' @export
#' @examples
#'
getPCAplot <- function(name,
                       pbm_dir,
                       matrix,
                       metadata,
                       cellID,
                       covariates,
                       label=FALSE,
                       loadings=FALSE,
                       label.repel=F,
                       size=3,
                       label.size=3,
                       selectedGenes=NULL) {

  pbm <- getPseudoBulkMatrix(name,
                             pbm_dir,
                             matrix,
                             metadata,
                             cellID,
                             covariates)

  matrix <- pbm$pseudo_bulk_matrix
  sampleCovs <- pbm$sampleCovs

  colsToRemove <- intersect(rownames(matrix), colnames(sampleCovs))
  matrix <- matrix[-match(colsToRemove, rownames(matrix)), ]

  require(ggfortify)

  if(length(grep("dataset", colnames(sampleCovs)))>0 & label==F) {
    shape="dataset"
  } else {
    shape=FALSE
  }

  matrix_list <- list()
  matrix_list[["complete"]] <- matrix

  if(!is.null(selectedGenes)) {
    matrix_list[["filtered"]] <- matrix[which(rownames(matrix) %in% selectedGenes), ]
  }

  plots <- list()

  for(i in 1:length(matrix_list)) {
    mymatrix <- matrix_list[[i]]
    mymatrix <- as.data.frame(t(mymatrix))

    pca <- prcomp(mymatrix, scale = TRUE)
    data <- as.data.frame(pca$x)
    data$dataset <- sampleCovs$dataset

    for (cov in covariates) {
      data$cov <- sampleCovs[, match(cov, colnames(sampleCovs))]

      p <- autoplot(pca, data = data, colour = "cov", shape=shape, label=label, loadings=loadings, size=size, label.size=label.size, label.repel=label.repel) +
                    theme_classic() +
                    theme(legend.position="right",
                          legend.key.width = unit(0.6, 'cm'),
                          legend.key.height = unit(0.3, 'cm'),
                          legend.text = element_text(size=10))  +
                    labs(color=gsub("_", " ", cov)) +
                    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                    ggtitle(paste0(nrow(matrix), " genes"))

      if(length(table(data$cov))>3) {
        p <- p + scale_colour_gradientn(colours = c("red", "blue"))
      } else {
        p <- p + scale_colour_manual(values=c("darkgreen", "lightgreen"))
      }

      plots[[paste0(names(matrix_list)[i], "_", cov)]] <- p
    }
  }

  return(plots)
}



#' getSeed it returns a seed for data sampling
#'
#' @param myseeds_list the list of seeds the user provide to reproduce the sampling and thus the models already created
#' @param seeds_list if myseeds_list is not provide, we create a list of seeds, the seeds_list, for the sampling
#' @return a seed
#' @export
#' @examples
#'
getSeed <- function(myseeds_list=NULL,
                    seeds_list) {

  if(is.null(myseeds_list)) {
    seed <- sample(1000:9999, 1)

    while(seed %in% seeds_list) {
      seed <- sample(1000:9999, 1)
    }
  } else {
    seed <- myseeds_list[[times]]
  }

  return(seed)

}



#' getSampling it split the donors in x folds as many times as specified.
#' It is a covariate-dependent sampling, which means we try to put a balanced proportion of each class in each fold when it is a categorical variable.
#' @param N the number of donors to split into x folds x times
#' @param nfolds the number of folds in which donors will be distributed
#' @param ntimes the number of times donors will be distributed into x folds
#' @param cov the covariate to take into account to get balanced classes in each fold
#' @param myseeds_list the user can specify a seed to reproduce one of the sampling and, therefore, one of the models already created
#' @return a list
#' @export
#' @examples
#'
getSampling <- function(N,
                        nfolds=3,
                        ntimes=10,
                        cov,
                        myseeds_list=NULL) {

  # If any covariate provide or the covariate is numerical, donors are randomly distribute in x folds x times
  if (is.null(cov) || is.numeric(cov)) {
    times <- 1
    foldid_list <- list()
    seeds_list <- list()

    while(times<ntimes+1) {
      seed <- getSeed(myseeds_list=myseeds_list, seeds_list=seeds_list)
      set.seed(seed)
      foldid <- sample(1:nfolds,size=N,replace=TRUE)
      all <- table(foldid)

      while(length(all)<nfolds) { # sometimes donors are distributed in only nfolds-1 folds, we repeat the process until we get the number of folds specified
        foldid <- sample(1:nfolds,size=N,replace=TRUE)
        all <- table(foldid)
      }
      foldid_list[[times]] <- foldid
      seeds_list[[times]] <- seed
      times <- times+1
    }
  }

  # If the covariate is categorical, donors are distributed so we the classes are balanced in each fold
  if (is.factor(cov)) {
    N1 <- which(cov==levels(cov)[1])
    N2 <- which(cov==levels(cov)[2])
    times <- 1
    foldid_list <- list()
    seeds_list <- list()

    while(times<ntimes+1) {
      n <- -1

      while(n<0) {
        seed <- getSeed(myseeds_list=myseeds_list, seeds_list=seeds_list)
        set.seed(seed)

        foldid <- sample(1:nfolds,size=N,replace=TRUE)
        all <- table(foldid)
        tab_N1 <- table(foldid[N1])

        if(length(all)==nfolds & length(tab_N1)==nfolds) { # if donors are distributed in x folds and class A is represented in all of the folds

          tab_N1 <- (tab_N1/all)*100
          index_1 <- c(which(tab_N1>70), which(tab_N1<30))

          if(length(index_1)>0) { # if the classes are not balanced for each fold, repeat the sampling
            n <- -1
          } else { # else, save this model
            n <- 0
          }
        } else { # else, we repeat the sampling
          n <- -1
        }
      }
      foldid_list[[times]] <- foldid
      seeds_list[[times]] <- seed
      times <- times + 1
    }
  }

  return(list(foldid_list=foldid_list,
              seeds_list=seeds_list))
}





#' getGLMNETmodels it applies cv.glmnet to find the most relevant genes to predict a covariate (i.e. diagnosis, age)
#'
#' @param data a pseudo-bulk expression matrix containing as many rows as genes + additional covariates
#' @param cov the covariate we want to predict
#' @param parallel if TRUE, the process of creating the models will be parallelized
#' @param ntimes_split the number of times we want to split the discovery dataset into train and test
#' @param ntimes_cv the number of times we want to split the train set into x folds for cross-validation purposes
#' @param nfolds_cv the number of folds we want to use for the cross validation (it depends on the number of donors we have)
#' @param myseeds_train_test the seed we used to split the discovery dataset into train and test to recreate a specific model of interest
#' @param myseeds_cv the seed we used to split the train set into x folds to recreate a specific model of interest
#' @return a list containing:
#' i) all the models
#' ii) a data frame containing the main features of all the models
#' iii) a data frame containing the genes selected as relevant in each model
#' iv) the pseudo-bulk matrix used for creating the models
#' v) the covariate we want to predict
#' vi) the name of the donors we used for training in each model created
#' @export
#' @examples
#'

getGLMNETmodels <- function(data,
                            cov,
                            parallel=TRUE,
                            ntimes_split=3,
                            ntimes_cv=3,
                            nfolds_cv=3,
                            myseeds_train_test=NULL,
                            myseeds_cv=NULL) {

  # Load libraries
  require(caret)
  require(glmnet)

  # Parallelization
  if(parallel==TRUE) {
    require(doParallel)
    registerDoParallel(cores=2)
  }

  # Linear regression or logistic regression?
  tab <- table(cov)

  if(length(tab)==2) {
    cov <- as.factor(as.character(cov))
    family <- "binomial"
    type <- "class"
    type.measure <- "auc"
  } else {
    cov <- as.numeric(cov)
    family <- "gaussian"
    type <- "response"
    type.measure <- "mse"
  }

  # Data split into train and test
  sampling <- getSampling(N=ncol(data),
                          nfolds=3,
                          ntimes=ntimes_split,
                          cov,
                          myseeds_list=myseeds_train_test)

  foldid_train_test <- sampling[["foldid_list"]]
  seeds_train_test <- sampling[["seeds_list"]]
  count <- 0

  # Create empty dataframes or list
  models_summary <- data.frame()
  genes_summary <- data.frame()
  models_list <- list()
  train_list <- list()

  # For each dataset split into train and test
  for (split_tt in 1:length(foldid_train_test)) {

    # 2/3 of the data are selected for training and 1/3 for test
    folds_to_select <- as.numeric(as.character(names(sort(table(foldid_train_test[[split_tt]]), decreasing=T))))[c(2,3)]
    samplesToSelect <- which(foldid_train_test[[split_tt]] %in% folds_to_select)
    cat("Number of samples for training:", length(samplesToSelect), "\n")
    cat("Number of samples for test:", ncol(data)-length(samplesToSelect), "\n")
    data.train <-  t(data[,samplesToSelect])
    data.test <-  t(data[,-samplesToSelect])
    cov.train <-  cov[samplesToSelect]
    cov.test <- cov[-samplesToSelect]
    train_list[[split_tt]] <- colnames(data)[samplesToSelect] # save the ids of the donors used for training

    # Split training set into x folds x times for cross-validation purposes
    sampling2 <- getSampling(N=nrow(data.train),
                             nfolds=nfolds_cv,
                             ntimes=ntimes_cv,
                             cov=cov.train,
                             myseeds_list=myseeds_cv)

    foldid_cv <- sampling2[["foldid_list"]]
    seeds_cv <- sampling2[["seeds_list"]]

    # For each set of folds created for cross-validation purposes
    for(split_cv in 1:length(foldid_cv)) {
      cat("- Train_test iter", split_tt, "CV iter", split_cv, "\n")

      cvfit<- glmnet::cv.glmnet(data.train,
                                cov.train,
                                alpha=1, # LASSO regularization mechanism
                                family = family,
                                keep=TRUE, # keep which donors belong to each fold in the model output
                                type_measure=type.measure,
                                parallel=parallel,
                                foldid=foldid_cv[[split_cv]]) # we use our own train split into x folds
      count <- count + 1

      # Save the coefficients of the model
      coeffs <- as.data.frame(t(as.matrix(stats::coef(cvfit,s="lambda.min"))))

      # Make predictions with both train and test to check the model
      predict.test <- stats::predict(cvfit, newx = data.test,type = type, s = "lambda.min")
      predict.train <- stats::predict(cvfit, newx = data.train,type = type, s = "lambda.min")

      # Save main features of the model
      lambda.min <- cvfit["lambda.min"]
      cvm <- cvfit$cvm[cvfit$index[1]]
      cvsd <- cvfit$cvsd[cvfit$index[1]]
      nzero <- cvfit$nzero[cvfit$index[1]]


      if(type=="response") {
        train_cor <- cor(predict.train, cov.train)
        train_R2 <- MLmetrics::R2_Score(predict.train, cov.train)
        train_rmse <- MLmetrics::RMSE(predict.train, cov.train)

        test_cor <- cor(predict.test, cov.test)
        test_R2 <- MLmetrics::R2_Score(predict.test, cov.test)
        test_rmse <- MLmetrics::RMSE(predict.test, cov.test)
        myRow <- unlist(c(seeds_train_test[[split_tt]], seeds_cv[[split_cv]], nrow(data.train), nrow(data.test), lambda.min, cvm, cvsd, nzero, train_cor, train_R2, train_rmse, test_cor, test_R2, test_rmse))
        names(myRow) <- c("seed_train_test", "seed_cv", "trainSamples", "testSamples", "lambda.min", "cvm", "cvsd", "nzero", "train_cor", "train_R2", "train_rmse", "test_cor", "test_R2", "test_rmse")
      } else {
        cm <- confusionMatrix(as.factor(as.vector(predict.test)),as.factor(cov.test), positive="1")
        accuracy <- as.vector(cm$overall[1])
        kappa <- as.vector(cm$overall[2])
        positive <- cm$positive
        sensitivity <- as.vector(cm$byClass[1])
        specificity <- as.vector(cm$byClass[2])
        ppv <- as.vector(cm$byClass[3])
        npv <- as.vector(cm$byClass[4])
        precision <- as.vector(cm$byClass[5])
        myRow <- c(seeds_train_test[[split_tt]], seeds_cv[[split_cv]], nrow(data.train), nrow(data.test), lambda.min, cvm, cvsd, nzero, accuracy, kappa, positive, sensitivity, specificity, ppv, npv, precision)
        names(myRow) <- c("seed_train_test", "seed_cv", "trainSamples", "testSamples", "lambda.min", "cvm", "cvsd", "nzero", "accuracy", "kappa", "positive", "sensitivity", "specificity", "ppv", "npv", "precision")
      }
      models_list[[count]] <- cvfit # save the model
      models_summary <- rbind(models_summary, myRow) # save the main statistics of this model as a row in a data frame
      colnames(models_summary) <- names(myRow)
      genes_summary <- rbind(genes_summary, cbind(seeds_train_test[[split_tt]], seeds_cv[[split_cv]], coeffs)) # save the relevant predictors selected in this model as a row in a data frame
    }

  }

  return(list(models=models_list,
              models_summary=models_summary,
              genes_summary=genes_summary,
              data=data,
              cov=cov,
              train_list=train_list))
}





#' getGLMNETmodelsPerCov it applies cv.glmnet to find the most relevant genes to predict a covariate (i.e. diagnosis, age)
#'
#' @param name the name of the file that contains the pseudo-bulk matrix if it has been already created
#' @param pbm_dir the directory where the pseudo-bulk matrix has been stored if it has been already created
#' @param matrix a single-cell gene expression matrix with as many rows as genes and as many columns as cells
#' @param metadata a data frame containing the metadata of each cell
#' @param cellID the name of the variable that represents the ID of each cell
#' @param name the name of the cell type corresponding to this expression matrix
#' @param discovery_dataset the name of the dataset the user want to use to create the models
#' @param outliers the id of the donors we want to exclude from the analysis because they are detected as outliers in the PCA
#' @param covariates the name of the covariates for which we want to create a model
#' @param models_dir the name of the directory where the models will be saved
#' @param parallel if TRUE, the process of creating the models will be parallelized
#' @param ntimes_split the number of times we want to split the discovery dataset into train and test
#' @param ntimes_cv the number of times we want to split the train set into x folds for cross-validation purposes
#' @param nfolds_cv the number of folds we want to use for the cross validation (it depends on the number of donors we have)
#' @param myseeds_train_test the seed we used to split the discovery dataset into train and test to recreate a specific model of interest
#' @param myseeds_cv the seed we used to split the train set into x folds to recreate a specific model of interest
#' @return all the models are saved in the models directory indicated by the user
#' @export
#' @examples
#'
getGLMNETmodelsPerCov <- function(name,
                                  pbm_dir,
                                  matrix=NULL,
                                  metadata=NULL,
                                  cellID=NULL,
                                  discovery_dataset=NULL,
                                  outliers=NULL,
                                  covariates=c("Donor_age_years", "Brain_Bank_Path_Dx", "Sex"),
                                  models_dir=getwd(),
                                  parallel=TRUE,
                                  ntimes_split=10,
                                  nfolds_cv=3,
                                  ntimes_cv=10,
                                  myseeds_train_test=NULL,
                                  myseeds_cv=NULL) {


  # Get the pseudo-bulk matrix
  pbm <- getPseudoBulkMatrix(name,
                             pbm_dir,
                             matrix,
                             metadata,
                             cellID,
                             covariates)

  # Select the discovery dataset to create the models
  if(is.null(discovery_dataset)) {
    tab <- table(pbm$sampleCovs$dataset)
    discovery_dataset <- names(tab)[which.max(tab)]
  }

  ids <- pbm$sampleCovs$sampleID[pbm$sampleCovs$dataset==discovery_dataset]

  # Decide if you want to discard outliers
  if(!is.null(outliers)) {
    ids <- setdiff(ids, outliers)
  }
  pbm$pseudo_bulk_matrix <- pbm$pseudo_bulk_matrix[,which(colnames(pbm$pseudo_bulk_matrix) %in% ids)]
  pbm$sampleCovs <- pbm$sampleCovs[which(pbm$sampleCovs$sampleID %in% ids), ]
  stopifnot(identical(colnames(pbm$pseudo_bulk_matrix), pbm$sampleCovs$sampleID))

  # For each covariate
  for (covariate in covariates) {

    # If it is a categorical variable, check if the classes are more or less balanced
    tab <- table(pbm$sampleCovs[, covariate])
    if (length(tab)==2) {
      proportion <- (max(tab)/sum(tab))*100
    } else {# it is a numerical variable
      proportion=0
    }

    if(proportion<60) {

      # Select the covariate
      matrix <- pbm$pseudo_bulk_matrix[-which(rownames(pbm$pseudo_bulk_matrix) %in% covariate), ]
      cov <- as.vector(unlist(pbm$pseudo_bulk_matrix[which(rownames(pbm$pseudo_bulk_matrix) %in% covariate), ]))

      # Create a set of models
      model <- getGLMNETmodels(data=matrix,
                               cov=cov,
                               parallel=parallel,
                               ntimes_split=ntimes_split,
                               nfolds_cv=nfolds_cv,
                               ntimes_cv=ntimes_cv,
                               myseeds_train_test=myseeds_train_test,
                               myseeds_cv=myseeds_cv)

      # Save the model
      cat("- Saving the model for the cell type ",  name, "and covariate", covariate, "\n")
      saveRDS(model, paste0(models_dir, "/", covariate, "_", name, ".rds"))
    } else {
      cat("The classes of the covariate", covariate, "are too unbalanced to create a model. We skip this covariate.\n")
    }
  }

  return()
}





