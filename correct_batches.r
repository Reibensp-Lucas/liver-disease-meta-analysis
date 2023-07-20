# ##############################################################################
#
## Example of batch effect correction
#
# ##############################################################################

#########################################################
## ATTENTION: You NEED to load the fit_q_plit.r script ##
#########################################################
# install.packages("getPass")
library(here)
# source('/g/scb/zeller/karcher/mi-eocrc/Nics_scripts/fit_q_plit.r')
source(here("fit_q_plit.r"))

#library("tidyverse")
#library("SIAMCAT")
library(SIAMCAT)
 require('devtools')
 require('getPass')
 require('git2r')
# #library('quadprog')
 # # devtools::install_git(
 #    'https://git.embl.de/grp-zeller/batch_effect_correction.git',
 #    credentials = git2r::cred_user_pass("reibensp", getPass::getPass())
 #  )
#library('BEATLE')
library(BEATLE)
library(tidyverse)


#library("doParallel")
#library("ConQuR")

#load('/Users/karcher/tmp/profiles_merged_with_metadata.rimage')
#load('/Users/karcher/16S.data.preparation.RData')
# load('/g/scb/zeller/karcher/mi-eocrc/Updated.Rdata/16S.data.preparation.RData')
load(here("16S.data.preparation.RData"))

check_inputs <- function(inputDataBefore, inputMeta, method){
  inputDataBefore <- inputDataBefore[, colnames(inputDataBefore) %in% rownames(inputMeta)]
  
  if (method %in% c("conqur")) {
    # conqur requires counts
    if (max(inputDataBefore) <= 1) {
      stop(str_c("input data need to be counts for chosen method ", method, "."))
    }
  } 
  
  if (method %in% c("BEATLE")) {
    # BEATLE requires (?) relative abundances.
    if (max(inputDataBefore) >= 1) {
      stop(str_c("input data need to be relative abundances for chosen method ", method, "."))
    }
    inputDataBefore <- as.matrix(inputDataBefore)
  }
  
  if (dim(inputDataBefore)[2] != dim(inputMeta)[1]) {
    #print(str_c())
    #print(dim(inputDataBefore))
    #print(dim(inputMeta))
    stop("Input data and metadata dont have same dimensions. Exiting...")
  }
  
  if (!all(colnames(inputDataBefore) == rownames(inputMeta))) {
    print(str_c("Input data and metadata need to be ordered with respect to each other"))
    print(dim(inputDataBefore))
    print(dim(inputMeta))
    stop("Exiting...")
  }
  return(list(inputDataBefore, inputMeta))
}

correctBatchesBEATLE <- function(inputDataBefore, inputMeta, method = "BEATLE", referenceData = NULL, referenceMeta = NULL){
  
  inputsChecked <- check_inputs(inputDataBefore = inputDataBefore, 
                                inputMeta = inputMeta,
                                method = method)
  inputDataBefore <- inputsChecked[[1]]
  inputMeta <- inputsChecked[[2]]
  
  inputsChecked <- check_inputs(inputDataBefore = referenceData, 
                                inputMeta = referenceMeta,
                                method = method)
  referenceData <- inputsChecked[[1]]
  referenceMeta <- inputsChecked[[2]]
  #return(list(inputDataBefore, inputMeta, referenceData, referenceMeta))
  
  # train Zeller model
  sc.obj.zeller <- siamcat(feat=referenceData, meta=referenceMeta,
                           label='condition', case= 'Case')
  sc.obj.zeller <- filter.features(sc.obj.zeller, filter.method = 'prevalence',
                                   cutoff=0.05)
  sc.obj.zeller <- create.data.split(sc.obj.zeller, num.folds = 5,
                                     num.resample = 5)
  sc.obj.zeller <- normalize.features(sc.obj.zeller, norm.method = 'log.std',
                                      norm.param=list(log.n0=1e-05, sd.min.q=0))
  sc.obj.zeller <- train.model(sc.obj.zeller, method='lasso')
  sc.obj.zeller <- make.predictions(sc.obj.zeller)
  sc.obj.zeller <- evaluate.predictions(sc.obj.zeller)
  #return(sc.obj.zeller)
  
  #model.interpretation.plot(sc.obj.zeller)
  
  # apply to Feng data
  # naive!
  sc.obj.feng <- siamcat(feat=inputDataBefore, meta=inputMeta,
                         label='condition', case= 'Case')
  sc.obj.feng <- make.predictions(sc.obj.zeller, sc.obj.feng)
  sc.obj.feng <- evaluate.predictions(sc.obj.feng)
  # model.evaluation.plot(sc.obj.feng)
  
  # batch effect correction!
  # first, make CV split for useful application of the batch effect correction
  sc.obj.feng <- create.data.split(sc.obj.feng,
                                   num.folds = 5, num.resample = 3)
  ref.data <- get.filt_feat.matrix(sc.obj.zeller)
  # data frame to hold the predictions on batch-effect-corrected data
  bec.pred <- matrix(
    NA, nrow=nrow(inputMeta), ncol=data_split(sc.obj.feng)$num.resample,
    dimnames = list(rownames(inputMeta),
                    paste0('pred_',
                           seq_len(data_split(sc.obj.feng)$num.resample))))
  
  transformed_data = list()
  
  # loop through the CV folds
  for (i in seq_len(data_split(sc.obj.feng)$num.resample)){
    message('Resampling round: ', i)
    for (j in seq_len(data_split(sc.obj.feng)$num.folds)){
      message('CV fold: ', j)
      
      train.data <- inputDataBefore[rownames(ref.data),
                              data_split(sc.obj.feng)$training.folds[[i]][[j]]]
      test.data <- inputDataBefore[rownames(ref.data),
                             data_split(sc.obj.feng)$test.folds[[i]][[j]]]
      meta.train <- inputMeta[colnames(train.data),]
      
      # learn the trafo from the controls in the training fold
      res <- correct.batch(
        target.data = train.data,
        reference.data = ref.data,
        log.n0 = 1e-05,
        ref.ctr.idx = rownames(referenceMeta[referenceMeta$condition=='CTR' | referenceMeta$condition=='Control',]),
        trg.ctr.idx = rownames(meta.train[meta.train$condition=='CTR' | meta.train$condition=='Control',]),
        diagnostic.plot = NULL, 
        verbose=1)
      
      # transform the test data
      new.test.data <- test.data
      for (f in rownames(ref.data)){
        trafo <- res$trafo[[f]]
        if (length(trafo) == 1){
          new.test.data[f,] <- log10(test.data[f,] + 1e-05)
          
        } else {
          new.data <- transf.ab(log10(test.data[f,] + 1e-05),
                                n=trafo$n, plit=trafo$plit)
          new.test.data[f,] <- new.data
        }
      }
      new.test.data <- 10^new.test.data
      new.test.data <- new.test.data - 1e-05
      transformed_data[[length(transformed_data) + 1]] <- new.test.data
      message('trafo done!')
      # make predctions on the test data
      sc.obj.test <- siamcat(feat=new.test.data, verbose=0)
      sc.obj.test <- make.predictions(sc.obj.zeller, sc.obj.test, verbose=0)
      preds <- rowMeans(pred_matrix(sc.obj.test))
      bec.pred[names(preds), i] <- preds
      message('predictions done!')
    }
  }
  
  transformed_data <- map(transformed_data, function(x) x %>% as.data.frame()  %>% rownames_to_column('sampleID') %>% pivot_longer(-sampleID))
  transformed_data <- do.call('rbind', transformed_data)
  transformed_data <- transformed_data %>%
    group_by(sampleID, name) %>%
    summarize(value = median(value)) %>%
    pivot_wider(id_cols = name, names_from = sampleID, values_from = value) %>%
    as.data.frame() %>%
    column_to_rownames('name') %>%
    t()
  
  sc.obj.feng.bec <- sc.obj.feng
  pred_matrix(sc.obj.feng.bec) <- bec.pred
  sc.obj.feng.bec <- evaluate.predictions(sc.obj.feng.bec)
  
  return(list('input_data_before' = inputDataBefore,
              'input_model_before' = sc.obj.feng,
              'input_data_after' = transformed_data, 
              'input_model_after' = sc.obj.feng.bec))
  
}

 # Beatle needs a reference dataset. I'm using a random one - Baxter - for now
 ref.meta <-  meta_16S %>%
   rownames_to_column('sampleID') %>%
   filter(cohort == "Baxter") %>%
   as.data.frame() %>%
   column_to_rownames('sampleID')

 ref.data <- data_16S[, colnames(data_16S) %in% (rownames(ref.meta))]
 ref.data <- ref.data[, match(rownames(ref.meta), colnames(ref.data))]

 input.meta <- meta_16S %>%
   #as.data.frame() %>%
   rownames_to_column('sampleID') %>%
   filter(cohort == "Fudan") %>%
   select(sampleID, condition, cohort) %>%
   as.data.frame() %>%
   column_to_rownames('sampleID')

 #input.data <- data_16S[, colnames(data_16S) %in% rownames(im)]
 input.data <- data_16S
 input.data <- input.data[, match(rownames(input.meta), colnames(input.data))]

# #########
# #########

 # ref.data <- models.lo.all@filt_feat[[1]]
 # ref.meta <- models.lo.all@label[[1]] %>%
 #   as.data.frame() %>%
 #   rename(condition = '.') %>%
 #   mutate(condition = ifelse(condition == "-1", "Control", "CRC")) %>%
 #   mutate(condition = factor(condition, levels = c("Control", "CRC"))) %>%
 #   mutate(dummy = "dummy")
 # 
 # input.data <- testProfile
 # set.seed(1)
 # input.meta <- data.frame(condition = ifelse(runif(dim(testProfile)[2]) > 0.5, "CRC", "Control")) %>%
 #   mutate(sampleID = colnames(testProfile)) %>%
 #   column_to_rownames('sampleID') %>%
 #   mutate(dummy = "dummy")
 # 
 # ref.data <- ref.data[, match(rownames(ref.meta), colnames(ref.data))]
 # input.data <- input.data[, match(rownames(input.meta), colnames(input.data))]
ref.meta <- meta_all %>% 
  filter(dataset_name == "Wang2017") %>%
   as.data.frame() %>%
   mutate(condition = ifelse(condition == 'Control', "Control" , "Case")) %>% 
   column_to_rownames('sampleID')
 
 ref.data <- wide_matrix[, colnames(wide_matrix) %in% (rownames(ref.meta))]
 ref.data <- ref.data[, match(rownames(ref.meta), colnames(ref.data))]
 
 input.meta <- meta_all %>%
   #as.data.frame() %>%
   #rownames_to_column('sampleID') %>%
   filter(dataset_name == "Yun2019") %>%
   select(sampleID, condition, dataset_name) %>%
   mutate(condition = ifelse(condition == 'Control', "Control" , "Case")) %>%
   as.data.frame() %>%
   column_to_rownames('sampleID')
 
 #input.data <- data_16S[, colnames(data_16S) %in% rownames(im)]
 input.data <- wide_matrix
 input.data <- input.data[, match(rownames(input.meta), colnames(input.data))]
 
 
 results <- correctBatchesBEATLE(inputMeta = input.meta, inputDataBefore = input.data, referenceData = ref.data, referenceMeta = ref.meta)

 write.table(results$input_data_after, 
            here('..', 'data', 'batch_corrected_counts', "Yun2019_counts.txt"),
             sep = "\t", row.names = TRUE, col.names = TRUE)
 # write.table(ref.data, here('..', 'data', 'batch_corrected_counts', "Wang2017_counts.txt"),
 #             sep = "\t", row.names = TRUE, col.names = TRUE)
 