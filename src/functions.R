
require("RhpcBLASctl")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

options(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/focal/2021-02-25+Y3JhbiwyOjc4NDM0OTs4NzIzN0YzQg"), download.file.method = 'libcurl')
options(BioC_mirror = 'https://packagemanager.rstudio.com/bioconductor')

#== General Purpose Functions ====

#' Data.frame prevalence filter function
#' It assumes that ANY numeric COLUMN represents a SAMPLE
#' Other columns will be ignored
#' @df Dataframe with samples in columns
#' @prev Percentage value of sample prevalence. Number of samples will be
#'   determined based on the number of samples found within df
filter_by_prevalance <- function(df, prev = 0.20) {

  #' Calculate number of samples with non-zero OTU count
  #' OTUs in rows
  prevalance <- function(df) {
    df %>%
      dplyr::select(where(is.numeric)) %>%
      mutate(prev = rowSums(. > 0)) %>%
      pull(prev)
  }

  nSamples <- df %>% dplyr::select(where(is.numeric)) %>% ncol
  thresh <- nSamples * prev

  tmp <- df %>%
    mutate(prev = prevalance(.), .after = 1) %>%
    # accumulate low prevalence counts
    mutate(prev = prev >= thresh)

  tmp.accum <- tmp %>%
    filter(prev == FALSE) %>%
    group_by(prev) %>%
    summarise(Feature = "LOW_PREV", across(where(is.numeric), sum)) %>%
    dplyr::select(-prev)

  res <- bind_rows(
    tmp %>% filter(prev == TRUE) %>% dplyr::select(-prev),
    tmp.accum
  )

  # sanity check
  select_num <- function(x) x %>% dplyr::select(where(is.numeric))
  stopifnot(all.equal(colSums(select_num(res)), colSums(select_num(df)),tolerance = 1e-10))

  return(res)
}



#' Data.frame filter function
#' @df Dataframe. One colum called 'Features' and samples in columns. All numeric
#' columns are considered abundance columns!
#' @to_keep Character(n). List of features to keep in df. Other features will
#' be summed up to 'LOW_PREV'.
keep_features <- function(df, to_keep) {

  tmp <- df %>% mutate(keep = Feature %in% to_keep, .after = "Feature")

  tmp.accum <- tmp %>%
    filter(!keep) %>%
    group_by(keep) %>%
    summarise(Feature = "LOW_PREV", across(where(is.numeric), sum)) %>%
    dplyr::select(-keep)

  res <- bind_rows(
    tmp %>% filter(keep) %>% dplyr::select(-keep),
    tmp.accum
  )

  # sanity check
  select_num <- function(x) x %>% dplyr::select(where(is.numeric))
  stopifnot(all.equal(colSums(select_num(res)), colSums(select_num(df)),tolerance = 1e-10))

  return(res)
}



roc_from_siam <- function(eval_dat) {

  au_label <- sprintf("auROC: %.3f", eval_dat$auroc)

  roc_ci_df <- tibble(
    sp = as.numeric(rownames(eval_dat$roc$ci)),
    se.low = eval_dat$roc$ci[,1],
    se.median = eval_dat$roc$ci[,2],
    se.high = eval_dat$roc$ci[,3],
  )


  roc_df <- tibble(
    Sensitivity = eval_dat$roc$sensitivities,
    Specificity = eval_dat$roc$specificities
  ) %>%
    mutate(
      TPR = 1 - Specificity
    ) %>%
    arrange(TPR, Sensitivity)


    # confidence interval
    ggplot() +
    geom_ribbon(aes(x = 1-se.median, ymin = se.low, ymax = se.high), data = roc_ci_df, alpha = 0.10) +
    # main ROC line +
    geom_step(aes(x = TPR, y = Sensitivity), data = roc_df, inherit.aes = F) +
    geom_abline(slope = 1, linetype = 3) +
    ggtitle("ROC Curve") +
    labs(x = "False positive rate",
         y = "True positive rate") +
    annotate("text", x = 0.60, y = 0.15, label = au_label) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(fill = NA, colour = "black"),
      axis.text = element_text(color = "black")
    )
}


score_F1 <- function(tp, fn) {
  tp / (tp + 0.5 * (fp + fn))
}


#' Evaluate prediction performance of a trained SIAMCAT model on novel data
#'
#' @param model_dat List. Must contain a character vector 'catos' and the trained
#'   siam model object in 'siam'
#' @param meta.test Meta-data for test cohort
#' @param profile.test Abundance profiles for test cohort
#' @param return_model If TRUE, the train-test SIAMCAT model object will be returned as well
validate_siam <- function(model_dat, meta.test, profile.test, return_model = FALSE, group = "Candida") {

  library("SIAMCAT")
  message(model_dat$fname)
  catos <- model_dat$catos[[1]]

  if (is.null(group) || is.na(group) || group == "") group <- model_dat$target

  prof <- profile.test %>%
    filter(Category %in% catos) %>%
    pull(Profile) %>%
    purrr::reduce(left_join) %>%
    column_to_rownames("Sample") %>%
    t()

  # Create SIAMCAT objects
  sc.train.obj <- model_dat$siam[[1]]
  sc.test.obj <- siamcat(feat = prof, meta = meta.test, label = group, case = "High")

  # Make predictions
  sc.test.obj <- make.predictions(
    siamcat = sc.train.obj,
    siamcat.holdout = sc.test.obj,
    normalize.holdout = TRUE)

  sc.test.obj <- evaluate.predictions(sc.test.obj)

  model.evaluation.plot('TRAIN-CV' = sc.train.obj,
                        'TEST'     = sc.test.obj,
                        colours    = c('dimgrey', 'orange'))

  # target value
  res <- list(train.eval = eval_data(sc.train.obj),
       test.eval  = eval_data(sc.test.obj))

  if (return_model) {
    res <- c(res, "mod.train"  = sc.train.obj, "mod.test"   = sc.test.obj)
  }

  return(res)
}


replace_train_classes <- function(model_dat, new_class_vec) {

  if (!all(sort(names(model_dat$siam[[1]]@label$info)) == sort(levels(new_class_vec))))
    stop("Class labels must be the same!")

  ### convert siam-labels & update
  mod_labels <- model_dat$siam[[1]]@label
  # make sure sample order is correct
  new_class_vec %<>% .[names(mod_labels$label)]
  # replace class labels by 1, -1 indices
  new_class_vec.siam <- mod_labels$info[new_class_vec] %>% set_names(names(new_class_vec))

  # update label object
  model_dat$siam[[1]]@label$label <- new_class_vec.siam
  return(model_dat)
}


normalize_clr <- function(prof) {
  res <- prof %>%
    pivot_longer(-Feature, names_to = "samples", values_to = "Abundance") %>%
    pivot_wider(names_from = "Feature", values_from = "Abundance") %>%
    column_to_rownames("samples") %>%
    # baysian zero-replacement
    zCompositions::cmultRepl() %>%
    compositions::clr() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Feature")
  
  # sanity check: sample means should be close to zero
  stopifnot(colMeans(res[,-1]) < 1e-10)
  return(res)
}


normalize_tss <- function(prof) {
  mat <- prof %>%
    pivot_longer(-Feature, names_to = "samples", values_to = "Abundance") %>%
    pivot_wider(names_from = "Feature", values_from = "Abundance") %>%
    column_to_rownames("samples")
  
  mat <- mat / rowSums(mat)
  
  mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Feature")
}



#== Feature Importance Functions ====

#' @param mod SIAMCAT model
#' @param min_robustness Minimum percentage of models that must contain a variable
siamcat_feature_importance <- function(mod, min_robustness = 0.50) {

  # extract features & weights
  # Re-implementation based on SIAMCAT output
  # this was done to include quantiles
  SIAMCAT::weight_matrix(mod) %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-1, names_to = "model", values_to = "weight") %>%
    # compute relative weights per model
    group_by(model) %>%
    mutate(rel.weight = weight / sum(abs(weight))) %>%
    # compute model weights
    group_by(Feature) %>%
    summarise(
      robustness = round( sum(abs(weight) > 0) / n(), 2),
      mean.weight = mean(weight),
      median.weight = median(weight),
      sd.weight = sd(weight),
      q25.weight = quantile(weight, 0.25),
      q75.weight = quantile(weight, 0.75),

      mean.rel.weight = mean(rel.weight),
      median.rel.weight = median(rel.weight),
      sd.rel.weight = sd(rel.weight),
      q25.rel.weight = quantile(rel.weight, 0.25),
      q75.rel.weight = quantile(rel.weight, 0.75),
    ) %>%
    # compute relative weights
    ungroup() %>%
    filter(robustness >= min_robustness) %>%
    arrange(desc(robustness))


  # SIAMCAT::feature_weights(mod) %>%
  #   rownames_to_column("Feature") %>%
  #   as_tibble() %>%
  #   rename("robustness" = "percentage") %>%
  #   # keep only frequently selected features
  #   filter(robustness >= min_robustness) %>%
  #   arrange(desc(robustness))
}



#' Load SIAMCAT model & compute average feature importance
#' @param fname Path to saved model file
#' @param min_robustness Minimum percentage of cross-validation runs to report a feature
load_siamcat_model_with_feat_imp <- function(fname, min_robustness = 0.25) {

  mod_tbl <- read_rds(fname)
  mod_tbl %>%
    # get feature weights & robustness
    mutate(feat_imp = list(siamcat_feature_importance(mod = siam[[1]], min_robustness = min_robustness)), .after = "eval") %>%
    # label encoding
    mutate(label_info = map(siam, ~ .x@label$info), .after = feat_imp) %>%
    select(-siam) %>%
    # compute accumulative explained weight
    mutate(weight_expl_mean   = map_dbl(feat_imp, ~ sum(abs(.x$mean.rel.weight)))) %>%
    mutate(weight_expl_median = map_dbl(feat_imp, ~ sum(abs(.x$median.rel.weight))))
}


#' Plot average feature importance of SIAMCAT models
#' @param model_dat
#' @param model_id
#' @param label_info
#' @param min_robustness
plot_importance <- function(model_dat, model_id, label_info, min_robustness = 0.5) {

  by <- "mean.rel.weight"
  expl_weight <- sum(abs(model_dat[[by]]))

  # Effect plot by feature weights
  dat <- model_dat %>%
    arrange(!!sym(by)) %>%
    mutate(Feature = fct_inorder(Feature)) %>%
    mutate(Group = factor(sign(!!sym(by)), levels = rev(label_info), labels = names(label_info))) %>%
    filter(robustness >= min_robustness) %>%
    mutate(robustness = sprintf("%.0f%s", 100 * robustness, "%"))

  max_x <- max(dat[[by]], dat$q75.rel.weight)
  min_x <- min(dat[[by]], dat$q25.rel.weight)
  pos_labels <- min_x - (max_x - min_x) * 0.05

  dat %>%
    ggplot(aes(x = !!sym(by), y = Feature, color = Group)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(xmin = q25.rel.weight, xmax = q75.rel.weight)) +
    theme_bw(base_size = 10) +
    scale_x_continuous(expand = expansion(c(0.1, 0.05))) +
    ggtitle(label = model_id, subtitle = sprintf("Explained Weight: %.2f%s", expl_weight, "%")) +
    geom_text(aes(label = robustness, x = pos_labels), size = 3, hjust = "right")
}




#== Data Augmentation ====

#'
#' @param siamcat SIAMCAT object, including a train-test slot!
#' @param nAugSamp Number of augmented samples to create PER GROUP. nAugSamp is further modified by param 'equalize class'.
#' For nAugSamp == 0, the 'siamcat' argument is returned without any modification.
#' @param equalize_class If TRUE, augmentation will result in equal number of samples per class. Thereby, 'nAugSamp' will
#'   be INCREASED for the SMALLER class (and LOWERED for the HIGHER class), i.e.
#'   nAugSamp * prop(class1) + nAugSamp * prop(class2) = nAugSamp
#' @example
#'   augment_siamcat(siamcat, 10, FALSE)
augment_siamcat = function(siamcat, nAugSamp, equalize_class = TRUE) {

  if (class(siamcat) != "siamcat") stop("siamcat is NOT a siamcat class object!")
  if (!("data_split" %in% names(attributes(siamcat)))) stop("SIAMCAT object must contain a data_split definition! Use SIAMCAT::create.data.split!")
  if (siamcat@label$type != "BINARY") stop("Augmentation does possible: it was only implemented using Binary class. Please augment function first!")

  if (nAugSamp == 0) return(siamcat)
  if (nAugSamp < 0) stop("Number of augmented samples must not be negative!")

  # Total-Sum scale abundance
  tss <- function(x) t(t(x) / colSums(x))

  siamcat@phyloseq@otu_table@.Data <- tss(siamcat@phyloseq@otu_table@.Data)

  # deep copy for later modification
  sc.obj.aug <- rlang::duplicate(siamcat)


  otu_table <- otu_table(sc.obj.aug@phyloseq)
  split_info <- sc.obj.aug@data_split$training.folds

  nRepeats <- length(split_info)
  for(rep_i in 1:length(split_info)) {
    for(k in 1:length(split_info[[rep_i]])) {

      # We perform augmentation within each training fold
      nFolds <- length(split_info[[rep_i]])
      train_samples <- split_info[[rep_i]][[k]]
      message(sprintf("Augment: Rep [%i / %i], fold [%i / %i]", rep_i, nRepeats, k, nFolds))


      # use label data is used to divide training samples by their classes
      class_df <- data.frame(class = sc.obj.aug@label$label) %>%
        rownames_to_column("Sample") %>%
        group_by(class) %>%
        nest() %>%
        ungroup() %>%
        rename("samples" = "data") %>%
        mutate(samples = map(samples, unlist, use.names = F)) %>%
        # remove samples that are not IN THIS FOLD
        mutate(samples = map(samples, ~ intersect(.x, train_samples))) %>%
        # separate samples by each class
        mutate(prof = map(samples, ~ otu_table[, .x, drop=F]))

      # Adapt nSamp to remove class imbalance
      class_df %<>% mutate(nAugSamp = nAugSamp, .after = class)
      if (equalize_class) {
        class_df %<>%
          mutate(classProp = map_dbl(samples, ~ 1 - length(.x) / length(train_samples))) %>%
          mutate(nAugSamp = round(nrow(.) * nAugSamp * classProp))
      }


      # create augmented samples PER CLASS
      class_augmented_df <- class_df %>%
        # for equalize_class == TRUE, nAugSamp can be zero for some classes
        filter(nAugSamp > 0) %>%
        mutate(augmented = map2(prof, nAugSamp, ~ augment_profile(otus = t(.x@.Data), nSamp = .y, add = F))) %>%
        # update names of augmented samples (make them unique)
        mutate(augmented = map2(augmented, class, ~ set_rownames(.x, rownames(.x) %>% str_replace("AUG", sprintf("AUG_%s_rep%.2i_fold%.2i_", .y, rep_i, k)) ))) %>%
        # TSS normalize
        mutate(augmented.tss = map(augmented, ~ .x / rowSums(.x))) %>%
        # Convert to phyloseq OTU table
        mutate(augmented.phy = map(augmented.tss, phyloseq::otu_table, taxa_are_rows = F))


      # update SIAMCAT object with new augmented training data & update sample information too
      # merge profiles from ALL classes together
      new_profiles <- t(reduce(class_augmented_df %>% pull(augmented.phy), rbind))
      # get names of new augmented samples
      class_augmented_df %<>% mutate(samples.augmented = map(augmented, rownames))
      new_samples <- unlist(class_augmented_df$samples.augmented)
      # create new class labels for augmented samples
      class_augmented_df %<>% mutate(labels.augmented = map2(class, samples.augmented, function(cl, samp) rep(cl, length(samp)) %>% set_names(samp)))
      new_labels <- unlist(class_augmented_df$labels.augmented)

      # update SIAMCAT object
      sc.obj.aug@phyloseq@otu_table@.Data                %<>% cbind(new_profiles)
      sc.obj.aug@data_split$training.folds[[rep_i]][[k]] %<>% c(new_samples) %>% unique()
      sc.obj.aug@label$label %<>% c(new_labels)

      # sanity checks
      stopifnot(new_samples %in% sc.obj.aug@data_split$training.folds[[rep_i]][[k]])
      stopifnot(colnames(new_profiles) %in% colnames(sc.obj.aug@phyloseq@otu_table))
      stopifnot(length(new_labels) == ncol(new_profiles))
    }
  }


  return(sc.obj.aug)
}


#' @param otus Samples must be in rows, features in columns
#' @param nSamp Number of augmented samples to return
#' @param add If TRUE, 'otus' will be returned alognside the new augmented samples
augment_profile <- function(otus, nSamp, add = FALSE) {

  exp <- min(3 + floor(log10(ncol(otus))), 7) # 10 million reads should be maximum (to revert RPKb)
  comm <- round( otus * 10^exp )

  set.seed(235)
  X <- SpiecEasi::synth_comm_from_counts(comm = comm, mar = 2, distr='negbin', n = nSamp)
  colnames(X) <- colnames(comm)
  rownames(X) <- sprintf("AUG%.3i", 1:nrow(X))

  if (add) {
    # merge old & new profiles
    return(rbind(comm, X))
  } else {
    return(X)
  }
}


#' @param siamcat SIAMCAT object with augmented samples
#' @param aug_pref Sample prefix used to indicate augmentation samples
#' @param only_labels if TRUE, only augmented labels will be removed. The profiles remain. This one is required to
#' successfully run SIAMCAT::make.predictions.
remove_augmented <- function(siamcat, aug_pref = "AUG_", only_labels = TRUE) {

  siamcat@label$label %<>% .[!str_detect(names(.), "AUG_")]

  if (!only_labels) {
    # remove augmented profiles as well
    siamcat@phyloseq@otu_table@.Data %<>% .[ , !str_detect(colnames(.), "AUG_") ]
    siamcat@norm_feat$norm.feat %<>% .[ , !str_detect(colnames(.), "AUG_") ]
  }
  return(siamcat)
}

