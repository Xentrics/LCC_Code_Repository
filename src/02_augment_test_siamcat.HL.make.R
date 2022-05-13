#!/usr/bin/env
library("magrittr")
library("tidyverse")
library("drake")

source("src/functions.R")


# Default ggplot2 theme
theme_my <-
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    axis.line.x = ggplot2::element_line(size = 0.8),
    axis.line.y = ggplot2::element_line(size = 0.8),
    axis.ticks = ggplot2::element_line(colour = "black", size = 0.8),
    axis.text = ggplot2::element_text(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

ggplot2::theme_set(theme_my)


covars = c("Gender", "BMI")

tTargets <- c("Candida", "Candida_albicans", "Candida_glabrata", "Candida_sake")



#== H/L Plans ====
get_train_plan <- function(max_expand = NULL) {

  drake::drake_plan(

    #== ABUNDANCE Profiles OF TRAIN SAMPLES (only to subset test features to accommodate variance drift and batch effects) ====
    tFeatures.hu2.strat_full.prev30.train = {

      to_remove <- c("UNMAPPED", "UNINTEGRATED", "UNGROUPED", "OTHER", "LOW_PREV")

      read_rds(file_in("data/profile.strat_full.main.rds")) %>%
        filter(Profiler %in% c("HU2", "Meta2")) %>%
        dplyr::select(-Summed, -counts_kept, -Profiler) %>%
        filter(Category %in% c("Bac", "MetaCyc", "EC", "Kpath", "Kterm")) %>%
        # remove description columns, if present
        mutate(Profile = map(Profile, ~ .x %>% dplyr::select(-any_of("Description")))) %>%
        # 30% prevalence filter
        mutate(Profile = map(Profile, filter_by_prevalance, prev = 0.30)) %>%
        # keep only the feature-vector
        mutate(Features_Train = map(Profile, ~ .x %>% dplyr::select(Feature) %>% filter(!(Feature %in% to_remove)) %>% pull(Feature))) %>%
        dplyr::select(-Profile)
    },
    

    #== ABUNDANCE PROFILES OF VALIDATION SAMPLES ====
    # Feature-Space of train & test was the same before prevalance filtering
    # Features are filtered by the ones kept in train
    tData.hu2.strat_full.test = {

      to_remove <- c("UNMAPPED", "UNINTEGRATED", "UNGROUPED", "OTHER", "LOW_PREV")

      transpose <- function(x) x %>% pivot_longer(-Feature, names_to = "Sample", values_to = "Abundance") %>% pivot_wider(names_from = "Feature", values_from = "Abundance")

      profs <- read_rds(file_in("data/profile.strat_full.test.rds")) %>%
        filter(Profiler %in% c("HU2", "Meta2")) %>%
        dplyr::select(-Summed, -counts_kept, -Profiler) %>%
        filter(Category %in% c("Bac", "MetaCyc", "EC", "Kpath")) %>%
        # remove description columns, if present
        mutate(Profile = map(Profile, ~ .x %>% dplyr::select(-any_of("Description")))) %>%
        # keep only selected samples
        mutate(Profile = map(Profile, ~ .x %>% dplyr::select(Feature, all_of(tSamples.test)))) %>%
        # add features used in train
        left_join(tFeatures.hu2.strat_full.prev30.train) %>%
        # filter features by those used in train
        mutate(Profile = map2(Profile, Features_Train, ~ keep_features(df = .x, to_keep = .y))) %>%
        # NOTE: in contrast to previous results, we remove unwanted features BEFORE normalization
        mutate(Profile = map(Profile, ~ .x %>% filter(!(Feature %in% to_remove)))) %>%
        # TSS normalize
        mutate(Profile = map(Profile, normalize_tss)) %>%
        # transpose features into columns
        mutate(Profile = map(Profile, transpose)) %>%
        # sort samples
        mutate(Profile = map(Profile, arrange, Sample))

      # sanity check: see if all features from train were present in test
      profs <- profs %>%
        mutate(Features_Missing = map2_int(Profile, Features_Train, ~ length(setdiff(.y, colnames(.x)))))

      stopifnot(profs$Features_Missing == 0)

      # target value
      profs %>% dplyr::select(-Features_Train, -Features_Missing)
    },

    # Names of test samples
    tSamples.test = {
      read_rds(file_in("data/meta_data_2021.test.rds")) %>%
        dplyr::rename("Sample" = "PatientID") %>%
        # remove incomplete entries (here: Resected cohort has not age)
        filter(!is.na(Age)) %>%
        # remove COPD & Resected cohort
        filter(Cohort == "Immuno") %>%
        filter(!preTreatment) %>%
        pull(Sample) %>%
        sort()
    },

    # Meta data of test samples
    tMeta.test = {
      tar <- read_rds(file_in("data/targets.test.rds")) %>%
        dplyr::rename("Sample" = "samples") %>%
        rename_with(~ .x %>% str_remove("fungi_") %>% str_remove("[sg]__"), contains("Candida"))
      meta <- read_rds(file_in("data/meta_data_2021.test.rds")) %>%
        dplyr::rename("Sample" = "PatientID")

      # Only keep samples with complete meta-data information
      meta <- tar %>%
        left_join(meta) %>%
        filter(Sample %in% tSamples.test) %>%
        arrange(Sample)

      # sanity check: make sure samples are matched with abundance profiles
      # remark: if we check here, we already know samples from all other profiles are identical
      stopifnot(tData.hu2.strat_full.test$Profile[[1]]$Sample == meta$Sample)

      # target value
      meta
    },


    #== Select Models for testing ====
    tFiles.strat_full.hu2 = {

      perf_meta <- tibble()
      
      # models from local computation
      if (file.exists("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.xls")) {
        perf_meta <- bind_rows(perf_meta, gdata::read.xls(file_in("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.xls")))
      }
      
      if (file.exists("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.top10_MINIMAL.xls")){
        perf_meta <- bind_rows(perf_meta, gdata::read.xls(file_in("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.top10_MINIMAL.xls")))
      }
      
      # models from the manuscript
      perf_meta <- bind_rows(perf_meta, gdata::read.xls(file_in("manuscript_models/performance_test_specific_summary.xls")))

      # select only the best models
      perf_meta <- bind_rows(
        # best performing
        perf_meta %>%
          filter(auROC > 0.77 | auPRC > 0.77) %>%
          select(target, fname) %>%
          filter(file.exists(fname)) %>%
          mutate(hash = map_chr(fname, ~ digest::digest(.x, file = TRUE))),
        # top 10 per target
        perf_meta %>%
          group_by(target) %>%
          slice_head(n = 10) %>%
          select(target, fname) %>%
          filter(file.exists(fname)) %>%
          mutate(hash = map_chr(fname, ~ digest::digest(.x, file = TRUE)))
      ) %>%
        filter(!duplicated(hash))

      perf_meta
    },


    #==== Load and evaluate models of test data ====
    tTest.strat_full.hu2 = target(
      {
        meta.test = tMeta.test %>%
          dplyr::select(Sample, all_of(tTargets), all_of(covars)) %>%
          column_to_rownames("Sample")

        evals <- validate_siam(
          model_dat = read_rds(tFiles.strat_full.hu2$fname),
          meta.test = meta.test,
          profile.test = tData.hu2.strat_full.test,
          group = NULL # use group from model directly
        )

        file.remove("Rplots.pdf", showWarnings = F)
        
        tFiles.strat_full.hu2 %>%
          mutate(evals = list(evals)) %>%
          mutate(
            train_auROC = map_dbl(evals, ~ .x[["train.eval"]][["auroc"]]),
            train_auPRC = map_dbl(evals, ~ .x[["train.eval"]][["auprc"]]),
            test_auROC  = map_dbl(evals, ~ .x[["test.eval"]][["auroc"]]),
            test_auPRC  = map_dbl(evals, ~ .x[["test.eval"]][["auprc"]])
          )
      },
      dynamic = map(tFiles.strat_full.hu2)
    ),


    # gather all test results
    tTest.summary = {
      tTest.strat_full.hu2 %>%
        mutate(SPACE = "HU2", MODE = "Strat-Full") %>%
        mutate(MODE = if_else(str_detect(fname, 'Virus'), 'Phage', MODE)) %>%
        arrange(target, desc(test_auPRC)) %>%
        relocate(fname, hash, .after = last_col())
    },

    # write out performance
    tTest.summary.out = target(
      {
        df <- tTest.summary %>% dplyr::select(-evals)

        df %<>%
          left_join(
            df %>%
              pivot_longer(contains("au"), names_to = "metric") %>%
              group_by(target, hash) %>%
              summarise(auMean = mean(value))
          ) %>%
          relocate(auMean, .after = target) %>%
          # order by average au performance
          arrange(target, desc(auMean)) %>%
          mutate(across(where(is.numeric), round, 3)) %>%
          arrange(target, desc(auMean))

          walk(
            "performance_test_summary.xls",
            ~ WriteXLS::WriteXLS(df, .x, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)
          )
      },
      format = "file"
    ),
  )

}


dplan <- get_train_plan()

cache <- storr::storr_rds(".drake_validate")

drake_config(
  plan = dplan,
  cache = cache,
  keep_going = FALSE,
  lock_envir = FALSE,
  lock_cache = TRUE,
  format = "qs",
  recover = TRUE,
  recoverable = TRUE
)

