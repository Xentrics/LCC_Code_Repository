#!/usr/bin/env

library("tidyverse")
library("drake")
library("magrittr")
library("SIAMCAT")

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



#== H/L Plans ====
get_train_plan_full_strat <- function(max_expand = NULL) {

  # to allow for more parallel execution,
  # each covar contains its own drake cache
  covars = c("Gender", "BMI")

  # generate targets
  rows_as_list <- function(x) lapply(1:nrow(x), function(i) x[i,])
  get_unified_names <- function(l) sapply(l, function(x) paste(sort(x), collapse = "_"))


  # Limit the number of categories/Features in order to not overload vita
  cc <- c("Bac", "EC", "MetaCyc", "Kpath")

  # column names of class to classify
  tTargets <- c("Candida", "Candida_albicans", "Candida_glabrata", "Candida_sake")

  drake::drake_plan(

    tData.strat_full.prev30 = {

      to_remove <- c("UNMAPPED", "UNINTEGRATED", "UNGROUPED", "OTHER", "LOW_PREV")

      transpose <- function(x) x %>% pivot_longer(-Feature, names_to = "Sample", values_to = "Abundance") %>% pivot_wider(names_from = "Feature", values_from = "Abundance")

      profs <- read_rds(file_in("data/profile.strat_full.main.rds")) %>%
        filter(Profiler %in% c("HU2", "Meta2", "ARG")) %>%
        dplyr::select(-Summed, -counts_kept, -Profiler) %>%
        filter(Category %in% c("Bac", "MetaCyc", "EC", "Kpath")) %>%
        # remove description columns, if present
        mutate(Profile = map(Profile, ~ .x %>% dplyr::select(-any_of("Description")))) %>%
        # 20% prevalence filter
        mutate(Profile = map(Profile, filter_by_prevalance, prev = 0.30)) %>%
        # NOTE: in contrast to previous results, we remove unwanted features BEFORE normalization
        mutate(Profile = map(Profile, ~ .x %>% filter(!(Feature %in% to_remove)))) %>%
        # TSS normalize
        mutate(Profile = map(Profile, normalize_tss)) %>%
        # transpose features into columns
        mutate(Profile = map(Profile, transpose)) %>%
        # sort samples
        mutate(Profile = map(Profile, arrange, Sample))

      # sanity check: all samples must be the same
      samples <- profs$Profile[[1]]$Sample
      mm <- profs %>%
        mutate(sample_mismatch = map_dbl(Profile, ~ sum(.x$Sample != samples))) %>%
        pull(sample_mismatch)
      stopifnot(mm == 0)

      # check passed!
      # target value
      profs
    },


    tMeta = {
      tar <- read_rds(file_in("data/targets.main.rds")) %>%
        dplyr::rename("Sample" = "samples") %>%
        rename_with(~ .x %>% str_remove("fungi_") %>% str_remove("[sg]__"), contains("Candida"))
      meta <- read_rds(file_in("data/meta_data_2021.main_imputed.rds")) %>%
        dplyr::rename("Sample" = "PatientID")

      # sanity check
      meta <- tar %>%
        left_join(meta) %>%
        arrange(Sample)

      # sanity check: make sure samples are matched with abundance profiles
      # remark: if we check here, we already know samples from all other profiles are identical
      stopifnot(tData.strat_full.prev30$Profile[[1]]$Sample == meta$Sample)

      # target value
      meta
    },




    #== Model Training purely using SIAMCAT functionality ====

    tData.strat_full.augment = target(
      {
        meta <- tMeta %>%
          dplyr::select(Sample, Gender, BMI, all_of(group)) %>%
          column_to_rownames("Sample")

        # Filter data & convert to matrix
        prof <- tData.strat_full.prev30 %>%
          filter(Category %in% catos) %>%
          pull(Profile) %>%
          purrr::reduce(left_join) %>%
          column_to_rownames("Sample") %>%
          t()

        # Create SIAMCAT object
        sc.obj <- siamcat(feat = prof, meta = meta, label = group, case = "High")

        # Create data-split slots
        set.seed(5345)
        sc.obj <-  create.data.split(
          siamcat = sc.obj,
          num.folds = 10,
          num.resample = 1
        )


        # Augment Data
        sc.obj <- tryCatch(
          {
            sc.obj <- augment_siamcat(siamcat = sc.obj, nAugSamp = nAug, equalize_class = TRUE)

            # Normalize data
            l0 <- if (is.null(log.n0)) {NULL} else {min(prof[prof!=0])/10}
            sc.obj <- normalize.features(
              siamcat = sc.obj,
              norm.method = norm_method,
              feature.type = "original",
              norm.param = list(
                log.n0 = l0,
                n.p = 2,
                norm.margin = 1
              )
            )
            sc.obj
          },
          error = function(e) {
            NA
          }
        )


        # target value
        sc.obj
      },
      transform = cross(
        group        = !!tTargets,
        catos        = !!cc,
        nAug         = c(0, 10, 20, 30, 40),
        norm_method  = c("log.unit", "rank.unit", "log.clr"),
        log.n0       = list("divby10"),
      ),

      trigger = trigger(command = FALSE, depend = FALSE),
    ),



    #==== Models using Bacteria-only full-stratified profiles ====

    models_simple.wi_feat.full_strat.prev30.ENET = target(
      {
        basedir  <- sprintf("results_species_prev30_augmented/models/HU2_full_strat_prev30_%s_%s", group, model_method)
        basename <- sprintf("%s/%s.%s_a%.2f.%s", basedir, norm_method, model_method, alpha, paste0(catos, collapse = "_"))
        if (perform.fs) basename <- sprintf("%s.FS_%s_%i", basename, method.fs, thres.fs)
        basename <- sprintf("%s.aug_%i", basename, nAug)
        dir.create(basedir, recursive = T, showWarnings = F)

        message(sprintf("Writing model to: [%s]", basename))

        # SIAMCAT object already contains normalized profiles, inluding CV split information & normalization
        sc.obj <- tData.strat_full.augment

        set.seed(5234)
        sc.obj <- tryCatch(
          {
            train.model(
              siamcat = sc.obj,
              method = model_method,
              param.set = list("alpha" = alpha),
              stratify = TRUE,
              modsel.crit = 'f1',
              perform.fs = perform.fs,
              param.fs = list(thres.fs = thres.fs, method.fs = method.fs, direction = "absolute")
            )
          }, error = function(e) {
            res <- tibble(
              eval         = list(auROC = 0, auPR = 0),
              target       = group,
              catos        = list(catos),
              norm_method  = norm_method,
              log.n0       = log.n0,
              model_method = model_method,
              nAug         = nAug,
              perform.fs   = perform.fs,
              method.fs    = method.fs,
              thres.fs     = thres.fs,
              fname        = NA
            )

            return(res)
          }
        )

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed

        # # access model
        # models <- models(sc.obj)
        # models[[1]]

        # make predictions within each CV
        # Note: Since we never actually predict augmented samples, they will not
        # get a prediction value. This will lead to a crash in make.predictions.
        # Likewise, we do not want the augmented samples to have direct influence
        # on our AUC. Hence, we kick them out here.

        sc.obj.noAug <- remove_augmented(siamcat = sc.obj, only_labels = TRUE) # remove only augmented labels, but keep profiles

        sc.obj <- tryCatch({
          make.predictions(sc.obj.noAug, verbose = 1)
        }, error = function(e) {
          res <- tibble(
            eval         = list(auROC = 0, auPR = 0),
            target       = group,
            catos        = list(catos),
            norm_method  = norm_method,
            log.n0       = log.n0,
            model_method = model_method,
            nAug         = nAug,
            perform.fs   = perform.fs,
            method.fs    = method.fs,
            thres.fs     = thres.fs,
            fname        = NA
          )

          return(res)
        })

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed
        # pred_matrix <- pred_matrix(sc.obj)

        sc.obj <-  evaluate.predictions(sc.obj)
        eval_dat <- eval_data(sc.obj)

        # evaluate using ROC and PR
        tryCatch(
          {
            #gg_roc <- roc_from_siam(eval_dat)
            #model.evaluation.plot(sc.obj, fn.plot = sprintf("%s.eval_plot.pdf", basename))

            model.interpretation.plot(sc.obj, sprintf("%s.interpret_plot.pdf", basename))
          }, error = function(e) {
            warning("Plotting failed.")
          }
        )


        # path to save model information
        fname <- sprintf("%s.rds", basename)

        # Merge Information
        res <- tibble(
          siam   = list(sc.obj),
          eval    = list(eval_dat),
          target       = group,
          #gg_roc = list(gg_roc),
          catos        = list(catos),
          norm_method  = norm_method,
          log.n0       = log.n0,
          model_method = model_method,
          nAug         = nAug,
          perform.fs   = perform.fs,
          method.fs    = method.fs,
          thres.fs     = thres.fs,
          fname        = fname
        )

        # save model to drive
        write_rds(x = res, file = fname)

        # remove model object from drake target to save stuff (that job is out-sourced to the rds file)
        res$siam <- NULL

        # target value
        res
      },

      transform = cross(
        tData.strat_full.augment,
        model_method = c("enet"),
        alpha        = c(0.1, 0.5, 1.0),
        perform.fs   = c(TRUE),
        method.fs    = c("AUC"),
        thres.fs     = c(40, 60, 80, 100)
      ),

      trigger = trigger(command = FALSE, depend = FALSE),
    ),


    models_simple.wo_feat.full_strat.prev30.ENET = target(
      {
        basedir  <- sprintf("results_species_prev30_augmented/models/HU2_full_strat_prev30_%s_%s", group, model_method)
        basename <- sprintf("%s/%s.%s_a%.2f.%s", basedir, norm_method, model_method, alpha, paste0(catos, collapse = "_"))
        if (perform.fs) basename <- sprintf("%s.FS_%s_%i", basename, method.fs, thres.fs)
        basename <- sprintf("%s.aug_%i", basename, nAug)
        dir.create(basedir, recursive = T, showWarnings = F)

        message(sprintf("Writing model to: [%s]", basename))

        # SIAMCAT object already contains normalized profiles, inluding CV split information & normalization
        sc.obj <- tData.strat_full.augment

        set.seed(5234)
        sc.obj <- tryCatch(
          {
            train.model(
              sc.obj,
              method = model_method,
              param.set = list("alpha" = alpha),
              stratify = TRUE,
              modsel.crit = 'f1',
              perform.fs = perform.fs,
              param.fs = list(thres.fs = thres.fs, method.fs = method.fs, direction = "absolute")
            )
          }, error = function(e) {
            res <- tibble(
              eval         = list(auROC = 0, auPR = 0),
              target       = group,
              catos        = list(catos),
              norm_method  = norm_method,
              log.n0       = log.n0,
              model_method = model_method,
              nAug         = nAug,
              perform.fs   = perform.fs,
              method.fs    = method.fs,
              thres.fs     = thres.fs,
              fname        = NA
            )

            return(res)
          }
        )

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed

        # # access model
        # models <- models(sc.obj)
        # models[[1]]

        # make predictions within each CV
        # Note: Since we never actually predict augmented samples, they will not
        # get a prediction value. This will lead to a crash in make.predictions.
        # Likewise, we do not want the augmented samples to have direct influence
        # on our AUC. Hence, we kick them out here.

        sc.obj.noAug <- remove_augmented(siamcat = sc.obj, only_labels = TRUE) # remove only augmented labels, but keep profiles

        sc.obj <- tryCatch({
          make.predictions(sc.obj.noAug, verbose = 1)
        }, error = function(e) {
          res <- tibble(
            eval         = list(auROC = 0, auPR = 0),
            target       = group,
            catos        = list(catos),
            norm_method  = norm_method,
            log.n0       = log.n0,
            model_method = model_method,
            nAug         = nAug,
            perform.fs   = perform.fs,
            method.fs    = method.fs,
            thres.fs     = thres.fs,
            fname        = NA
          )

          return(res)
        })

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed
        # pred_matrix <- pred_matrix(sc.obj)

        sc.obj <-  evaluate.predictions(sc.obj)
        eval_dat <- eval_data(sc.obj)

        # evaluate using ROC and PR
        tryCatch(
          {
            #gg_roc <- roc_from_siam(eval_dat)
            #model.evaluation.plot(sc.obj, fn.plot = sprintf("%s.eval_plot.pdf", basename))

            model.interpretation.plot(sc.obj, sprintf("%s.interpret_plot.pdf", basename))
          }, error = function(e) {
            warning("Plotting failed.")
          }
        )


        # path to save model information
        fname <- sprintf("%s.rds", basename)

        # Merge Information
        res <- tibble(
          siam   = list(sc.obj),
          eval   = list(eval_dat),
          target       = group,
          #gg_roc = list(gg_roc),
          catos  = list(catos),
          norm_method = norm_method,
          log.n0      = log.n0,
          model_method = model_method,
          nAug         = nAug,
          perform.fs   = perform.fs,
          method.fs    = method.fs,
          thres.fs     = thres.fs,
          fname        = fname
        )

        # save model to drive
        write_rds(x = res, file = fname)

        # remove model object from drake target to save stuff (that job is out-sourced to the rds file)
        res$siam <- NULL

        # target value
        res
      },

      transform = cross(
        tData.strat_full.augment,
        model_method = c("enet"),
        alpha        = c(0.1, 0.5, 1.0),
        perform.fs   = c(FALSE),
        method.fs    = c(NA_character_),
        thres.fs     = c(NA_integer_)
      ),

      trigger = trigger(command = FALSE, depend = FALSE),
    ),


    models_simple.wi_feat.full_strat.prev30 = target(
      {
        basedir  <- sprintf("results_species_prev30_augmented/models/HU2_full_strat_prev30_%s_%s", group, model_method)
        basename <- sprintf("%s/%s.%s.%s", basedir, norm_method, model_method, paste0(catos, collapse = "_"))
        if (perform.fs) basename <- sprintf("%s.FS_%s_%i", basename, method.fs, thres.fs)
        basename <- sprintf("%s.aug_%i", basename, nAug)
        dir.create(basedir, recursive = T, showWarnings = F)

        message(sprintf("Writing model to: [%s]", basename))

        # SIAMCAT object already contains normalized profiles, inluding CV split information & normalization
        sc.obj <- tData.strat_full.augment

        set.seed(5234)
        sc.obj <- tryCatch(
          {
            train.model(
              siamcat = sc.obj,
              method = model_method,
              stratify = TRUE,
              modsel.crit = 'f1',
              perform.fs = perform.fs,
              param.fs = list(thres.fs = thres.fs, method.fs = method.fs, direction = "absolute")
            )
          }, error = function(e) {
            res <- tibble(
              eval         = list(auROC = 0, auPR = 0),
              target       = group,
              catos        = list(catos),
              norm_method  = norm_method,
              log.n0       = log.n0,
              model_method = model_method,
              nAug         = nAug,
              perform.fs   = perform.fs,
              method.fs    = method.fs,
              thres.fs     = thres.fs,
              fname        = NA
            )

            return(res)
          }
        )

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed

        # # access model
        # models <- models(sc.obj)
        # models[[1]]

        # make predictions within each CV
        # Note: Since we never actually predict augmented samples, they will not
        # get a prediction value. This will lead to a crash in make.predictions.
        # Likewise, we do not want the augmented samples to have direct influence
        # on our AUC. Hence, we kick them out here.

        sc.obj.noAug <- remove_augmented(siamcat = sc.obj, only_labels = TRUE) # remove only augmented labels, but keep profiles

        sc.obj <- tryCatch({
          make.predictions(sc.obj.noAug, verbose = 1)
        }, error = function(e) {
          res <- tibble(
            eval         = list(auROC = 0, auPR = 0),
            target       = group,
            catos        = list(catos),
            norm_method  = norm_method,
            log.n0       = log.n0,
            model_method = model_method,
            nAug         = nAug,
            perform.fs   = perform.fs,
            method.fs    = method.fs,
            thres.fs     = thres.fs,
            fname        = NA
          )

          return(res)
        })

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed
        # pred_matrix <- pred_matrix(sc.obj)

        sc.obj <-  evaluate.predictions(sc.obj)
        eval_dat <- eval_data(sc.obj)

        # evaluate using ROC and PR
        tryCatch(
          {
            #gg_roc <- roc_from_siam(eval_dat)
            #model.evaluation.plot(sc.obj, fn.plot = sprintf("%s.eval_plot.pdf", basename))

            model.interpretation.plot(sc.obj, sprintf("%s.interpret_plot.pdf", basename))
          }, error = function(e) {
            warning("Plotting failed.")
          }
        )


        # path to save model information
        fname <- sprintf("%s.rds", basename)

        # Merge Information
        res <- tibble(
          siam   = list(sc.obj),
          eval    = list(eval_dat),
          target       = group,
          #gg_roc = list(gg_roc),
          catos        = list(catos),
          norm_method  = norm_method,
          log.n0       = log.n0,
          model_method = model_method,
          nAug         = nAug,
          perform.fs   = perform.fs,
          method.fs    = method.fs,
          thres.fs     = thres.fs,
          fname        = fname
        )

        # save model to drive
        write_rds(x = res, file = fname)

        # remove model object from drake target to save stuff (that job is out-sourced to the rds file)
        res$siam <- NULL

        # target value
        res
      },

      transform = cross(
        tData.strat_full.augment,
        model_method = c("lasso_ll"),
        perform.fs   = c(TRUE),
        method.fs    = c("AUC"),
        thres.fs     = c(40, 60, 80, 100)
      ),

      trigger = trigger(command = FALSE, depend = FALSE),
    ),


    models_simple.wo_feat.full_strat.prev30 = target(
      {
        basedir  <- sprintf("results_species_prev30_augmented/models/HU2_full_strat_prev30_%s_%s", group, model_method)
        basename <- sprintf("%s/%s.%s.%s", basedir, norm_method, model_method, paste0(catos, collapse = "_"))
        if (perform.fs) basename <- sprintf("%s.FS_%s_%i", basename, method.fs, thres.fs)
        basename <- sprintf("%s.aug_%i", basename, nAug)
        dir.create(basedir, recursive = T, showWarnings = F)

        message(sprintf("Writing model to: [%s]", basename))

        # SIAMCAT object already contains normalized profiles, inluding CV split information & normalization
        sc.obj <- tData.strat_full.augment

        set.seed(5234)
        sc.obj <- tryCatch(
          {
            train.model(
              sc.obj,
              method = model_method,
              stratify = TRUE,
              modsel.crit = 'f1',
              perform.fs = perform.fs,
              param.fs = list(thres.fs = thres.fs, method.fs = method.fs, direction = "absolute")
            )
          }, error = function(e) {
            res <- tibble(
              eval         = list(auROC = 0, auPR = 0),
              target       = group,
              catos        = list(catos),
              norm_method  = norm_method,
              log.n0       = log.n0,
              model_method = model_method,
              nAug         = nAug,
              perform.fs   = perform.fs,
              method.fs    = method.fs,
              thres.fs     = thres.fs,
              fname        = NA
            )

            return(res)
          }
        )

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed

        # # access model
        # models <- models(sc.obj)
        # models[[1]]

        # make predictions within each CV
        # Note: Since we never actually predict augmented samples, they will not
        # get a prediction value. This will lead to a crash in make.predictions.
        # Likewise, we do not want the augmented samples to have direct influence
        # on our AUC. Hence, we kick them out here.

        sc.obj.noAug <- remove_augmented(siamcat = sc.obj, only_labels = TRUE) # remove only augmented labels, but keep profiles

        sc.obj <- tryCatch({
          make.predictions(sc.obj.noAug, verbose = 1)
        }, error = function(e) {
          res <- tibble(
            eval         = list(auROC = 0, auPR = 0),
            target       = group,
            catos        = list(catos),
            norm_method  = norm_method,
            log.n0       = log.n0,
            model_method = model_method,
            nAug         = nAug,
            perform.fs   = perform.fs,
            method.fs    = method.fs,
            thres.fs     = thres.fs,
            fname        = NA
          )

          return(res)
        })

        if (!isS4(sc.obj)) return(sc.obj) # if TRUE, model training has failed
        # pred_matrix <- pred_matrix(sc.obj)

        sc.obj <-  evaluate.predictions(sc.obj)
        eval_dat <- eval_data(sc.obj)

        # evaluate using ROC and PR
        tryCatch(
          {
            #gg_roc <- roc_from_siam(eval_dat)
            #model.evaluation.plot(sc.obj, fn.plot = sprintf("%s.eval_plot.pdf", basename))

            model.interpretation.plot(sc.obj, sprintf("%s.interpret_plot.pdf", basename))
          }, error = function(e) {
            warning("Plotting failed.")
          }
        )


        # path to save model information
        fname <- sprintf("%s.rds", basename)

        # Merge Information
        res <- tibble(
          siam   = list(sc.obj),
          eval   = list(eval_dat),
          target       = group,
          #gg_roc = list(gg_roc),
          catos  = list(catos),
          norm_method = norm_method,
          log.n0      = log.n0,
          model_method = model_method,
          nAug         = nAug,
          perform.fs   = perform.fs,
          method.fs    = method.fs,
          thres.fs     = thres.fs,
          fname        = fname
        )

        # save model to drive
        write_rds(x = res, file = fname)

        # remove model object from drake target to save stuff (that job is out-sourced to the rds file)
        res$siam <- NULL

        # target value
        res
      },

      transform = cross(
        tData.strat_full.augment,
        model_method = c("lasso_ll"),
        perform.fs   = c(FALSE),
        method.fs    = c(NA_character_),
        thres.fs     = c(NA_integer_)
      ),

      trigger = trigger(command = FALSE, depend = FALSE),
    ),


    model_perf.full_strat = target(
      {
        bind_rows(
          models_simple.wo_feat.full_strat.prev30.ENET,
          models_simple.wi_feat.full_strat.prev30.ENET,
          models_simple.wo_feat.full_strat.prev30,
          models_simple.wi_feat.full_strat.prev30
        ) %>%
          select(-any_of("siam")) %>%
          # unlist remaining columns
          ungroup() %>%
          mutate(auROC = map_dbl(eval, ~ ifelse('auroc' %in% names(.x), .x$auroc, NA)), .before = fname) %>%
          mutate(auPRC = map_dbl(eval, ~ ifelse('auprc' %in% names(.x), .x$auprc, NA)), .after = auROC) %>%
          mutate(catos = map_chr(catos, ~ paste(sort(.x), collapse = ", "))) %>%
          arrange(desc(auROC)) %>%
          dplyr::select(-any_of(c("siam", "eval")))
      },
      transform = combine(models_simple.wo_feat.full_strat.prev30.ENET,
                          models_simple.wi_feat.full_strat.prev30.ENET,
                          models_simple.wo_feat.full_strat.prev30,
                          models_simple.wi_feat.full_strat.prev30)
    ),



    model_perf.full_strat.out = target(
      {
        model_perf.formatted <- model_perf.full_strat %>% mutate(across(where(is.numeric), ~ round(.x, digits = 2)))

        walk("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.xls",
             ~ WriteXLS::WriteXLS(model_perf.formatted, ExcelFileName = .x, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1))
      },
      format = "file"
    ),



    model_perf.full_strat.top10.out = target(
      {
        model_perf.formatted <- model_perf.full_strat %>%
          mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>%
          group_by(target) %>%
          slice_head(n = 5)

        walk("results_species_prev30_augmented/performance_AUGMENTED_HU2_full_strat_prev30.top10.xls",
             ~ WriteXLS::WriteXLS(model_perf.formatted, ExcelFileName = .x, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1))
      },
      format = "file"
    ),

    #== END ====
  )
}


dplan <- get_train_plan_full_strat()


dir.create("results_species_prev30_augmented")


options(clustermq.scheduler = "multicore")
library("clustermq")

cache <- storr::storr_rds(".drake_hu2_species")

drake_config(
  plan = dplan,
  cache = cache,
  keep_going = FALSE,
  lock_envir = FALSE,
  lock_cache = TRUE,
  parallelism = "clustermq",
  jobs = getOption("Ncpus", parallel::detectCores(logical = F) / 2),
  verbose = 2,
  format = "qs",
  recover = TRUE,
  recoverable = TRUE
)

