# Public code: supply in-memory objects using R/inputs.R; file I/O and directory settings removed.
# Audit changes and unresolved scientific issues: docs/CODE_AUDIT.md.
ML.Dev.Prog.Sig_modified <- function(train_data, list_train_vali_Data, candidate_genes = NULL, unicox.filter.for.candi = NULL, 
    unicox_p_cutoff = NULL, mode = NULL, single_ml = NULL, alpha_for_Enet = NULL, direction_for_stepcox = NULL, 
    double_ml1 = NULL, double_ml2 = NULL, nodesize = NULL, seed = NULL) {
    if (is.null(alpha_for_Enet) == T) {
        alpha_for_Enet <- 0.1
    }
    else {
        alpha_for_Enet <- alpha_for_Enet
    }
    if (is.null(direction_for_stepcox) == T) {
        direction_for_stepcox <- "both"
    }
    else {
        direction_for_stepcox <- direction_for_stepcox
    }
    if (T) {
        library(Matrix)
        library(survival)
        library(randomForestSRC)
        library(glmnet)
        library(plsRcox)
        library(superpc)
        library(gbm)
        library(CoxBoost)
        library(survivalsvm)
        library(dplyr)
        library(tibble)
        library(BART)
        library(miscTools)
        library(compareC)
        library(ggplot2)
        library(ggsci)
        library(tidyr)
        library(ggbreak)
        library(mixOmics)
        library(data.table)
        invisible(NULL)
        options(stringsAsFactors = FALSE)
    }
    if (T) {
        SigUnicox <- function(gene_list, inputSet, unicox_pcutoff) {
            print("Starting the data preprocess")
            print("Rejecting a null value")
            library(survival)
            inputSet[is.na(inputSet)] = 0
            inputSet <- inputSet %>% as.data.frame()
            inputSet$OS.time <- as.numeric(inputSet$OS.time)
            inputSet <- inputSet[inputSet$OS.time > 0, ]
            gene_list <- gsub("-", ".", gene_list)
            gene_list <- gsub("_", ".", gene_list)
            colnames(inputSet)[4:ncol(inputSet)] <- gsub("-", ".", colnames(inputSet)[4:ncol(inputSet)])
            colnames(inputSet)[4:ncol(inputSet)] <- gsub("_", ".", colnames(inputSet)[4:ncol(inputSet)])
            print("Gets the intersection of genelist and expression profile")
            comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)], gene_list)
            print("Processing the  input representation matrix")
            inputSet <- inputSet[, c("ID", "OS.time", "OS", comsa1)]
            inputSet[, c(1:2)] <- apply(inputSet[, c(1:2)], 2, as.factor)
            inputSet[, c(2:ncol(inputSet))] <- apply(inputSet[, c(2:ncol(inputSet))], 2, as.numeric)
            inputSet <- as.data.frame(inputSet)
            print("Data preprocessing completed")
            display.progress = function(index, totalN, breakN = 20) {
                if (index%%ceiling(totalN/breakN) == 0) {
                  cat(paste(round(index * 100/totalN), "% ", sep = ""))
                }
            }
            print("Stating the univariable cox regression")
            unicox <- data.frame()
            for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
                display.progress(index = i, totalN = ncol(inputSet[, 4:ncol(inputSet)]))
                gene <- colnames(inputSet[, 4:ncol(inputSet)])[i]
                tmp <- data.frame(expr = as.numeric(inputSet[, 4:ncol(inputSet)][, i]), futime = inputSet$OS.time, 
                  fustat = inputSet$OS, stringsAsFactors = F)
                cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
                coxSummary <- summary(cox)
                unicox <- rbind.data.frame(unicox, data.frame(gene = gene, HR = as.numeric(coxSummary$coefficients[, 
                  "exp(coef)"])[1], z = as.numeric(coxSummary$coefficients[, "z"])[1], pvalue = as.numeric(coxSummary$coefficients[, 
                  "Pr(>|z|)"])[1], lower = as.numeric(coxSummary$conf.int[, 3][1]), upper = as.numeric(coxSummary$conf.int[, 
                  4][1]), stringsAsFactors = F), stringsAsFactors = F)
            }
            print("Finished the univariable cox regression")
            selgene <- unicox[which(unicox$pvalue < unicox_pcutoff), "gene"]
            return(selgene)
        }
    }
    rf_nodesize <- nodesize
    seed <- seed
    list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
        colnames(x) = gsub("-", ".", colnames(x))
        return(x)
    })
    candidate_genes = gsub("-", ".", candidate_genes)
    colnames(train_data) = gsub("-", ".", colnames(train_data))
    common_feature = c("ID", "OS.time", "OS", candidate_genes)
    for (i in names(list_train_vali_Data)) {
        common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
    }
    message(paste0("---the number of the raw candidate genes is ", length(candidate_genes), " ---"))
    message(paste0("---the number of the common feature is ", length(common_feature) - 3, " ---"))
    returnIDtoRS = function(rs.table.list, rawtableID) {
        for (i in names(rs.table.list)) {
            rs.table.list[[i]]$ID = rawtableID[[i]]$ID
            rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select("ID", everything())
        }
        return(rs.table.list)
    }
    if (!is.na(rf_nodesize) & !is.na(seed) & mode %in% c("all", "single", "double") & identical(c("ID", "OS.time", 
        "OS"), colnames(train_data)[1:3]) & length(candidate_genes) > 0 & identical(c("ID", "OS.time", "OS"), common_feature[1:3]) & 
        length(common_feature) > 3) {
        message("--- Data preprocessing ---")
        common_feature = c("ID", "OS.time", "OS", candidate_genes)
        for (i in names(list_train_vali_Data)) {
            common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
        }
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x[, common_feature]
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, as.numeric)
            return(x)
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
            return(x)
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x[, c(2:3)] <- apply(x[, c(2:3)], 2, as.numeric)
            return(x)
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x <- x[!is.na(x$OS.time) & !is.na(x$OS), ]
            return(x)
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x <- x[x$OS.time > 0, ]
            return(x)
        })
        list_train_vali_Data <- lapply(list_train_vali_Data, function(x) {
            x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(x) {
                x[is.na(x)] = mean(x, na.rm = T)
                return(x)
            })
            return(x)
        })
        train_data = train_data[, common_feature]
        train_data[, -c(1:3)] <- apply(train_data[, -c(1:3)], 2, as.numeric)
        train_data[, c(1:2)] <- apply(train_data[, c(1:2)], 2, as.factor)
        train_data[, c(2:3)] <- apply(train_data[, c(2:3)], 2, as.numeric)
        if (is.null(unicox_p_cutoff)) {
            unicox_p_cutoff = 0.05
        }
        else {
            unicox_p_cutoff = unicox_p_cutoff
        }
        if (is.null(unicox.filter.for.candi)) {
            unicox.filter.for.candi = T
        }
        else {
            unicox.filter.for.candi = unicox.filter.for.candi
        }
        if (unicox.filter.for.candi) {
            cd.gene = common_feature[-c(1:3)]
            cd.gene1 = SigUnicox(gene_list = cd.gene, inputSet = train_data, unicox_pcutoff = unicox_p_cutoff)
            message(paste0("---the number of the final unicox filtered candidate genes is ", length(cd.gene1), 
                " ---"))
            print(cd.gene1)
            cd.gene2 = c(common_feature[1:3], cd.gene1)
            common_feature = cd.gene2
            train_data = as.data.frame(train_data)[, common_feature]
        }
        else {
            message(paste0("---the number of the final not unicox filtered candidate genes is ", length(common_feature) - 
                3, " ---"))
        }
        est_dd <- as.data.frame(train_data)[, common_feature[-1]]
        val_dd_list <- lapply(list_train_vali_Data, function(x) {
            x[, common_feature[-1]]
        })
        pre_var = common_feature[-c(1:3)]
        if (mode == "all") {
            result <- data.frame()
            ml.res = list()
            riskscore = list()
            message("---1-1 RSF ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
            })
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- "RSF"
            result <- rbind(result, cc)
            ml.res[["RSF"]] = fit
            riskscore[["RSF"]] = rs
            message("---1-2.RSF + CoxBoost ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                  trace = TRUE, start.penalty = 500, parallel = T)
                cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), maxstepno = 500, 
                  K = 10, type = "verweij", penalty = pen$penalty)
                fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                  penalty = pen$penalty)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                    2], type = "lp")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "CoxBoost")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "CoxBoost")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "CoxBoost")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-3.RSF + Enet ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                for (alpha in seq(0.1, 0.9, 0.1)) {
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                      s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("RSF + ", "Enet", "[α=", alpha, "]")
                  result <- rbind(result, cc)
                  ml.res[[paste0("RSF + ", "Enet", "[α=", alpha, "]")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("RSF + ", "Enet", "[α=", alpha, "]")]] = rs
                }
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-4.RSF + GBM ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                best <- which.min(fit$cv.error)
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "GBM")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "GBM")]] = list(fit = fit, best = best)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "GBM")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message("---1-5.RSF + Lasso ---")
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "Lasso")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "Lasso")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "Lasso")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-6.RSF + plsRcox ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                  nt = 10, verbose = FALSE)
                fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "plsRcox")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "plsRcox")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "plsRcox")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-7.RSF + Ridge ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "Ridge")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "Ridge")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "Ridge")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-8.RSF + StepCox ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                for (direction in c("both", "backward", "forward")) {
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("RSF + ", "StepCox", "[", direction, "]")
                  result <- rbind(result, cc)
                  ml.res[[paste0("RSF + ", "StepCox", "[", direction, "]")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("RSF + ", "StepCox", "[", direction, "]")]] = rs
                }
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-9.RSF + SuperPC ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                  2)])
                set.seed(seed)
                fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                repeat {
                  tryCatch({
                    cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                      max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                    break
                  }, error = function(e) {
                    cat("Error:", conditionMessage(e), "\n")
                    cat("Retrying...\n")
                    invisible(NULL)
                  })
                }
                rs <- lapply(val_dd_list2, function(w) {
                  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                    2)])
                  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                    ])], n.components = 1)
                  rr <- as.numeric(ff$v.pred)
                  rr2 <- cbind(w[, 1:2], RS = rr)
                  return(rr2)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "SuperPC")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "SuperPC")]] = list(fit, cv.fit)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "SuperPC")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---1-10.RSF + survival-SVM ---")
            set.seed(seed)
            fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                importance = T, proximity = T, forest = T, seed = seed)
            imp_scores <- fit$importance
            rid <- names(imp_scores)[imp_scores > 0]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("RSF + ", "survival-SVM")
                result <- rbind(result, cc)
                ml.res[[paste0("RSF + ", "survival-SVM")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("RSF + ", "survival-SVM")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
            }
            message("---2.Enet ---")
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            for (alpha in seq(0.1, 0.9, 0.1)) {
                set.seed(seed)
                message(paste0("--- 2.Enet", "[α=", alpha, "] ---"))
                fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Enet", "[α=", alpha, "]")
                result <- rbind(result, cc)
                ml.res[[paste0("Enet", "[α=", alpha, "]")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Enet", "[α=", alpha, "]")]] = rs
            }
            message("---3.StepCox ---")
            for (direction in c("both", "backward", "forward")) {
                message(paste0("---3.StepCox", "[", direction, "]---"))
                fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("StepCox", "[", direction, "]")
                result <- rbind(result, cc)
                ml.res[[paste0("StepCox", "[", direction, "]")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("StepCox", "[", direction, "]")]] = rs
            }
            direction = "both"
            if (T) {
                fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
                rid <- names(coef(fit))
                if (length(rid) > 1) {
                  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                    x[, c("OS.time", "OS", rid)]
                  })
                  set.seed(seed)
                  message(paste0("---3.StepCox", "[", direction, "]", " + CoxBoost ---"))
                  pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                    2)]), trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                    maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                      2], type = "lp")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + CoxBoost")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = rs
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  for (alpha in seq(0.1, 0.9, 0.1)) {
                    set.seed(seed)
                    message(paste0("--- 3.StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "] ---"))
                    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                        s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = rs
                  }
                  set.seed(seed)
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + GBM ---"))
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                  best <- which.min(fit$cv.error)
                  set.seed(seed)
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + GBM")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + GBM")]] = list(fit = fit, best = best)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + GBM")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Lasso ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Lasso")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + plsRcox ---"))
                  set.seed(seed)
                  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                    nt = 10, verbose = FALSE)
                  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + plsRcox")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Ridge ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Ridge")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + RSF ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + RSF")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + RSF")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + RSF")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + SuperPC ---"))
                  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                    featurenames = colnames(est_dd2)[-c(1, 2)])
                  set.seed(seed)
                  fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                  repeat {
                    tryCatch({
                      cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                        max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                      break
                    }, error = function(e) {
                      cat("Error:", conditionMessage(e), "\n")
                      cat("Retrying...\n")
                      invisible(NULL)
                    })
                  }
                  rs <- lapply(val_dd_list2, function(w) {
                    test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                      2)])
                    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                      ])], n.components = 1)
                    rr <- as.numeric(ff$v.pred)
                    rr2 <- cbind(w[, 1:2], RS = rr)
                    return(rr2)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + SuperPC")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = rs
                }
                else {
                  warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                }
            }
            direction = "backward"
            if (T) {
                fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
                rid <- names(coef(fit))
                if (length(rid) > 1) {
                  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                    x[, c("OS.time", "OS", rid)]
                  })
                  set.seed(seed)
                  message(paste0("---3.StepCox", "[", direction, "]", " + CoxBoost ---"))
                  pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                    2)]), trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                    maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                      2], type = "lp")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + CoxBoost")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = rs
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  for (alpha in seq(0.1, 0.9, 0.1)) {
                    set.seed(seed)
                    message(paste0("--- 3.StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "] ---"))
                    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                        s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = rs
                  }
                  set.seed(seed)
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + GBM ---"))
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                  best <- which.min(fit$cv.error)
                  set.seed(seed)
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + GBM")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + GBM")]] = list(fit = fit, best = best)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + GBM")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Lasso ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Lasso")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + plsRcox ---"))
                  set.seed(seed)
                  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                    nt = 10, verbose = FALSE)
                  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + plsRcox")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Ridge ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Ridge")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + RSF ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + RSF")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + RSF")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + RSF")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + SuperPC ---"))
                  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                    featurenames = colnames(est_dd2)[-c(1, 2)])
                  set.seed(seed)
                  fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                  repeat {
                    tryCatch({
                      cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                        max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                      break
                    }, error = function(e) {
                      cat("Error:", conditionMessage(e), "\n")
                      cat("Retrying...\n")
                      invisible(NULL)
                    })
                  }
                  rs <- lapply(val_dd_list2, function(w) {
                    test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                      2)])
                    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                      ])], n.components = 1)
                    rr <- as.numeric(ff$v.pred)
                    rr2 <- cbind(w[, 1:2], RS = rr)
                    return(rr2)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + SuperPC")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + survival-SVM ---"))
                  fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + survival-SVM")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + survival-SVM")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + survival-SVM")]] = rs
                }
                else {
                  warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                }
            }
            direction = "forward"
            if (T) {
                fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
                rid <- names(coef(fit))
                if (length(rid) > 1) {
                  est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                  val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                    x[, c("OS.time", "OS", rid)]
                  })
                  set.seed(seed)
                  message(paste0("---3.StepCox", "[", direction, "]", " + CoxBoost ---"))
                  pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                    2)]), trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                    maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                      2], type = "lp")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + CoxBoost")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + CoxBoost")]] = rs
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  for (alpha in seq(0.1, 0.9, 0.1)) {
                    set.seed(seed)
                    message(paste0("--- 3.StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "] ---"))
                    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                        s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction, "]", " + Enet", "[α=", alpha, "]")]] = rs
                  }
                  set.seed(seed)
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + GBM ---"))
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                  best <- which.min(fit$cv.error)
                  set.seed(seed)
                  fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                    interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + GBM")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + GBM")]] = list(fit = fit, best = best)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + GBM")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Lasso ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Lasso")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Lasso")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + plsRcox ---"))
                  set.seed(seed)
                  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                    nt = 10, verbose = FALSE)
                  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + plsRcox")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + plsRcox")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + Ridge ---"))
                  x1 <- as.matrix(est_dd2[, rid])
                  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                      2)]), s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + Ridge")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + Ridge")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + RSF ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + RSF")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + RSF")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + RSF")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + SuperPC ---"))
                  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                    featurenames = colnames(est_dd2)[-c(1, 2)])
                  set.seed(seed)
                  fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                  repeat {
                    tryCatch({
                      cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                        max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                      break
                    }, error = function(e) {
                      cat("Error:", conditionMessage(e), "\n")
                      cat("Retrying...\n")
                      invisible(NULL)
                    })
                  }
                  rs <- lapply(val_dd_list2, function(w) {
                    test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                      2)])
                    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                      ])], n.components = 1)
                    rr <- as.numeric(ff$v.pred)
                    rr2 <- cbind(w[, 1:2], RS = rr)
                    return(rr2)
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + SuperPC")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + SuperPC")]] = rs
                  message(paste0("--- 3.StepCox", "[", direction, "]", " + survival-SVM ---"))
                  fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("StepCox", "[", direction, "]", " + survival-SVM")
                  result <- rbind(result, cc)
                  ml.res[[paste0("StepCox", "[", direction, "]", " + survival-SVM")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("StepCox", "[", direction, "]", " + survival-SVM")]] = rs
                }
                else {
                  warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                }
            }
            message("---4.CoxBoost ---")
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                  2], type = "lp")))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("CoxBoost")
            result <- rbind(result, cc)
            ml.res[[paste0("CoxBoost")]] = fit
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("CoxBoost")]] = rs
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                for (alpha in seq(0.1, 0.9, 0.1)) {
                  set.seed(seed)
                  message(paste0("--- 4.CoxBoost", " + Enet", "[α=", alpha, "] ---"))
                  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                      s = fit$lambda.min)))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("CoxBoost", " + Enet", "[α=", alpha, "]")
                  result <- rbind(result, cc)
                  ml.res[[paste0("CoxBoost", " + Enet", "[α=", alpha, "]")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("CoxBoost", " + Enet", "[α=", alpha, "]")]] = rs
                }
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message(paste0("--- 4.CoxBoost + ", "GBM ---"))
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                best <- which.min(fit$cv.error)
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "GBM")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "GBM")]] = list(fit = fit, best = best)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "GBM")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message(paste0("--- 4.CoxBoost + ", "Lasso ---"))
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "Lasso")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "Lasso")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "Lasso")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message(paste0("--- 4.CoxBoost + ", "plsRcox ---"))
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                  nt = 10, verbose = FALSE)
                fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "plsRcox")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "plsRcox")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "plsRcox")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message(paste0("--- 4.CoxBoost + ", "Ridge ---"))
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                x1 <- as.matrix(est_dd2[, rid])
                x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "Ridge")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "Ridge")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "Ridge")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                for (direction in c("both", "backward", "forward")) {
                  message(paste0("--- 4.CoxBoost + ", "StepCox", "[", direction, "] ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("CoxBoost + ", "StepCox", "[", direction, "]")
                  result <- rbind(result, cc)
                  ml.res[[paste0("CoxBoost + ", "StepCox", "[", direction, "]")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("CoxBoost + ", "StepCox", "[", direction, "]")]] = rs
                }
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            set.seed(seed)
            message(paste0("--- 4.CoxBoost + ", "SuperPC ---"))
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                  2)])
                set.seed(seed)
                fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                repeat {
                  tryCatch({
                    cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                      max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                    break
                  }, error = function(e) {
                    cat("Error:", conditionMessage(e), "\n")
                    cat("Retrying...\n")
                    invisible(NULL)
                  })
                }
                rs <- lapply(val_dd_list2, function(w) {
                  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                    2)])
                  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                    ])], n.components = 1)
                  rr <- as.numeric(ff$v.pred)
                  rr2 <- cbind(w[, 1:2], RS = rr)
                  return(rr2)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "SuperPC")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "SuperPC")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            message(paste0("--- 4.CoxBoost + ", "survival-SVM ---"))
            set.seed(seed)
            pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), trace = TRUE, 
                start.penalty = 500, parallel = T)
            cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                K = 10, type = "verweij", penalty = pen$penalty)
            fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                penalty = pen$penalty)
            rid <- as.data.frame(coef(fit))
            rid$id <- rownames(rid)
            rid <- rid[which(rid$`coef(fit)` != 0), "id"]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost + ", "survival-SVM")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost + ", "survival-SVM")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost + ", "survival-SVM")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
            }
            message("---5.plsRcox ---")
            set.seed(seed)
            cv.plsRcox.res = cv.plsRcox(list(x = est_dd[, pre_var], time = est_dd$OS.time, status = est_dd$OS), 
                nt = 10, verbose = FALSE)
            fit <- plsRcox(est_dd[, pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("plsRcox")
            result <- rbind(result, cc)
            ml.res[[paste0("plsRcox")]] = fit
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("plsRcox")]] = rs
            message("---6.superpc ---")
            data <- list(x = t(est_dd[, -c(1, 2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 
                2)])
            set.seed(seed)
            fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
            repeat {
                tryCatch({
                  cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                    max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                  break
                }, error = function(e) {
                  cat("Error:", conditionMessage(e), "\n")
                  cat("Retrying...\n")
                  invisible(NULL)
                })
            }
            rs <- lapply(val_dd_list, function(w) {
                test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                  2)])
                ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                  ])], n.components = 1)
                rr <- as.numeric(ff$v.pred)
                rr2 <- cbind(w[, 1:2], RS = rr)
                return(rr2)
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("SuperPC")
            result <- rbind(result, cc)
            ml.res[[paste0("SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("SuperPC")]] = rs
            message("---7.GBM ---")
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph", n.trees = 10000, 
                interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
            best <- which.min(fit$cv.error)
            set.seed(seed)
            fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph", n.trees = best, 
                interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("GBM")
            result <- rbind(result, cc)
            ml.res[[paste0("GBM")]] = list(fit = fit, best = best)
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("GBM")]] = rs
            message("---8.survivalsvm ---")
            fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd, gamma.mu = 1)
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("survival - SVM")
            result <- rbind(result, cc)
            ml.res[[paste0("survival - SVM")]] = fit
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("survival - SVM")]] = rs
            message("---9.Ridge ---")
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
            cv.fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", type.measure = "class")
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                  s = cv.fit$lambda.min)))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("Ridge")
            result <- rbind(result, cc)
            ml.res[[paste0("Ridge")]] = list(fit = fit, cv.fit = cv.fit)
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("Ridge")]] = rs
            message("---10.Lasso ---")
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            rs <- lapply(val_dd_list, function(x) {
                cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                  s = fit$lambda.min)))
            })
            cc <- data.frame(Cindex = sapply(rs, function(x) {
                as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
            })) %>% rownames_to_column("ID")
            cc$Model <- paste0("Lasso")
            result <- rbind(result, cc)
            ml.res[[paste0("Lasso")]] = fit
            rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
            riskscore[[paste0("Lasso")]] = rs
            message(paste0("--- 10.Lasso + ", "CoxBoost ---"))
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                  trace = TRUE, start.penalty = 500, parallel = T)
                cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), maxstepno = 500, 
                  K = 10, type = "verweij", penalty = pen$penalty)
                fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                  penalty = pen$penalty)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                    2], type = "lp")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso + ", "CoxBoost")
                result <- rbind(result, cc)
                ml.res[[paste0("Lasso + ", "CoxBoost")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso + ", "CoxBoost")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            message(paste0("--- 10.Lasso + ", "GBM ---"))
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                best <- which.min(fit$cv.error)
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso + ", "GBM")
                result <- rbind(result, cc)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso + ", "GBM")]] = rs
                ml.res[[paste0("Lasso + ", "GBM")]] = list(fit = fit, best = best)
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            message(paste0("--- 10.Lasso + ", "plsRcox ---"))
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                  nt = 10, verbose = FALSE)
                fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso + ", "plsRcox")
                result <- rbind(result, cc)
                ml.res[[paste0("Lasso + ", "plsRcox")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso + ", "plsRcox")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            message(paste0("--- 10.Lasso + ", "RSF ---"))
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            rid <- rid[-1]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                set.seed(seed)
                fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                  importance = T, proximity = T, forest = T, seed = seed)
                rs <- lapply(val_dd_list2, function(x) {
                  cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso", " + RSF")
                result <- rbind(result, cc)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso + ", "RSF")]] = rs
                ml.res[[paste0("Lasso", " + RSF")]] = fit
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                for (direction in c("both", "backward", "forward")) {
                  message(paste0("--- 10.Lasso + ", "StepCox", "[", direction, "] ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction)
                  rs <- lapply(val_dd_list2, function(x) {
                    cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                  })
                  cc <- data.frame(Cindex = sapply(rs, function(x) {
                    as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                  })) %>% rownames_to_column("ID")
                  cc$Model <- paste0("Lasso + ", "StepCox", "[", direction, "]")
                  result <- rbind(result, cc)
                  ml.res[[paste0("Lasso + ", "StepCox", "[", direction, "]")]] = fit
                  rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                  riskscore[[paste0("Lasso + ", "StepCox", "[", direction, "]")]] = rs
                }
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            message(paste0("--- 10.Lasso + ", "superPC ---"))
            x1 <- as.matrix(est_dd[, pre_var])
            x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
            set.seed(seed)
            fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
            fit$lambda.min
            myCoefs <- coef(fit, s = "lambda.min")
            rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
            if (length(rid) > 1) {
                est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                  x[, c("OS.time", "OS", rid)]
                })
                data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, featurenames = colnames(est_dd2)[-c(1, 
                  2)])
                set.seed(seed)
                fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                repeat {
                  tryCatch({
                    cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                      max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                    break
                  }, error = function(e) {
                    cat("Error:", conditionMessage(e), "\n")
                    cat("Retrying...\n")
                    invisible(NULL)
                  })
                }
                rs <- lapply(val_dd_list2, function(w) {
                  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                    2)])
                  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                    ])], n.components = 1)
                  rr <- as.numeric(ff$v.pred)
                  rr2 <- cbind(w[, 1:2], RS = rr)
                  return(rr2)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso + ", "SuperPC")
                result <- rbind(result, cc)
                ml.res[[paste0("Lasso + ", "SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso + ", "SuperPC")]] = rs
            }
            else {
                warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
            }
            result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
            return(result.cindex.fit)
        }
        else if (mode == "single") {
            result <- data.frame()
            ml.res = list()
            riskscore = list()
            if (single_ml %in% "RSF") {
                message("--- Start RSF ---")
                set.seed(seed)
                fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                  importance = T, proximity = T, forest = T, seed = seed)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- "RSF"
                result <- rbind(result, cc)
                ml.res[["RSF"]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[["RSF"]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "Enet") {
                message("--- Start Enet ---")
                x1 <- as.matrix(est_dd[, pre_var])
                x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha_for_Enet, nfolds = 10)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Enet", "[α=", alpha_for_Enet, "]")
                result <- rbind(result, cc)
                ml.res[[paste0("Enet", "[α=", alpha_for_Enet, "]")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Enet", "[α=", alpha_for_Enet, "]")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "StepCox" & direction_for_stepcox %in% c("both", "backward", "forward")) {
                message("--- Start StepCox ---")
                set.seed(seed)
                fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]")
                result <- rbind(result, cc)
                ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "CoxBoost") {
                message("--- Start CoxBoost ---")
                set.seed(seed)
                pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                  trace = TRUE, start.penalty = 500, parallel = T)
                cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                  K = 10, type = "verweij", penalty = pen$penalty)
                fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                  penalty = pen$penalty)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                    2], type = "lp")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("CoxBoost")
                result <- rbind(result, cc)
                ml.res[[paste0("CoxBoost")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("CoxBoost")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "plsRcox") {
                message("--- Start plsRcox ---")
                set.seed(seed)
                cv.plsRcox.res = cv.plsRcox(list(x = est_dd[, pre_var], time = est_dd$OS.time, status = est_dd$OS), 
                  nt = 10, verbose = T)
                fit <- plsRcox(est_dd[, pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("plsRcox")
                result <- rbind(result, cc)
                ml.res[[paste0("plsRcox")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("plsRcox")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "superpc") {
                message("--- Start superpc ---")
                set.seed(seed)
                data <- list(x = t(est_dd[, -c(1, 2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 
                  2)])
                set.seed(seed)
                fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                repeat {
                  tryCatch({
                    cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                      max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                    break
                  }, error = function(e) {
                    cat("Error:", conditionMessage(e), "\n")
                    cat("Retrying...\n")
                    invisible(NULL)
                  })
                }
                rs <- lapply(val_dd_list, function(w) {
                  test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                    2)])
                  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                    ])], n.components = 1)
                  rr <- as.numeric(ff$v.pred)
                  rr2 <- cbind(w[, 1:2], RS = rr)
                  return(rr2)
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("SuperPC")
                result <- rbind(result, cc)
                ml.res[[paste0("SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("SuperPC")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "GBM") {
                message("--- Start GBM ---")
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph", n.trees = 10000, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                best <- which.min(fit$cv.error)
                set.seed(seed)
                fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd, distribution = "coxph", n.trees = best, 
                  interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("GBM")
                result <- rbind(result, cc)
                ml.res[[paste0("GBM")]] = list(fit = fit, best = best)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("GBM")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "survivalsvm") {
                message("--- Start survivalsvm ---")
                set.seed(seed)
                fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd, gamma.mu = 1)
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("survival - SVM")
                result <- rbind(result, cc)
                ml.res[[paste0("survival - SVM")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("survival - SVM")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "Ridge") {
                message("--- Start Ridge ---")
                x1 <- as.matrix(est_dd[, pre_var])
                x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                set.seed(seed)
                fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
                cv.fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", type.measure = "class")
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = cv.fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Ridge")
                result <- rbind(result, cc)
                ml.res[[paste0("Ridge")]] = list(fit = fit, cv.fit = cv.fit)
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Ridge")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            else if (single_ml %in% "Lasso") {
                message("--- Start Lasso ---")
                x1 <- as.matrix(est_dd[, pre_var])
                x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                set.seed(seed)
                fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                rs <- lapply(val_dd_list, function(x) {
                  cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 2)]), 
                    s = fit$lambda.min)))
                })
                cc <- data.frame(Cindex = sapply(rs, function(x) {
                  as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                })) %>% rownames_to_column("ID")
                cc$Model <- paste0("Lasso")
                result <- rbind(result, cc)
                ml.res[[paste0("Lasso")]] = fit
                rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                riskscore[[paste0("Lasso")]] = rs
                result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                return(result.cindex.fit)
            }
            message("--- Finish and return the fit ---")
        }
        else if (mode == "double" & double_ml1 %in% c("RSF", "StepCox", "CoxBoost", "Lasso") & double_ml2 %in% 
            c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", "survivalsvm", "Ridge", "Lasso")) {
            result <- data.frame()
            ml.res = list()
            riskscore = list()
            if (double_ml1 == "RSF" & double_ml2 %in% c("Enet", "StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", 
                "survivalsvm", "Ridge", "Lasso")) {
                if (double_ml2 == "Enet" & alpha_for_Enet %in% seq(0.1, 0.9, 0.1)) {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ", alpha_for_Enet, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    for (alpha in seq(0.1, 0.9, 0.1)) {
                      set.seed(seed)
                      fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
                      rs <- lapply(val_dd_list2, function(x) {
                        cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 
                          2)]), s = fit$lambda.min)))
                      })
                      cc <- data.frame(Cindex = sapply(rs, function(x) {
                        as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                      })) %>% rownames_to_column("ID")
                      cc$Model <- paste0("RSF + ", "Enet", "[α=", alpha_for_Enet, "]")
                      result <- rbind(result, cc)
                      ml.res[[paste0("RSF + ", "Enet", "[α=", alpha_for_Enet, "]")]] = fit
                      rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                      riskscore[[paste0("RSF + ", "Enet", "[α=", alpha_for_Enet, "]")]] = rs
                      result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                      return(result.cindex.fit)
                    }
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "StepCox" & direction_for_stepcox %in% c("both", "backward", "forward")) {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ", direction_for_stepcox, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction_for_stepcox)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "StepCox", "[", direction_for_stepcox, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "StepCox", "[", direction_for_stepcox, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "StepCox", "[", direction_for_stepcox, "]")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "CoxBoost") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                      2)]), trace = TRUE, start.penalty = 500, parallel = T)
                    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                      penalty = pen$penalty)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                        2], type = "lp")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "CoxBoost")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "CoxBoost")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "CoxBoost")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "plsRcox") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                      nt = 10, verbose = FALSE)
                    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "plsRcox")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "plsRcox")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "plsRcox")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "superpc") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                      featurenames = colnames(est_dd2)[-c(1, 2)])
                    set.seed(seed)
                    fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                    repeat {
                      tryCatch({
                        cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                          max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                        break
                      }, error = function(e) {
                        cat("Error:", conditionMessage(e), "\n")
                        cat("Retrying...\n")
                        invisible(NULL)
                      })
                    }
                    rs <- lapply(val_dd_list2, function(w) {
                      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                        2)])
                      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                        ])], n.components = 1)
                      rr <- as.numeric(ff$v.pred)
                      rr2 <- cbind(w[, 1:2], RS = rr)
                      return(rr2)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "SuperPC")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "SuperPC")]] = list(fit, cv.fit)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "SuperPC")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "GBM") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                    best <- which.min(fit$cv.error)
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "GBM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "GBM")]] = list(fit = fit, best = best)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "GBM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "survivalsvm") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "survival-SVM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "survival-SVM")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "survival-SVM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Ridge") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "Ridge")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "Ridge")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "Ridge")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Lasso") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                    importance = T, proximity = T, forest = T, seed = seed)
                  rid <- var.select(object = fit, conservative = "high")
                  rid <- rid$topvars
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("RSF + ", "Lasso")
                    result <- rbind(result, cc)
                    ml.res[[paste0("RSF + ", "Lasso")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("RSF + ", "Lasso")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2")
                  }
                }
            }
            else if (double_ml1 == "StepCox" & double_ml2 %in% c("CoxBoost", "Enet", "RSF", "plsRcox", "superpc", 
                "GBM", "survivalsvm", "Ridge", "Lasso") & direction_for_stepcox %in% c("both", "backward", "forward")) {
                if (double_ml2 == "CoxBoost") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    message(paste0("---3.StepCox", "[", direction_for_stepcox, "]", " + CoxBoost ---"))
                    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                      2)]), trace = TRUE, start.penalty = 500, parallel = T)
                    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                      penalty = pen$penalty)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                        2], type = "lp")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + CoxBoost")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + CoxBoost")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + CoxBoost")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Enet") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha_for_Enet, nfolds = 10)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                        s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + Enet", "[α=", alpha_for_Enet, 
                      "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Enet", "[α=", alpha_for_Enet, 
                      "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Enet", "[α=", alpha_for_Enet, 
                      "]")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "RSF") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                      importance = T, proximity = T, forest = T, seed = seed)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + RSF")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + RSF")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + RSF")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "plsRcox") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                      nt = 10, verbose = FALSE)
                    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + plsRcox")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + plsRcox")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + plsRcox")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "superpc") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                      featurenames = colnames(est_dd2)[-c(1, 2)])
                    set.seed(seed)
                    fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                    repeat {
                      tryCatch({
                        cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                          max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                        break
                      }, error = function(e) {
                        cat("Error:", conditionMessage(e), "\n")
                        cat("Retrying...\n")
                        invisible(NULL)
                      })
                    }
                    rs <- lapply(val_dd_list2, function(w) {
                      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                        2)])
                      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                        ])], n.components = 1)
                      rr <- as.numeric(ff$v.pred)
                      rr2 <- cbind(w[, 1:2], RS = rr)
                      return(rr2)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + SuperPC")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + SuperPC")]] = list(fit = fit, 
                      cv.fit = cv.fit)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + SuperPC")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "GBM") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                    best <- which.min(fit$cv.error)
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + GBM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + GBM")]] = list(fit = fit, best = best)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + GBM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "survivalsvm") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + survival-SVM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + survival-SVM")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + survival-SVM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Ridge") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + Ridge")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Ridge")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Ridge")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Lasso") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " ", direction_for_stepcox, " + ", double_ml2, " ---"))
                  fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction_for_stepcox)
                  rid <- names(coef(fit))
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("StepCox", "[", direction_for_stepcox, "]", " + Lasso")
                    result <- rbind(result, cc)
                    ml.res[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Lasso")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("StepCox", "[", direction_for_stepcox, "]", " + Lasso")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by StepCox, the first machine learning algorithm, is less than 2")
                  }
                }
            }
            else if (double_ml1 == "CoxBoost" & double_ml2 %in% c("Enet", "StepCox", "plsRcox", "superpc", "GBM", 
                "survivalsvm", "Ridge", "Lasso")) {
                if (double_ml2 == "Enet" & alpha_for_Enet %in% seq(0.1, 0.9, 0.1)) {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ", alpha_for_Enet, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha_for_Enet, nfolds = 10)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "link", newx = as.matrix(x[, -c(1, 2)]), 
                        s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost", " + Enet", "[α=", alpha_for_Enet, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost", " + Enet", "[α=", alpha_for_Enet, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost", " + Enet", "[α=", alpha_for_Enet, "]")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "GBM") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                    best <- which.min(fit$cv.error)
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "GBM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "GBM")]] = list(fit = fit, best = best)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "GBM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Lasso") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "Lasso")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "Lasso")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "Lasso")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "plsRcox") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                      nt = 10, verbose = FALSE)
                    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "plsRcox")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "plsRcox")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "plsRcox")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "Ridge") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    x1 <- as.matrix(est_dd2[, rid])
                    x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
                    set.seed(seed)
                    fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 0, type.measure = "class")
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "response", newx = as.matrix(x[, -c(1, 
                        2)]), s = fit$lambda.min)))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "Ridge")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "Ridge")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "Ridge")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "StepCox" & direction_for_stepcox %in% c("both", "backward", "forward")) {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ", direction_for_stepcox, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction_for_stepcox)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "StepCox", "[", direction_for_stepcox, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "StepCox", "[", direction_for_stepcox, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "StepCox", "[", direction_for_stepcox, "]")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "superpc") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  set.seed(seed)
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                      featurenames = colnames(est_dd2)[-c(1, 2)])
                    set.seed(seed)
                    fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                    repeat {
                      tryCatch({
                        cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                          max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                        break
                      }, error = function(e) {
                        cat("Error:", conditionMessage(e), "\n")
                        cat("Retrying...\n")
                        invisible(NULL)
                      })
                    }
                    rs <- lapply(val_dd_list2, function(w) {
                      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                        2)])
                      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                        ])], n.components = 1)
                      rr <- as.numeric(ff$v.pred)
                      rr2 <- cbind(w[, 1:2], RS = rr)
                      return(rr2)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "SuperPC")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "SuperPC")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "survivalsvm") {
                  set.seed(seed)
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), 
                    trace = TRUE, start.penalty = 500, parallel = T)
                  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), maxstepno = 500, 
                    K = 10, type = "verweij", penalty = pen$penalty)
                  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                    penalty = pen$penalty)
                  rid <- as.data.frame(coef(fit))
                  rid$id <- rownames(rid)
                  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("CoxBoost + ", "survival-SVM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("CoxBoost + ", "survival-SVM")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("CoxBoost + ", "survival-SVM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by CoxBoost, the first machine learning algorithm, is less than 2")
                  }
                }
            }
            else if (double_ml1 == "Lasso" & double_ml2 %in% c("StepCox", "CoxBoost", "plsRcox", "superpc", "GBM", 
                "survivalsvm", "RSF")) {
                if (double_ml2 == "CoxBoost") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    pen <- optimCoxBoostPenalty(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 
                      2)]), trace = TRUE, start.penalty = 500, parallel = T)
                    cv.res <- cv.CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
                    fit <- CoxBoost(est_dd2[, "OS.time"], est_dd2[, "OS"], as.matrix(est_dd2[, -c(1, 2)]), stepno = cv.res$optimal.step, 
                      penalty = pen$penalty)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1], newstatus = x[, 
                        2], type = "lp")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "CoxBoost")
                    result <- rbind(result, cc)
                    ml.res[[paste0("Lasso + ", "CoxBoost")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "CoxBoost")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "GBM") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = 10000, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
                    best <- which.min(fit$cv.error)
                    set.seed(seed)
                    fit <- gbm(formula = Surv(OS.time, OS) ~ ., data = est_dd2, distribution = "coxph", n.trees = best, 
                      interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 8)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = "link")))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "GBM")
                    result <- rbind(result, cc)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "GBM")]] = rs
                    ml.res[[paste0("Lasso + ", "GBM")]] = list(fit = fit, best = best)
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "plsRcox") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), 
                      nt = 10, verbose = FALSE)
                    fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "plsRcox")
                    result <- rbind(result, cc)
                    ml.res[[paste0("Lasso + ", "plsRcox")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "plsRcox")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "RSF") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  rid <- rid[-1]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    set.seed(seed)
                    fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd2, ntree = 1000, nodesize = rf_nodesize, splitrule = "logrank", 
                      importance = T, proximity = T, forest = T, seed = seed)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso", " + RSF")
                    result <- rbind(result, cc)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "RSF")]] = rs
                    ml.res[[paste0("Lasso", " + RSF")]] = fit
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "StepCox" & direction_for_stepcox %in% c("both", "backward", "forward")) {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ", direction_for_stepcox, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit <- step(coxph(Surv(OS.time, OS) ~ ., est_dd2), direction = direction_for_stepcox)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = predict(fit, type = "risk", newdata = x))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "StepCox", "[", direction_for_stepcox, "]")
                    result <- rbind(result, cc)
                    ml.res[[paste0("Lasso + ", "StepCox", "[", direction_for_stepcox, "]")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "StepCox", "[", direction_for_stepcox, "]")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "superpc") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS, 
                      featurenames = colnames(est_dd2)[-c(1, 2)])
                    set.seed(seed)
                    fit <- superpc.train(data = data, type = "survival", s0.perc = 0.5)
                    repeat {
                      tryCatch({
                        cv.fit <- superpc.cv(fit, data, n.threshold = 20, n.fold = 10, n.components = 3, min.features = 2, 
                          max.features = nrow(data$x), compute.fullcv = TRUE, compute.preval = TRUE)
                        break
                      }, error = function(e) {
                        cat("Error:", conditionMessage(e), "\n")
                        cat("Retrying...\n")
                        invisible(NULL)
                      })
                    }
                    rs <- lapply(val_dd_list2, function(w) {
                      test <- list(x = t(w[, -c(1, 2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1, 
                        2)])
                      ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, 
                        ])], n.components = 1)
                      rr <- as.numeric(ff$v.pred)
                      rr2 <- cbind(w[, 1:2], RS = rr)
                      return(rr2)
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "SuperPC")
                    result <- rbind(result, cc)
                    ml.res[[paste0("Lasso + ", "SuperPC")]] = list(fit = fit, cv.fit = cv.fit)
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "SuperPC")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
                else if (double_ml2 == "survivalsvm") {
                  message(paste0("--- ", double_ml1, " + ", double_ml2, " ---"))
                  x1 <- as.matrix(est_dd[, pre_var])
                  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
                  set.seed(seed)
                  fit = cv.glmnet(x1, x2, nfold = 10, family = "cox", alpha = 1, type.measure = "class")
                  fit$lambda.min
                  myCoefs <- coef(fit, s = "lambda.min")
                  rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0)]
                  if (length(rid) > 1) {
                    est_dd2 <- train_data[, c("OS.time", "OS", rid)]
                    val_dd_list2 <- lapply(list_train_vali_Data, function(x) {
                      x[, c("OS.time", "OS", rid)]
                    })
                    fit = survivalsvm(Surv(OS.time, OS) ~ ., data = est_dd2, gamma.mu = 1)
                    rs <- lapply(val_dd_list2, function(x) {
                      cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))
                    })
                    cc <- data.frame(Cindex = sapply(rs, function(x) {
                      as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])
                    })) %>% rownames_to_column("ID")
                    cc$Model <- paste0("Lasso + ", "survival-SVM")
                    result <- rbind(result, cc)
                    ml.res[[paste0("Lasso + ", "survival-SVM")]] = fit
                    rs = returnIDtoRS(rs.table.list = rs, rawtableID = list_train_vali_Data)
                    riskscore[[paste0("Lasso + ", "survival-SVM")]] = rs
                    result.cindex.fit = list(Cindex.res = result, ml.res = ml.res, riskscore = riskscore, Sig.genes = pre_var)
                    return(result.cindex.fit)
                  }
                  else {
                    warning("The number of seleted candidate gene by Lasso, the first machine learning algorithm, is less than 2")
                  }
                }
            }
        }
        message("--- The analysis has been completed ---")
    }
    else {
        print("Please provide the full parameters, verify that the column names of the supplied queues match the conditions\n          (ID, OS.time, OS for the first through third columns, respectively) ,  verify there exist common genes in\n          all cohorts")
    }
}

do_bulk_TimiGP <- function(dat_exp, pdata, dataset = "Bindea2013_Cancer", log = FALSE, GMNorm = FALSE, QV = TRUE) {
    library(ggplot2)
    library(TimiGP)
    library(qs)
    if (dataset == "Bindea2013_Cancer") {
        prad_restore("input_068")
        geneset <- CellType_Bindea2013_cancer
    }
    else if (dataset == "Charoentong2017") {
        prad_restore("input_069")
        geneset <- CellType_Charoentong2017_Bindea2013_Xu2018_Immune
    }
    else if (dataset == "Newman2015") {
        prad_restore("input_070")
        geneset <- CellType_Newman2015_LM22
    }
    else {
        stop("dataset must in one of Bindea2013_Cancer Charoentong2017 Newman2015")
    }
    marker <- unique(geneset$Gene)
    info <- TimiCheckEvent(pdata)
    com_sam <- intersect(rownames(info), colnames(dat_exp))
    info <- info[com_sam, ]
    dat_exp <- dat_exp[, com_sam]
    rna <- TimiPrePropress(marker = marker, rna = dat_exp, cohort = rownames(info), log = log, GMNorm = GMNorm)
    mps <- TimiGenePair(rna)
    res <- TimiCOX(mps = mps, info = info, p.adj = "BH")
    mps <- res[[1]]
    cox_res <- res[[2]]
    NET <- tryCatch(expr = {
        TimiGeneNetwork(resdata = cox_res, dataset = dataset, export = TRUE, path = "./")
    }, error = function(e) {
        message("TimiGeneNetwork 执行失败: ", e$message)
        return(NULL)
    }, warning = function(w) {
        message("警告: ", w$message)
        invokeRestart("muffleWarning")
    })
    cell_pair <- TimiCellPair(geneset = geneset, core = 1)
    background <- TimiBG(marker.pair = row.names(cox_res))
    if (QV) {
        GP <- rownames(cox_res)[which(cox_res$QV < 0.05)]
    }
    else {
        GP <- rownames(cox_res)[which(cox_res$PV < 0.05)]
    }
    res <- TimiEnrich(gene = GP, background = background, geneset = cell_pair, p.adj = "BH", core = 1)
    score <- TimiFS(res)
    res_TimiGP <- list(res = res, score = score)
    NET <- tryCatch(expr = {
        TimiGeneNetwork(resdata = res, dataset = dataset, export = TRUE, path = "./")
    }, error = function(e) {
        message("TimiGeneNetwork 执行失败: ", e$message)
        return(NULL)
    }, warning = function(w) {
        message("警告: ", w$message)
        invokeRestart("muffleWarning")
    })
    return(res_TimiGP)
    TimiDotplot(resdata = res, select = c(1:10))
    TimiCellChord(resdata = res, dataset = dataset)
    TimiGeneChord(resdata = res, select = 1)
    TimiFSBar(score, select = 1:10)
}

do_bulk_gsva_ssgsea <- function(dat_exp, geneset, id) {
    library(tidyverse)
    library(GSVA)
    gene_set <- list(pathway = geneset)
    names(gene_set)[1] <- id
    gsvaParam <- gsvaParam(as.matrix(dat_exp), geneSets = gene_set, minSize = 1)
    my_gsva <- gsva(gsvaParam)
    ssgseaParam <- ssgseaParam(as.matrix(dat_exp), geneSets = gene_set, minSize = 1)
    my_ssgsea <- gsva(ssgseaParam)
    dat <- list(res_gsva = my_gsva, res_ssgsea = my_ssgsea)
    dat1 <- dat$res_gsva %>% t() %>% as.data.frame()
    dat2 <- dat$res_ssgsea %>% t() %>% as.data.frame()
    dat <- list(res_gsva = dat1, res_ssgsea = dat2)
    return(dat)
}

do_cor_test_col <- function(data, variable, methods = "pearson") {
    x = as.numeric(data[, variable])
    varible_others = colnames(data)[colnames(data) != variable]
    outTab = data.frame()
    for (varible_another in varible_others) {
        y = as.numeric(data[, varible_another])
        corT = cor.test(x, y, method = methods)
        cor = corT$estimate
        pvalue = corT$p.value
        outTab = rbind(outTab, data.frame(variable = variable, varible_others = varible_another, cor = cor, pvalue = pvalue))
    }
    return(outTab)
}

do_sc_AUCell <- function(scRNA, markers, cores_use = 5, name_use) {
    library(Seurat)
    library(AUCell)
    markers = intersect(markers, rownames(scRNA))
    geneSets = list(Custom = markers)
    cells_rankings = AUCell_buildRankings(as.matrix(GetAssayData(object = scRNA, assay = "RNA")), nCores = cores_use)
    merge_row <- function(data1, data2) {
        samesample = intersect(rownames(data1), rownames(data2))
        data1 = data1[samesample, , drop = FALSE]
        data2 = data2[samesample, , drop = FALSE]
        data3 = cbind(data1, data2)
        return(data3)
    }
    cells_AUC = AUCell_calcAUC(geneSets, cells_rankings, nCores = cores_use, aucMaxRank = nrow(cells_rankings) * 
        0.05)
    AUC_Score = cells_AUC@assays@data@listData[["AUC"]]
    AUC_Score = t(AUC_Score)
    colnames(AUC_Score) = "AUC"
    scRNA@meta.data = merge_row(scRNA@meta.data, AUC_Score)
    colnames(scRNA@meta.data)[ncol(scRNA@meta.data)] = name_use
    return(scRNA)
}

do_sc_CytoTRACE2 <- function(sc, species, version) {
    library(CytoTRACE2)
    library(tidyverse)
    library(Seurat)
    if (version == "V5") {
        sc <- JoinLayers(sc)
        sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")
    }
    else {
        sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")
    }
    res <- cytotrace2(input = sc, ncores = 1, species = species, is_seurat = T, slot_type = "counts", seed = 5201314, 
        batch_size = 10000, smooth_batch_size = 1000, parallelize_models = TRUE, parallelize_smoothing = TRUE, 
        max_pcs = 200, full_model = FALSE)
    return(res)
}

do_sc_GIS <- function(sc = sc, version = "v5") {
    library(genomicInstability)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    if (version == "v5") {
        sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")
    }
    else {
        sc <- sc
    }
    counts_matrix <- GetAssayData(sc, assay = "RNA", slot = "counts") %>% as.matrix()
    gene.df <- bitr(rownames(counts_matrix), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
    counts_matrix <- counts_matrix[gene.df$SYMBOL, ]
    rownames(counts_matrix) <- gene.df$ENTREZID
    cnv <- inferCNV(as.matrix(counts_matrix))
    cnv <- genomicInstabilityScore(cnv)
    dat <- data.frame(row.names = names(cnv[["gis"]]), GIS = cnv[["gis"]])
    sc <- AddMetaData(sc, metadata = dat)
    return(sc)
}

do_sc_Pipeline <- function(scRNA, method = "LogNormalize", resolution = 1, integration = NULL) {
    if (method == "LogNormalize") {
        scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
        scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 3000)
        scRNA <- ScaleData(scRNA)
        scRNA <- RunPCA(scRNA, npcs = 30, pc.genes = VariableFeatures(object = scRNA), verbose = FALSE)
        pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        co1 <- which(cumu > 90 & pct < 5)[1]
        co1
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
        co2
        pcSelect <- min(co1, co2)
        print(pcSelect)
        if (is.null(integration)) {
            scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "pca", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "pca", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else if (integration == "harmony") {
            library(harmony)
            scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", reduction.use = "pca", assay.use = "RNA", 
                max.iter.harmony = 20)
            scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "harmony", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "harmony", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else if (integration == "css") {
            library(simspec)
            scRNA <- cluster_sim_spectrum(object = scRNA, label_tag = "orig.ident", spectrum_type = "corr_ztransform", 
                corr_method = "spearman")
            scRNA <- FindNeighbors(scRNA, reduction = "css", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "css", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "css", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else {
            stop("Invalid integration method specified. Please choose from NULL, 'harmony', or 'css'.")
        }
    }
    else if (method == "SCT") {
        scRNA <- SCTransform(scRNA)
        scRNA <- RunPCA(scRNA, assay = "SCT", npcs = 30, pc.genes = VariableFeatures(object = scRNA), verbose = FALSE)
        pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        co1 <- which(cumu > 90 & pct < 5)[1]
        co1
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
        co2
        pcSelect <- min(co1, co2)
        print(pcSelect)
        if (is.null(integration)) {
            scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "pca", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "pca", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else if (integration == "harmony") {
            library(harmony)
            scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", reduction.use = "pca", assay.use = "RNA", 
                max.iter.harmony = 20)
            scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "harmony", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "harmony", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else if (integration == "css") {
            library(simspec)
            scRNA <- cluster_sim_spectrum(object = scRNA, label_tag = "orig.ident", spectrum_type = "corr_ztransform", 
                corr_method = "spearman")
            scRNA <- FindNeighbors(scRNA, reduction = "css", dims = 1:pcSelect)
            scRNA <- FindClusters(object = scRNA, reduction = "css", resolution = resolution)
            scRNA <- RunUMAP(object = scRNA, reduction = "css", dims = 1:pcSelect, check_duplicates = FALSE)
            return(scRNA)
        }
        else {
            stop("Invalid integration method specified. Please choose from NULL, 'harmony', or 'css'.")
        }
    }
}

do_sc_ROE <- function(sc, col_celltype, col_group, col_sample, levels) {
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(tidyverse)
    library(circlize)
    library(Startrac)
    library(ggplot2)
    library(tictoc)
    library(ggpubr)
    library(sscVis)
    dat_use <- sc@meta.data
    colnames(dat_use)[colnames(dat_use) == col_celltype] <- "majorCluster"
    colnames(dat_use)[colnames(dat_use) == col_group] <- "loc"
    colnames(dat_use)[colnames(dat_use) == col_sample] <- "col_sample"
    if (!is.null(levels)) {
        dat_use$loc <- factor(dat_use$loc, levels = levels)
    }
    else {
        dat_use <- dat_use
    }
    R_oe <- calTissueDist(dat_use, byPatient = FALSE, colname.cluster = "majorCluster", colname.patient = "col_sample", 
        colname.tissue = "loc", method = "chisq", min.rowSum = 0)
    col_fun <- colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)), c("#00599F", "#f6f8e6", "#D01910"))
    mycolor <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", 
        "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF")
    sample_groups <- paste0("Group", 1:ncol(R_oe))
    names(sample_groups) <- colnames(R_oe)
    col_Group <- setNames(mycolor[1:ncol(R_oe)], sample_groups)
    col_ha <- HeatmapAnnotation(Group = sample_groups[colnames(R_oe)], col = list(Group = col_Group), show_legend = FALSE)
    Heatmap(as.matrix(R_oe), top_annotation = col_ha, rect_gp = gpar(col = "black", lwd = 1), show_heatmap_legend = TRUE, 
        cluster_rows = TRUE, cluster_columns = FALSE, row_names_side = "right", column_names_side = "top", show_column_names = TRUE, 
        show_row_names = TRUE, col = col_fun, row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 12), 
        heatmap_legend_param = list(title = "Ro/e", at = c(0, max(R_oe)), labels = c("Min", "Max")), cell_fun = function(j, 
            i, x, y, width, height, fill) {
            value <- R_oe[i, j]
            symbol <- if (value == 0) {
                "−"
            }
            else if (value > 0 & value < 0.2) {
                "+/−"
            }
            else if (value >= 0.2 & value <= 0.8) {
                "+"
            }
            else if (value > 0.8 & value <= 1) {
                "++"
            }
            else if (value > 1) {
                "+++"
            }
            grid.text(symbol, x, y, gp = gpar(fontsize = 12, col = "black"))
        })
}

do_sc_bulk_Progeny <- function(dat_exp, sc = NULL, dat_type = "Bulk", org = "Human") {
    library(Seurat)
    library(progeny)
    library(tidyverse)
    if (dat_type == "Bulk") {
        progeny_res <- progeny(as.matrix(dat_exp), scale = TRUE, organism = org, top = 100, perm = 1, verbose = FALSE)
        progeny_res <- as.data.frame(progeny_res)
        return(progeny_res)
    }
    else if (dat_type == "sc") {
        dat_exp <- GetAssayData(sc, assay = "RNA", slot = "data") %>% as.matrix()
        progeny_res <- progeny(dat_exp, scale = TRUE, organism = org, top = 100, perm = 1, verbose = FALSE)
        progeny_res <- as.data.frame(progeny_res)
        sc <- AddMetaData(sc, progeny_res)
        return(sc)
    }
}

do_sc_p_n_gene <- function(sc, gene, expression) {
    library(Seurat)
    df <- FetchData(sc, vars = gene)
    df$status <- ifelse(df[, 1] > expression, paste0(gene, "+"), paste0(gene, "-"))
    sc <- AddMetaData(object = sc, metadata = df)
    sc$status <- factor(sc$status, levels = c(paste0(gene, "+"), paste0(gene, "-")))
    Idents(sc) <- "status"
    print(table(sc@meta.data$status))
    return(sc)
}

lzq_GSEA.integrated_1 <- function(genes, gene.type = "SYMBOL", organism = "Human", GO.ont = "BP", KEGG.use.internal.data = FALSE, 
    perform.WikiPathways = TRUE, perform.Reactome = TRUE, perform.MsigDB = TRUE, MsigDB.category = "H", perform.disease.ontoloty = TRUE, 
    perform.Cancer.Gene.Network = TRUE, perform.DisGeNET = TRUE, perform.CellMarker = TRUE, perform.CMAP = TRUE, 
    pvalue.cutoff = 0.05, padjust.method = "BH", min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3) {
    lzq_updateSymbolforDL_1 <- function(data, unmapGene_keep = F) {
        ann <- HGNChelper::checkGeneSymbols(rownames(data), unmapped.as.na = unmapGene_keep) %>% stats::na.omit() %>% 
            dplyr::select(-Approved)
        colnames(ann) <- c("ID", "Symbol")
        if (sum(grepl("///", ann$Symbol)) > 0) {
            ann1 <- ann[!grepl("///", ann$Symbol), ]
            ann2 <- ann[grepl("///", ann$Symbol), ]
            ann2 <- lapply(ann2$Symbol, function(x) {
                l <- unlist(strsplit(x, " /// "))
                tmp <- ann2[ann2$Symbol == x, ] %>% dplyr::distinct(ID, Symbol, .keep_all = T)
                tmp2 <- lapply(l, function(d) {
                  kk <- tmp
                  kk$Symbol <- d
                  return(kk)
                }) %>% Reduce(rbind, .)
                return(tmp2)
            }) %>% Reduce(rbind, .) %>% dplyr::distinct(ID, Symbol, .keep_all = T)
            ann <- rbind(ann1, ann2) %>% dplyr::distinct(ID, Symbol, .keep_all = T)
        }
        data2 <- merge(ann, data, by.x = 1, by.y = 0)[, -1] %>% dplyr::group_by(Symbol) %>% dplyr::summarise_all(mean) %>% 
            tibble::column_to_rownames("Symbol")
        return(data2)
    }
    d <- as.data.frame(genes)
    if (gene.type == "SYMBOL") {
        cat("+++ Updating gene symbols...")
        cat("\n")
        d2 <- suppressWarnings(lzq_updateSymbolforDL_1(d, unmapGene_keep = T)) %>% dplyr::arrange(dplyr::desc(genes))
    }
    if (organism == "Human") {
        OrgDb <- "org.Hs.eg.db"
    }
    if (organism == "Mouse") {
        OrgDb <- "org.Mm.eg.db"
    }
    if (gene.type != "ENTREZID") {
        cat(paste0("+++ Transforming ", gene.type, " to ENTREZID..."))
        cat("\n")
        d3 <- suppressWarnings(clusterProfiler::bitr(rownames(d2), fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))
        d3 <- merge(d3, d2, by.x = 1, by.y = 0) %>% dplyr::arrange(dplyr::desc(genes))
        genes <- d3$genes
        names(genes) <- d3$ENTREZID
    }
    else {
        genes <- d2$genes
        names(genes) <- rownames(d2)
    }
    cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
    cat("\n")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::gseGO(geneList = genes, OrgDb = OrgDb, keyType = "ENTREZID", 
        ont = GO.ont, pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
        maxGSSize = max.Geneset.Size, eps = 1e-100, seed = T, verbose = F)))
    ego <- DOSE::setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
    cat("+++ Performing KEGG enrichment...")
    cat("\n")
    if (organism == "Human") {
        KEGG.organism <- "hsa"
    }
    if (organism == "Mouse") {
        KEGG.organism <- "mmu"
    }
    ekegg <- suppressWarnings(suppressMessages(clusterProfiler::gseKEGG(geneList = genes, organism = KEGG.organism, 
        minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, 
        eps = 1e-100, seed = T, verbose = F, use_internal_data = KEGG.use.internal.data)))
    ekegg <- DOSE::setReadable(ekegg, OrgDb = OrgDb, keyType = "ENTREZID")
    res <- list(GO = ego, KEGG = ekegg)
    if (perform.WikiPathways) {
        cat("+++ Performing WikiPathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            WikiPathways.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            WikiPathways.organism <- "Mus musculus"
        }
        eWP <- suppressWarnings(clusterProfiler::gseWP(geneList = genes, organism = WikiPathways.organism, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, 
            seed = T, verbose = F))
        res[["WikiPathways"]] <- DOSE::setReadable(eWP, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.Reactome) {
        cat("+++ Performing Reactome pathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            Reactome.organism <- "human"
        }
        if (organism == "Mouse") {
            Reactome.organism <- "mouse"
        }
        eRP <- suppressWarnings(ReactomePA::gsePathway(geneList = genes, organism = Reactome.organism, minGSSize = min.Geneset.Size, 
            maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, 
            seed = T, verbose = F))
        res[["Reactome"]] <- DOSE::setReadable(eRP, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.disease.ontoloty) {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::gseDO(geneList = genes, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = T, verbose = F))
        res[["DiseaseOntology"]] <- DOSE::setReadable(eDO, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.Cancer.Gene.Network) {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::gseNCG(geneList = genes, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = T, verbose = F))
        res[["CancerGeneNetwork"]] <- DOSE::setReadable(eNCG, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.DisGeNET) {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::gseDGN(geneList = genes, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = TRUE, verbose = FALSE))
        res[["DisGeNET"]] <- DOSE::setReadable(eDGN, OrgDb = OrgDb, keyType = "ENTREZID")
    }
    if (perform.CellMarker) {
        cat("+++ Performing CellMarker enrichment...")
        cat("\n")
        if (organism == "Human") {
            cell_marker_data <- suppressMessages(prad_input("input_071"))
        }
        if (organism == "Mouse") {
            cell_marker_data <- suppressMessages(prad_input("input_072"))
        }
        cells <- cell_marker_data %>% dplyr::select(cellName, geneID) %>% dplyr::mutate(geneID = strsplit(geneID, 
            ", ")) %>% tidyr::unnest(cols = geneID)
        eCM <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = TRUE, verbose = FALSE, 
            TERM2GENE = cells))
        res[["CellMarker.Enrichment"]] <- list(Enrich.results = DOSE::setReadable(eCM, OrgDb = OrgDb, keyType = "ENTREZID"), 
            Genesets = cells)
    }
    if (perform.MsigDB) {
        cat(paste0("+++ Performing MsigDB-", MsigDB.category, " enrichment..."))
        cat("\n")
        if (organism == "Human") {
            MsigDB.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            MsigDB.organism <- "Mus musculus"
        }
        mall <- msigdbr::msigdbr(species = MsigDB.organism)
        if (MsigDB.category == "All" | MsigDB.category == "ALL") {
            mg <- mall %>% dplyr::select(gs_name, ncbi_gene)
        }
        else {
            mg <- msigdbr::msigdbr(species = MsigDB.organism, category = MsigDB.category) %>% dplyr::select(gs_name, 
                ncbi_gene)
        }
        eMSIG <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = T, verbose = F, 
            TERM2GENE = mg))
        res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = DOSE::setReadable(eMSIG, OrgDb = OrgDb, 
            keyType = "ENTREZID"), Genesets = mg)
    }
    if (perform.CMAP) {
        if (organism != "Human") {
            stop("CMAP only supports organism = Human!")
        }
        cat("+++ Performing CMAP enrichment...")
        cat("\n")
        eCMAP <- suppressWarnings(clusterProfiler::GSEA(geneList = genes, minGSSize = CMAP.min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method, eps = 1e-100, seed = T, verbose = F, 
            TERM2GENE = CMAPfromDSEATM[, seq_len(2)]))
        res[["CMAP"]] <- list(Enrich.results = DOSE::setReadable(eCMAP, OrgDb = OrgDb, keyType = "ENTREZID"), Genesets = CMAPfromDSEATM)
    }
    nsig <- ifelse(is.null(nrow(res$GO)), 0, nrow(res$GO)) + ifelse(is.null(nrow(res$KEGG)), 0, nrow(res$KEGG)) + 
        ifelse(is.null(nrow(res$Module.KEGG)), 0, nrow(res$Module.KEGG)) + ifelse(is.null(nrow(res$WikiPathways)), 
        0, nrow(res$WikiPathways)) + ifelse(is.null(nrow(res$Reactome)), 0, nrow(res$Reactome)) + ifelse(is.null(nrow(res$DiseaseOntology)), 
        0, nrow(res$DiseaseOntology)) + ifelse(is.null(nrow(res$CancerGeneNetwork)), 0, nrow(res$CancerGeneNetwork)) + 
        ifelse(is.null(nrow(res$DisGeNET)), 0, nrow(res$DisGeNET)) + ifelse(is.null(nrow(res$CellMarker.Enrichment$Enrich.results)), 
        0, nrow(res$CellMarker.Enrichment$Enrich.results)) + ifelse(is.null(nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)), 
        0, nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)) + ifelse(is.null(nrow(res$CMAP$Enrich.results)), 
        0, nrow(res$CMAP$Enrich.results))
    cat(paste0("+++ ", nsig, " significant terms were detected..."))
    cat("\n")
    cat("+++ Done!")
    return(res)
}

lzq_ORA.integrated_1 <- function(genes, background.genes = NULL, gene.type = "SYMBOL", organism = "Human", GO.ont = "BP", 
    KEGG.use.internal.data = FALSE, perform.WikiPathways = TRUE, perform.Reactome = TRUE, perform.MsigDB = TRUE, 
    MsigDB.category = "H", perform.disease.ontoloty = TRUE, perform.Cancer.Gene.Network = TRUE, perform.DisGeNET = TRUE, 
    perform.CellMarker = TRUE, perform.CMAP = TRUE, pvalue.cutoff = 0.05, qvalue.cutoff = 0.05, padjust.method = "BH", 
    min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3) {
    lzq_updateSymbol_1 <- function(genes, unmapGene_keep = F) {
        ann <- HGNChelper::checkGeneSymbols(genes, unmapped.as.na = unmapGene_keep) %>% stats::na.omit() %>% dplyr::select(-Approved)
        colnames(ann) <- c("ID", "Symbol")
        if (sum(grepl("///", ann$Symbol)) > 0) {
            ann1 <- ann[!grepl("///", ann$Symbol), ]
            ann2 <- ann[grepl("///", ann$Symbol), ]
            ann2 <- lapply(ann2$Symbol, function(x) {
                l <- unlist(strsplit(x, " /// "))
                tmp <- ann2[ann2$Symbol == x, ] %>% dplyr::distinct(ID, Symbol, .keep_all = T)
                tmp2 <- lapply(l, function(d) {
                  kk <- tmp
                  kk$Symbol <- d
                  return(kk)
                }) %>% Reduce(rbind, .)
                return(tmp2)
            }) %>% Reduce(rbind, .) %>% dplyr::distinct(ID, Symbol, .keep_all = T)
            ann <- rbind(ann1, ann2) %>% dplyr::distinct(ID, Symbol, .keep_all = T)
        }
        return(ann)
    }
    if (gene.type == "SYMBOL") {
        cat("+++ Updating gene symbols...")
        cat("\n")
        genes <- suppressWarnings(lzq_updateSymbol_1(genes = genes, unmapGene_keep = T)[[2]])
        if (!is.null(background.genes)) {
            background.genes <- suppressWarnings(lzq_updateSymbol_1(genes = stats::na.omit(background.genes), unmapGene_keep = T)[[2]])
        }
    }
    if (organism == "Human") {
        OrgDb <- "org.Hs.eg.db"
    }
    if (organism == "Mouse") {
        OrgDb <- "org.Mm.eg.db"
    }
    if (gene.type != "ENTREZID") {
        cat(paste0("+++ Transforming ", gene.type, " to ENTREZID..."))
        cat("\n")
        genes <- suppressWarnings(clusterProfiler::bitr(genes, fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))[, 
            "ENTREZID"]
        if (!is.null(background.genes)) {
            background.genes <- suppressWarnings(clusterProfiler::bitr(background.genes, fromType = gene.type, 
                toType = "ENTREZID", OrgDb = OrgDb))[, "ENTREZID"]
        }
    }
    cat(paste0("+++ Performing GO-", GO.ont, " enrichment..."))
    cat("\n")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::enrichGO(gene = genes, OrgDb = OrgDb, keyType = "ENTREZID", 
        ont = GO.ont, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, 
        minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, readable = TRUE, universe = background.genes)))
    if (!is.null(ego)) {
        ego <- lzq_getEF(ego)
    }
    cat("+++ Performing KEGG enrichment...")
    cat("\n")
    if (organism == "Human") {
        KEGG.organism <- "hsa"
    }
    if (organism == "Mouse") {
        KEGG.organism <- "mmu"
    }
    ekegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichKEGG(gene = genes, organism = KEGG.organism, 
        pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, 
        maxGSSize = max.Geneset.Size, use_internal_data = KEGG.use.internal.data, universe = background.genes)))
    ekegg <- DOSE::setReadable(ekegg, OrgDb = OrgDb, keyType = "ENTREZID")
    if (!is.null(ekegg)) {
        ekegg <- lzq_getEF(ekegg)
    }
    res <- list(GO = ego, KEGG = ekegg)
    if (perform.WikiPathways) {
        cat("+++ Performing WikiPathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            WikiPathways.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            WikiPathways.organism <- "Mus musculus"
        }
        eWP <- suppressWarnings(clusterProfiler::enrichWP(gene = genes, organism = WikiPathways.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            universe = background.genes))
        eWP <- DOSE::setReadable(eWP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eWP)) {
            res[["WikiPathways"]] <- lzq_getEF(eWP)
        }
        else {
            res[["WikiPathways"]] <- eWP
        }
    }
    if (perform.Reactome) {
        cat("+++ Performing Reactome pathways enrichment...")
        cat("\n")
        if (organism == "Human") {
            Reactome.organism <- "human"
        }
        if (organism == "Mouse") {
            Reactome.organism <- "mouse"
        }
        eRP <- suppressWarnings(ReactomePA::enrichPathway(gene = genes, organism = Reactome.organism, pvalueCutoff = pvalue.cutoff, 
            qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, 
            readable = TRUE, universe = background.genes))
        eRP <- DOSE::setReadable(eRP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eRP)) {
            res[["Reactome"]] <- lzq_getEF(eRP)
        }
        else {
            res[["Reactome"]] <- eRP
        }
    }
    if (perform.disease.ontoloty) {
        cat("+++ Performing Disease Ontoloty enrichment...")
        cat("\n")
        eDO <- suppressWarnings(DOSE::enrichDO(gene = genes, ont = "HDO", pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, readable = TRUE, 
            universe = background.genes))
        eDO <- DOSE::setReadable(eDO, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eDO)) {
            res[["DiseaseOntology"]] <- lzq_getEF(eDO)
        }
        else {
            res[["DiseaseOntology"]] <- eDO
        }
    }
    if (perform.Cancer.Gene.Network) {
        cat("+++ Performing Cancer Gene Network enrichment...")
        cat("\n")
        eNCG <- suppressWarnings(DOSE::enrichNCG(gene = genes, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, readable = T, 
            universe = background.genes))
        eNCG <- DOSE::setReadable(eNCG, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eNCG)) {
            res[["CancerGeneNetwork"]] <- lzq_getEF(eNCG)
        }
        else {
            res[["CancerGeneNetwork"]] <- eNCG
        }
    }
    if (perform.DisGeNET) {
        cat("+++ Performing DisGeNET enrichment...")
        cat("\n")
        eDGN <- suppressWarnings(DOSE::enrichDGN(gene = genes, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, readable = T, 
            universe = background.genes))
        eDGN <- DOSE::setReadable(eDGN, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eDGN)) {
            res[["DisGeNET"]] <- lzq_getEF(eDGN)
        }
        else {
            res[["DisGeNET"]] <- eDGN
        }
    }
    if (perform.CellMarker) {
        cat("+++ Performing CellMarker enrichment...")
        cat("\n")
        if (organism == "Human") {
            cell_marker_data <- suppressMessages(prad_input("input_073"))
        }
        if (organism == "Mouse") {
            cell_marker_data <- suppressMessages(prad_input("input_074"))
        }
        cells <- cell_marker_data %>% dplyr::select(cellName, geneID) %>% dplyr::mutate(geneID = strsplit(geneID, 
            ", ")) %>% tidyr::unnest(cols = geneID)
        eCM <- suppressWarnings(clusterProfiler::enricher(gene = genes, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, TERM2GENE = cells, 
            universe = background.genes))
        eCM <- DOSE::setReadable(eCM, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eCM)) {
            res[["CellMarker.Enrichment"]] <- list(Enrich.results = lzq_getEF(eCM), Genesets = cells)
        }
        else {
            res[["CellMarker.Enrichment"]] <- list(Enrich.results = eCM, Genesets = cells)
        }
    }
    if (perform.MsigDB) {
        cat(paste0("+++ Performing MsigDB-", MsigDB.category, " enrichment..."))
        cat("\n")
        if (organism == "Human") {
            MsigDB.organism <- "Homo sapiens"
        }
        if (organism == "Mouse") {
            MsigDB.organism <- "Mus musculus"
        }
        mall <- msigdbr::msigdbr(species = MsigDB.organism)
        if (MsigDB.category == "All" | MsigDB.category == "ALL") {
            mg <- mall %>% dplyr::select(gs_name, ncbi_gene)
        }
        else {
            mg <- msigdbr::msigdbr(species = MsigDB.organism, category = MsigDB.category) %>% dplyr::select(gs_name, 
                ncbi_gene)
        }
        eMSIG <- suppressWarnings(clusterProfiler::enricher(gene = genes, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size, TERM2GENE = mg, 
            universe = background.genes))
        eMSIG <- DOSE::setReadable(eMSIG, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eMSIG)) {
            res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = lzq_getEF(eMSIG), Genesets = mg)
        }
        else {
            res[[paste0("MsigDB.", MsigDB.category)]] <- list(Enrich.results = eMSIG, Genesets = mg)
        }
    }
    if (perform.CMAP) {
        if (organism != "Human") {
            stop("CMAP only supports organism = Human!")
        }
        cat("+++ Performing CMAP enrichment...")
        cat("\n")
        eCMAP <- suppressWarnings(clusterProfiler::enricher(gene = genes, pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff, 
            pAdjustMethod = padjust.method, minGSSize = CMAP.min.Geneset.Size, maxGSSize = max.Geneset.Size, TERM2GENE = CMAPfromDSEATM[, 
                seq_len(2)], universe = background.genes))
        eCMAP <- DOSE::setReadable(eCMAP, OrgDb = OrgDb, keyType = "ENTREZID")
        if (!is.null(eCMAP)) {
            res[["CMAP"]] <- list(Enrich.results = lzq_getEF(eCMAP), Genesets = CMAPfromDSEATM)
        }
        else {
            res[["CMAP"]] <- list(Enrich.results = eCMAP, Genesets = CMAPfromDSEATM)
        }
    }
    nsig <- ifelse(is.null(nrow(res$GO)), 0, nrow(res$GO)) + ifelse(is.null(nrow(res$KEGG)), 0, nrow(res$KEGG)) + 
        ifelse(is.null(nrow(res$Module.KEGG)), 0, nrow(res$Module.KEGG)) + ifelse(is.null(nrow(res$WikiPathways)), 
        0, nrow(res$WikiPathways)) + ifelse(is.null(nrow(res$Reactome)), 0, nrow(res$Reactome)) + ifelse(is.null(nrow(res$DiseaseOntology)), 
        0, nrow(res$DiseaseOntology)) + ifelse(is.null(nrow(res$CancerGeneNetwork)), 0, nrow(res$CancerGeneNetwork)) + 
        ifelse(is.null(nrow(res$DisGeNET)), 0, nrow(res$DisGeNET)) + ifelse(is.null(nrow(res$CellMarker.Enrichment$Enrich.results)), 
        0, nrow(res$CellMarker.Enrichment$Enrich.results)) + ifelse(is.null(nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)), 
        0, nrow(res[[paste0("MsigDB.", MsigDB.category)]]$Enrich.results)) + ifelse(is.null(nrow(res$CMAP$Enrich.results)), 
        0, nrow(res$CMAP$Enrich.results))
    cat(paste0("+++ ", nsig, " significant terms were detected..."))
    cat("\n")
    cat("+++ Done!")
    return(res)
}

plot_ORA_Bar <- function(dat, color = "#D01910") {
    library(ggplot2)
    dat_plot <- dat
    dat_plot <- dat_plot %>% dplyr::arrange(pvalue) %>% dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
    p <- ggplot(dat_plot) + geom_bar(aes(x = -log10(pvalue), y = Description), stat = "identity", fill = paste0(color, 
        "50"), width = 0.8, position = position_nudge(y = 0.2)) + geom_text(aes(x = 0.05, y = Description, label = str_wrap(Description, 
        width = 60)), hjust = 0, vjust = 1, size = 3.5, lineheight = 1, nudge_y = 0.35, color = "black") + scale_x_continuous(name = "", 
        limits = c(0, max(-log10(dat_plot$pvalue)) * 1.05), expand = c(0, 0.15)) + theme_minimal() + theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(color = "black", size = 12), axis.title = element_text(color = "black", size = 12), 
        axis.title.y = element_blank(), axis.line.x = element_line(color = "black", linewidth = 0.8), axis.ticks.x = element_line(color = "black"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), panel.grid = element_blank(), 
        plot.margin = margin(12, 12, 12, 12), legend.position = "none") + labs(x = "", title = NULL)
    print(p)
}

plot_Time_ROC <- function(dat, target, year = TRUE) {
    library(timeROC)
    library(ggplot2)
    if (year) {
        dat <- dat
    }
    else {
        stop("时间变量必须以年为单位")
    }
    if (all(c("RFS", "RFS.time") %in% colnames(dat))) {
        print("数据具有所分析变量，运行下面代码")
    }
    else {
        stop("数据框缺少RFS或RFS.time列，请检查")
    }
    marker <- dat[, target, drop = T]
    timeroc <- timeROC(T = dat$RFS.time, delta = dat$RFS, marker = marker, cause = 1, weighting = "cox", times = c(1, 
        3, 5), iid = F)
    dat1 <- timeroc$TP %>% as.data.frame()
    dat2 <- timeroc$FP %>% as.data.frame()
    t1 <- data.frame(tp = dat1$`t=1`, fp = dat2$`t=1`) %>% mutate(time = "1-Year")
    t2 <- data.frame(tp = dat1$`t=3`, fp = dat2$`t=3`) %>% mutate(time = "3-Year")
    t3 <- data.frame(tp = dat1$`t=5`, fp = dat2$`t=5`) %>% mutate(time = "5-Year")
    tmp <- rbind(t1, t2, t3)
    p <- ggplot(tmp, aes(fp, tp, color = time)) + geom_line(size = 1) + labs(x = "1-Specificity", y = "Sensitivity", 
        color = NULL) + theme_bw(base_rect_size = 1.5) + geom_abline(slope = 1, color = "grey70") + ggtitle("TCGA Cohort") + 
        theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 15, colour = "darkred", 
            face = "bold"), panel.grid = element_blank(), panel.grid.major = element_line(color = "white", linetype = "dashed"), 
            panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, size = 18, 
                colour = "darkred", face = "bold"), legend.text = element_text(size = 14), legend.position = c(0.99, 
                0.013), legend.justification = c(1, 0)) + scale_color_manual(values = c("#2a9d8f", "#fcbf49", "#eb5e28")) + 
        scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01))
    print(p)
}

plot_expr_cor_Dot_plus <- function(dat, col_cor = "cor", col_pvalue = "pvalue") {
    library(tidyverse)
    library(ggplot2)
    library(dplyr)
    dat_use <- dat
    colnames(dat_use)[colnames(dat_use) == col_pvalue] <- "pvalue"
    colnames(dat_use)[colnames(dat_use) == col_cor] <- "cor"
    add_significance <- function(pvalue) {
        case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**", pvalue < 0.05 ~ "*", TRUE ~ "ns")
    }
    dat_use$Signif <- add_significance(dat_use$pvalue)
    dat_use$Sign <- ifelse(dat_use$cor > 0, "Pos", "Neg")
    order_pos <- dat_use %>% dplyr::filter(cor > 0) %>% dplyr::arrange(Sign, cor) %>% dplyr::select(varible_others) %>% 
        unlist()
    order_neg <- dat_use %>% dplyr::filter(cor <= 0) %>% dplyr::arrange(Sign, desc(cor)) %>% dplyr::select(varible_others) %>% 
        unlist() %>% rev()
    dat_use$varible_others <- factor(dat_use$varible_others, levels = c(order_neg, order_pos))
    pos_limits <- max(dat_use$cor)
    neg_limits <- min(dat_use$cor)
    p <- ggplot(data = dat_use) + geom_segment(aes(x = 0, xend = cor, y = varible_others, yend = varible_others, 
        color = Sign), linewidth = 1, show.legend = F) + geom_point(aes(x = cor, y = varible_others, color = Sign), 
        size = 8) + geom_point(aes(x = cor, y = varible_others, color = Sign), size = 10, shape = 21, fill = NA) + 
        geom_text(data = dat_use %>% dplyr::filter(Sign == "Pos"), aes(x = -0.03, y = varible_others, label = varible_others, 
            colour = Sign), hjust = "right", size = 5, show.legend = F) + geom_text(data = dat_use %>% dplyr::filter(Sign == 
        "Neg"), aes(x = 0.03, y = varible_others, label = varible_others, colour = Sign), hjust = "left", size = 5, 
        show.legend = F) + geom_text(data = dat_use %>% dplyr::filter(Sign == "Pos"), aes(x = cor + 0.03, y = varible_others, 
        label = Signif, color = Sign), size = 8, vjust = 0.7, hjust = 0, show.legend = F) + geom_text(data = dat_use %>% 
        dplyr::filter(Sign == "Neg"), aes(x = cor - 0.03, y = varible_others, label = Signif, color = Sign), size = 8, 
        vjust = 0.7, hjust = 1, show.legend = F) + geom_vline(xintercept = 0, linewidth = 1) + scale_x_continuous(expand = expansion(mult = c(0.1, 
        0.1)), limits = c(neg_limits * 1.1, pos_limits * 1.1), breaks = c(round(neg_limits * 1.1, digits = 1), 
        round(neg_limits * 0.7, digits = 1), round(neg_limits * 0.3, digits = 1), 0, round(pos_limits * 0.3, digits = 1), 
        round(pos_limits * 0.7, digits = 1), round(pos_limits * 1.1, digits = 1))) + labs(x = "Correlation Coefficient", 
        y = "") + scale_color_manual(values = c(Pos = "#BC3C29FF", Neg = "#0072B5FF"), breaks = unique(dat_use$Sign)) + 
        theme_bw() + theme(axis.text.y = element_blank(), axis.text.x = element_text(color = "#000000", size = 18), 
        axis.ticks.y = element_blank(), axis.title.x = element_text(color = "#000000", size = 20), panel.border = element_rect(color = "#000000", 
            linewidth = 1.5), panel.grid.major = element_line(color = "white", linewidth = 0.7), panel.grid.minor = element_line(color = "white", 
            linewidth = 0.7), legend.text = element_text(color = "#000000", size = 12), legend.title = element_text(color = "#000000", 
            size = 15))
    print(p)
}

plot_expr_cor_Scatter_1 <- function(data, var1, var2, method = "pearson") {
    library(ggpubr)
    library(ggExtra)
    data$shape = ifelse(data[, var1] > median(data[, var1]) & data[, var2] > median(data[, var2]), "Mixed", ifelse(data[, 
        var1] > median(data[, var1]) & data[, var2] <= median(data[, var2]), var1, ifelse(data[, var1] <= median(data[, 
        var1]) & data[, var2] > median(data[, var2]), var2, "Quiescent")))
    data$shape = factor(data$shape, levels = levels(factor(data$shape)))
    shape = levels(factor(data$shape))
    length = length(levels(factor(data$shape)))
    bioCol = c("#BC3C29FF", "#0072B5FF", "#7876B1FF", "#6F99ADFF")
    p1 = ggplot(data, aes(data[, var1], data[, var2])) + xlab(var1) + ylab(var2) + geom_point(aes(colour = shape, 
        shape = shape)) + scale_color_manual(values = bioCol[1:length]) + scale_shape_manual(values = c(2, 8, 12, 
        18)) + geom_smooth(method = "lm", formula = y ~ x) + theme_bw() + stat_cor(method = method, aes(x = data[, 
        var1], y = data[, var2])) + theme(legend.position = "none")
    p2 = ggMarginal(p1, type = "densigram", xparams = list(fill = scales::alpha("#E18727FF", 0.5)), yparams = list(fill = scales::alpha("#20854EFF", 
        0.5)))
    print(p2)
}

plot_expr_diff_Bar <- function(data, Group = "Group", variable, levels = NULL, alphas = 0.5, mycolor = c("#BC3C29FF", 
    "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", "#E64B35FF")) {
    xmax = max(data[, variable])
    xmin = min(data[, variable])
    ymax = 3
    ymin = 0
    y = (ymax - ymin) * (data[, variable] - xmin)/(xmax - xmin) + ymin
    data[, variable] = y
    data[, variable] = unlist(as.numeric(data[, variable]))
    data[, "expression"] = data[, variable]
    colnames(data)[colnames(data) == Group] = "Group"
    if (is.null(levels)) {
        data[, "Group"] = factor(data[, "Group"], levels = unique(data[, "Group"]))
    }
    else {
        data[, "Group"] = factor(data[, "Group"], levels = levels)
    }
    if (length(unique(data[, "Group"])) == 2) {
        p_wilcox = wilcox.test(data[which(data$Group == unique(data[, "Group"])[1]), "expression"], data[which(data$Group == 
            unique(data[, "Group"])[2]), "expression"])$p.value
        p_wilcox = ifelse(p_wilcox < 0.001, "p<0.001", round(p_wilcox, 3))
        data_summary = data %>% group_by(Group) %>% summarise(mean_expression = mean(expression, na.rm = TRUE), 
            se_expression = sd(expression, na.rm = TRUE)/sqrt(n()), n = n())
        p_wilcox_df = data.frame(Group1 = unique(data[, "Group"])[1], Group2 = unique(data[, "Group"])[2], p_value = c(p_wilcox))
        p_wilcox_df$x_mid1 = c(1)
        p_wilcox_df$x_mid2 = c(2)
        p_wilcox_df$y_value = c(max(data_summary$mean_expression) * 1.05)
        histogram = ggplot(data_summary, aes(x = Group, y = mean_expression, fill = Group)) + geom_col(position = "dodge", 
            width = 0.4) + geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), 
            width = 0.2, position = position_dodge(width = 0.9)) + ylim(0, max(data_summary$mean_expression) * 
            1.1) + xlab("") + ylab(paste0(variable, " Mean Expression")) + labs(fill = "Group") + theme_minimal() + 
            theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), 
                axis.text = element_text(size = 9.5, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, 
                  vjust = 1)) + scale_color_manual(values = scales::alpha(mycolor, alphas)) + scale_fill_manual(values = scales::alpha(mycolor, 
            alphas)) + geom_segment(data = p_wilcox_df, aes(x = x_mid1, xend = x_mid2, y = y_value, yend = y_value), 
            colour = "black", size = 0.5, lineGroup = "dashed", inherit.aes = FALSE) + geom_text(data = p_wilcox_df, 
            aes(x = (x_mid1 + x_mid2)/2, y = y_value, label = p_value), vjust = -0.5, size = 3.5, check_overlap = TRUE, 
            inherit.aes = FALSE)
        print(histogram)
    }
    else if (length(unique(data[, "Group"])) == 3) {
        p_wilcox_1v2 = wilcox.test(data[which(data$Group == unique(data[, "Group"])[1]), "expression"], data[which(data$Group == 
            unique(data[, "Group"])[2]), "expression"])$p.value
        p_wilcox_1v2 = ifelse(p_wilcox_1v2 < 0.001, "p<0.001", round(p_wilcox_1v2, 3))
        p_wilcox_1v3 = wilcox.test(data[which(data$Group == unique(data[, "Group"])[1]), "expression"], data[which(data$Group == 
            unique(data[, "Group"])[3]), "expression"])$p.value
        p_wilcox_1v3 = ifelse(p_wilcox_1v3 < 0.001, "p<0.001", round(p_wilcox_1v3, 3))
        p_wilcox_2v3 = wilcox.test(data[which(data$Group == unique(data[, "Group"])[2]), "expression"], data[which(data$Group == 
            unique(data[, "Group"])[3]), "expression"])$p.value
        p_wilcox_2v3 = ifelse(p_wilcox_2v3 < 0.001, "p<0.001", round(p_wilcox_2v3, 3))
        data_summary = data %>% group_by(Group) %>% summarise(mean_expression = mean(expression, na.rm = TRUE), 
            se_expression = sd(expression, na.rm = TRUE)/sqrt(n()), n = n())
        p_wilcox_df = data.frame(Group1 = c(unique(data[, "Group"])[1], unique(data[, "Group"])[1], unique(data[, 
            "Group"])[2]), Group2 = c(unique(data[, "Group"])[2], unique(data[, "Group"])[3], unique(data[, "Group"])[3]), 
            p_value = c(p_wilcox_1v2, p_wilcox_1v3, p_wilcox_2v3))
        p_wilcox_df$x_mid1 = c(1, 1, 2)
        p_wilcox_df$x_mid2 = c(2, 3, 3)
        p_wilcox_df$y_value = c(max(data_summary$mean_expression) * 1.05, max(data_summary$mean_expression) * 1.2, 
            max(data_summary$mean_expression) * 1.35)
        histogram = ggplot(data_summary, aes(x = Group, y = mean_expression, fill = Group)) + geom_col(position = "dodge", 
            width = 0.4) + geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), 
            width = 0.2, position = position_dodge(width = 0.9)) + ylim(0, max(data_summary$mean_expression) * 
            1.35) + xlab("") + ylab(paste0(variable, " Mean Expression")) + labs(fill = "Group") + theme_minimal() + 
            theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), 
                axis.text = element_text(size = 9.5, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, 
                  vjust = 1)) + scale_color_manual(values = scales::alpha(mycolor, alphas)) + scale_fill_manual(values = scales::alpha(mycolor, 
            alphas)) + geom_segment(data = p_wilcox_df, aes(x = x_mid1, xend = x_mid2, y = y_value, yend = y_value), 
            colour = "black", size = 0.5, lineGroup = "dashed", inherit.aes = FALSE) + geom_text(data = p_wilcox_df, 
            aes(x = (x_mid1 + x_mid2)/2, y = y_value, label = p_value), vjust = -0.5, size = 3.5, check_overlap = TRUE, 
            inherit.aes = FALSE)
        print(histogram)
    }
    else {
        p_kruskal = kruskal.test(expression ~ Group, data = data)$p.value
        p_kruskal = ifelse(p_kruskal < 0.001, "p<0.001", round(p_kruskal, 3))
        data_summary = data %>% group_by(Group) %>% summarise(mean_expression = mean(expression, na.rm = TRUE), 
            se_expression = sd(expression, na.rm = TRUE)/sqrt(n()), n = n())
        histogram = ggplot(data_summary, aes(x = Group, y = mean_expression, fill = Group)) + geom_col(position = "dodge", 
            width = 0.4) + geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), 
            width = 0.2, position = position_dodge(width = 0.9)) + ylim(0, max(data_summary$mean_expression) * 
            1.1) + xlab("") + ylab(paste0(variable, " Mean Expression")) + labs(fill = "Group") + theme_minimal() + 
            theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), 
                axis.text = element_text(size = 9.5, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, 
                  vjust = 1)) + scale_color_manual(values = scales::alpha(mycolor, alphas)) + scale_fill_manual(values = scales::alpha(mycolor, 
            alphas)) + geom_text(aes(x = (1 + length(unique(data[, "Group"])))/2, y = max(data_summary$mean_expression) * 
            1.05, label = p_kruskal), vjust = -0.5, size = 3.5, check_overlap = TRUE, inherit.aes = FALSE)
        print(histogram)
    }
}

plot_expr_diff_Ridge <- function(data, Group = "Group", variable, levels = NULL, test.methods = "wilcox.test", 
    alphas = 0.5, mycolor = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", 
        "#EE4C97FF", "#E64B35FF")) {
    suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))
    suppressPackageStartupMessages(library(ggExtra, quietly = TRUE))
    suppressPackageStartupMessages(library(aplot, quietly = TRUE))
    data[, variable] = unlist(as.numeric(data[, variable]))
    data[, "expression"] = data[, variable]
    colnames(data)[colnames(data) == Group] = "Group"
    if (test.methods == "wilcox.test") {
        p.val = wilcox.test(expression ~ Group, data = data)
    }
    else if (test.methods == "kruskal.test") {
        p.val = kruskal.test(expression ~ Group, data = data)
    }
    else if (test.methods == "t.test") {
        p.val = t.test(expression ~ Group, data = data)
    }
    p.lab = paste0("P", ifelse(p.val$p.value < 0.001, " < 0.001", paste0(" = ", round(p.val$p.value, 3))))
    if (is.null(levels)) {
        data[, "Group"] = factor(data[, "Group"], levels = unique(data[, "Group"]))
    }
    else {
        data[, "Group"] = factor(data[, "Group"], levels = levels)
    }
    p_top = ggplot(data, aes(x = expression, color = Group, fill = Group)) + geom_density() + scale_color_manual(values = alpha(mycolor, 
        alphas)) + scale_fill_manual(values = alpha(mycolor, alphas)) + theme_classic() + xlab(paste0("Estimated Expression of ", 
        variable)) + ylab(NULL) + theme(legend.position = "none", legend.title = element_blank(), axis.text.x = element_text(size = 12, 
        color = "black"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_rug()
    p_bot = ggplot(data, aes(Group, expression, fill = Group)) + geom_boxplot(aes(col = Group)) + scale_color_manual(values = alpha(mycolor, 
        alphas)) + scale_fill_manual(values = alpha(mycolor, alphas)) + xlab(NULL) + ylab("Estimated Expression") + 
        theme_void() + theme(legend.position = "none", legend.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 11, color = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + annotate(geom = "text", x = 1.25, hjust = 1, y = max(data[, "expression"]) * 
        1, size = 4, angle = 270, fontface = "bold", label = p.lab) + coord_flip()
    dat = ggplot_build(p_bot)$data[[1]]
    p_bot = p_bot + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), color = "white", 
        inherit.aes = F)
    p_bot
    p <- p_top %>% insert_bottom(p_bot, height = 0.4)
    print(p)
}

plot_expr_diff_Violin <- function(data, Group = "Group", variable = "NINJ1", colors = c(Normal = "#00599F", Tumor = "#D01910"), 
    method = "t.test") {
    library(ggplot2)
    library(ggpubr)
    if (!is.data.frame(data)) 
        stop("data必须为数据框")
    if (!Group %in% names(data)) 
        stop(paste(Group, "不存在于数据框中"))
    if (!variable %in% names(data)) 
        stop(paste(variable, "不存在于数据框中"))
    data[[Group]] <- factor(data[[Group]])
    group_levels <- levels(data[[Group]])
    if (length(group_levels) != 2) 
        stop("分组变量必须包含两个水平")
    if (is.null(names(colors))) {
        if (length(colors) == 2) {
            names(colors) <- group_levels
        }
        else {
            stop("颜色参数需为命名向量或长度为2的未命名向量")
        }
    }
    else {
        if (!all(group_levels %in% names(colors))) {
            stop("颜色名称需与分组水平完全匹配: ", paste(group_levels, collapse = ", "))
        }
        colors <- colors[group_levels]
    }
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[Group]], y = .data[[variable]], fill = .data[[Group]])) + 
        ggplot2::geom_violin(trim = FALSE, alpha = 0.8, width = 0.7, color = NA) + ggplot2::geom_boxplot(width = 0.2, 
        color = "black", alpha = 0, outlier.shape = NA) + ggplot2::geom_jitter(size = 2.5, shape = 21, color = "black", 
        position = ggplot2::position_jitter(width = 0.1, height = 0), alpha = 0.8) + ggplot2::scale_fill_manual(values = colors) + 
        ggpubr::stat_compare_means(method = method, comparisons = list(group_levels), label = "p.signif", size = 5, 
            tip.length = 0) + ggplot2::labs(y = paste0(variable, " Expression Level"), x = "") + ggplot2::theme_classic(base_size = 12) + 
        ggplot2::theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), 
            axis.text = ggplot2::element_text(color = "black", size = 11), axis.title.y = ggplot2::element_text(size = 12, 
                margin = ggplot2::margin(r = 10)), axis.line = ggplot2::element_line(linewidth = 0.5), plot.margin = ggplot2::unit(c(5, 
                15, 5, 5), "pt"))
    print(p)
}

plot_expr_multidiff_Box <- function(data, variable, Group, Expression, color_use = c("#00599F", "#D01910")) {
    data$expr = data[, Expression]
    data$Type = data[, Group]
    data$variable = data[, variable]
    pvalues = sapply(unique(data$variable), function(x) {
        res = wilcox.test(as.numeric(expr) ~ Type, data = subset(data, variable == x))
        res$p.value
    })
    pv = data.frame(cancer_type = unique(data$variable), pvalue = pvalues)
    pv$sigcode = cut(pv$pvalue, c(0, 1e-04, 0.001, 0.01, 0.05, 1), labels = c("****", "***", "**", "*", ""))
    count_N = data %>% group_by(variable, Type) %>% tally
    count_N$n = paste("n =", count_N$n)
    ylabname = Expression
    boxplot_res = ggplot(data, aes(x = variable, y = expr, color = Type, fill = Type)) + geom_point(shape = 21, 
        size = 2, position = position_jitterdodge(), alpha = 0.35) + geom_boxplot(notch = F, outlier.size = -1, 
        color = "black", lwd = 0.8, alpha = 0.7) + theme_classic() + scale_color_manual(values = color_use) + scale_fill_manual(values = color_use) + 
        theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + theme(axis.text.x = element_text(colour = "black", 
        size = 11, angle = 45, hjust = 0.5, vjust = 0.5), legend.position = "bottom") + geom_text(aes(x = cancer_type, 
        y = max(data[, Expression]), label = pv$sigcode), data = pv, inherit.aes = F) + xlab("") + ylab("") + geom_text(data = count_N, 
        aes(label = n, y = min(data[, Expression]) * 1.3, color = Type), position = position_dodge2(0.9), size = 3, 
        angle = 90, hjust = 0) + ylim(min(data[, Expression]) * 1.4, max(data[, Expression])) + theme(axis.text.x = element_text(angle = 45, 
        hjust = 1, vjust = 1), legend.position = "bottom", panel.border = element_rect(colour = "black", fill = NA, 
        size = 0.6))
    res_diff_pancancer = list(original_data = data, pvalue_res = pv, boxplot_res = boxplot_res)
    res_diff_pancancer
}

plot_heatmap_1 <- function(hplot_input, bar_color_use = "steelblue", legend_heatmap_use, legend_topann_use, cellwidth_use = 1, 
    cellheight_use = 0.5, fontsize_use = 10, n = 50, topann_color_use = c("#BC3C29FF", "#0072B5FF", "#E18727FF", 
        "#20854EFF", "#7876B1FF"), heatmap_name_use = "heatmap_name_use") {
    library(tidyverse)
    library(ComplexHeatmap)
    hplot <- function(hplot_input, right_annotation_use, topann_color_use, bar_color_use, cellwidth, cellheight, 
        legend_heatmap_use, legend_topann_use, fontsize_use) {
        col_ha <- columnAnnotation(legend_topann_use = colnames(hplot_input), col = list(legend_topann_use = topann_color_use), 
            show_annotation_name = F)
        row_ha <- rowAnnotation(bar = anno_barplot(right_annotation_use, bar_width = 0.8, border = FALSE, gp = gpar(fill = bar_color_use, 
            col = NA), add_numbers = TRUE, numbers_offset = unit(-10, "mm"), axis_param = list(labels_rot = 0), 
            numbers_gp = gpar(fontsize = fontsize_use, col = "white"), width = unit(3, "cm")), show_annotation_name = FALSE)
        Heatmap(as.matrix(hplot_input), name = legend_heatmap_use, right_annotation = row_ha, top_annotation = col_ha, 
            col = c("#27879f", "#e7efc7", "#ef7740"), rect_gp = gpar(col = "black", lwd = 1), cluster_columns = FALSE, 
            cluster_rows = FALSE, show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left", width = unit(cellwidth * 
                ncol(hplot_input) + 2, "cm"), height = unit(cellheight * nrow(hplot_input), "cm"), column_split = factor(colnames(hplot_input), 
                levels = colnames(hplot_input)), column_title = NULL, cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(label = format(hplot_input[i, j], digits = 2, nsmall = 2), x, y, gp = gpar(fontsize = 10))
            })
    }
    for (i in 1:ncol(hplot_input)) {
        hplot_input[, i] <- as.numeric(hplot_input[, i])
    }
    df <- hplot_input %>% mutate(mean = apply(hplot_input, 1, mean)) %>% arrange(., desc(mean)) %>% dplyr::slice_head(n = n)
    right_annotation_use <- round(df$mean, 2)
    hplot_input <- df[, -ncol(df)]
    topann_color_use <- setNames(topann_color_use[1:length(colnames(hplot_input))], colnames(hplot_input))
    hm <- hplot(hplot_input, right_annotation_use, topann_color_use, bar_color_use, cellwidth_use, cellheight_use, 
        legend_heatmap_use, legend_topann_use, fontsize_use)
    existing_color_mapping <- hm@top_annotation@anno_list[["legend_topann_use"]]@color_mapping
    colors <- existing_color_mapping@colors
    new_color_mapping <- ColorMapping(levels = colnames(hplot_input), colors = colors[colnames(hplot_input)], name = legend_topann_use)
    hm@top_annotation@anno_list[["legend_topann_use"]]@color_mapping <- new_color_mapping
    invisible(NULL)
    draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
    invisible(invisible(NULL))
}

plot_sc_CellType_Cortplot <- function(sc, idents, version = "v4") {
    library(RColorBrewer)
    library(corrplot)
    library(Seurat)
    if (version == "v5") {
        sc <- JoinLayers(sc)
        sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")
    }
    else {
        sc[["RNA"]] <- as(sc[["RNA"]], Class = "Assay")
    }
    Idents(sc) <- idents
    Markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.585, verbose = TRUE, max.cells.per.ident = 2000)
    Markers <- subset(Markers[grep("^RP[L|S]", Markers$gene, ignore.case = FALSE, invert = TRUE), ], subset = p_val_adj < 
        0.05)
    Markers_av <- AverageExpression(sc, group.by = idents, features = unique(Markers$gene), assays = "RNA")
    Markers_av <- Markers_av$RNA
    gene_cell_exp <- t(scale(t(Markers_av), scale = TRUE, center = TRUE))
    cell_cor <- cor(as.matrix(gene_cell_exp), method = "spearman")
    corrplot(cell_cor, method = "circle", order = "hclust", hclust.method = "ward.D", type = "full", addrect = length(unique(Idents(sc))), 
        col = rev(brewer.pal(n = 8, name = "RdYlBu")), tl.cex = 0.8, tl.col = "black", cl.cex = 0.5, cl.pos = "r", 
        cl.length = 5, cl.ratio = 0.1, mar = c(3, 0, 2, 0), outline = TRUE, addgrid.col = "lightgray")
}

plot_sc_Markers_Simplified <- function(scRNA, label = "seurat_clusters", Immune_cells = c("PTPRC", "CD53", "CD3E", 
    "CD19", "CD14", "CD68"), Cardiacs = c("TTN", "MYH7", "MYH6", "TNNT2"), Epithelials = c("EPCAM", "KRT19", "CDH1", 
    "KRT18"), Fibroblasts = c("DCN", "COL1A1", "COL1A2", "THY1"), Endothelials = c("PECAM1", "CDH5", "VWF", "TEK"), 
    SMCs = c("ACTA2", "TAGLN", "PDGFRB", "MYH11"), Mast_cells = c("GATA2", "KIT", "MS4A2", "TPSB2"), Oligodendrocytes = c("OLIG1", 
        "OLIG2", "MOG", "CLDN11"), Pericyte = c("KCNJ8", "RGS5"), Plasmas = c("IGHG1", "MZB1", "TNFRSF17", "XBP1"), 
    Cyclings = c("STMN1", "MKI67", "TOP2A", "CDK1"), mycolor = c("#BC3C29FF", "#208A42", "#0072B5FF", "#E18727FF", 
        "#20854EFF", "#7876B1FF", "#6F99ADFF", "#4DBBD5FF", "#00A087FF", "#6A6599FF", "#79AF97FF")) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(tidyverse)
        library(cowplot)
    })
    categories <- c("Immune_cells", "Cardiacs", "Epithelials", "Fibroblasts", "Endothelials", "SMCs", "Mast_cells", 
        "Oligodendrocytes", "Pericyte", "Plasmas", "Cyclings")
    selected_markers <- list()
    selected_labels <- list()
    for (category in categories) {
        markers <- get(category)
        available_markers <- intersect(rownames(scRNA), markers)
        if (length(available_markers) > 0) {
            selected_markers[[category]] <- available_markers
            selected_labels[[category]] <- rep(category, length(available_markers))
        }
    }
    all_markers <- unlist(selected_markers)
    all_labels <- unlist(selected_labels)
    data.usage <- DotPlot(scRNA, features = all_markers, group.by = label)$data
    data.anno <- data.frame(features.plot = all_markers, label = all_labels)
    df.plot <- plyr::join(data.usage, data.anno)
    p <- ggplot(df.plot, aes(x = features.plot, y = as.numeric(id), size = pct.exp, color = avg.exp.scaled)) + 
        geom_point() + scale_size("% detected", range = c(0, 10)) + scale_colour_gradient2("Average\nexpression", 
        low = c("#2F5B89", "#6B9AB7", "lightblue"), mid = "white", high = c("#F39B7FFF", "#ED5E57", "#A80C3A"), 
        midpoint = 0) + cowplot::theme_cowplot() + ylab("") + xlab("") + theme_bw() + scale_y_continuous(breaks = 1:length(levels(df.plot$id)), 
        labels = levels(df.plot$id), sec.axis = dup_axis()) + facet_grid(~label, scales = "free_x", space = "free") + 
        theme_classic() + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5, color = "black", 
        face = "bold"), axis.text.y = element_text(size = 10, color = "skyblue", face = "bold"), axis.title.x = element_text(size = 14, 
        colour = "black", hjust = 0.5), axis.ticks.y = element_blank(), axis.text.y.right = element_blank(), axis.ticks.x = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.2), panel.spacing = unit(0, "mm"), strip.text.x = element_text(size = 10, 
            face = "bold", color = "#FFFFFF", vjust = 0.5, margin = margin(b = 3, t = 3)), strip.background = element_rect(colour = "black", 
            fill = "black", size = 1))
    g <- ggplot_gtable(ggplot_build(p))
    strips <- which(grepl("strip-", g$layout$name))
    cols <- mycolor[1:length(unique(df.plot$label))] %>% as.vector()
    for (i in seq_along(strips)) {
        k <- which(grepl("rect", g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
        l <- which(grepl("titleGrob", g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
        g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
    }
    plot(g)
}

plot_sc_UMAP_Feature_3 <- function(scRNA, marker, slot_use = "counts", reduction_use = "umap", size = 0.8) {
    library(Seurat)
    library(ggplot2)
    library(Nebulosa)
    library(ggnetwork)
    color = c(alpha("lightgrey", 0.1), "#330066", "#336699", "#66CC66", "#FFCC33", "red")
    plot = plot_density(scRNA, features = marker, slot = slot_use, reduction = reduction_use, pal = "magma", raster = TRUE, 
        size = size) + theme_blank() + theme(panel.background = element_rect(fill = "black")) & theme(legend.frame = element_rect(colour = "black"), 
        legend.ticks = element_line(colour = "black", linewidth = 0), legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.8, 
            "cm"), legend.title = element_text(color = "black", face = "bold", size = 8))
    print(plot)
    plot2 = plot_density(scRNA, features = marker, slot = slot_use, reduction = reduction_use, pal = "magma", raster = TRUE, 
        size = size) + theme_blank() + theme(legend.frame = element_rect(colour = "black"), legend.ticks = element_line(colour = "black", 
        linewidth = 0), legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.8, "cm"), legend.title = element_text(color = "black", 
        face = "bold", size = 8))
    print(plot2)
}

plot_sc_UMAP_cell <- function(scRNA, plotname = "CellType", reduction = "umap", label = "celltype", alphas = 0.8, 
    mycolor = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", 
        "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", 
        "#7E6148FF", "#B09C85FF", "#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", "#008280FF", "#BB0021FF", 
        "#5F559BFF", "#A20056FF", "#808180FF", "#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", 
        "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", 
        "#6A6599FF", "#80796BFF", "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2", 
        "#b5bd61", "#17becf", "#aec7e8")) {
    library(ggplot2)
    library(tidyverse)
    library(ggrepel)
    umap = scRNA@reductions[[reduction]]@cell.embeddings %>% as.data.frame() %>% cbind(Cluster = scRNA@meta.data[, 
        label])
    colnames(umap) = c("UMAP_1", "UMAP_2", "Cluster")
    cell_type_med = umap %>% dplyr::group_by(Cluster) %>% dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    UMAP_1_min = min(umap$UMAP_1)
    UMAP_2_min = min(umap$UMAP_2)
    p = ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) + geom_point(size = 0.1, alpha = alphas) + scale_color_manual(values = mycolor) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
            axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), 
            plot.background = element_rect(fill = "white")) + theme(legend.title = element_blank(), legend.key = element_rect(fill = "white"), 
        legend.text = element_text(size = 20), legend.key.size = unit(1, "cm")) + guides(color = guide_legend(override.aes = list(size = 5))) + 
        annotate("segment", x = UMAP_1_min, y = UMAP_2_min, xend = UMAP_1_min + 3, yend = UMAP_2_min, colour = "black", 
            size = 1, arrow = arrow(length = unit(0.3, "cm"))) + annotate("segment", x = UMAP_1_min, y = UMAP_2_min, 
        xend = UMAP_1_min, yend = UMAP_2_min + 3, colour = "black", size = 1, arrow = arrow(length = unit(0.3, 
            "cm"))) + annotate("text", x = UMAP_1_min + 1.5, y = UMAP_2_min - 1, label = "UMAP_1", color = "black", 
        size = 3, fontface = "bold") + annotate("text", x = UMAP_1_min - 1, y = UMAP_2_min + 1.5, label = "UMAP_2", 
        color = "black", size = 3, fontface = "bold", angle = 90) + geom_label_repel(aes(label = Cluster), fontface = "bold", 
        data = cell_type_med, point.padding = unit(0.5, "lines")) + theme(legend.position = "none")
    print(p)
}

plot_survival_KM_1 <- function(data, minprop = 0.3, GeneName, time, event, Timeunit = "year") {
    suppressPackageStartupMessages(library(survival))
    suppressPackageStartupMessages(library(survminer))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(ggplot2))
    data <- data %>% dplyr::select(c(time, event, GeneName))
    if (Timeunit == "month") {
        t <- data %>% dplyr::pull(time) %>% as.numeric()
        t <- t/30.4
        data[, 1] <- t
    }
    else if (Timeunit == "day") {
        t <- data %>% dplyr::pull(time) %>% as.numeric()
        t <- t/365
        data[, 1] <- t
    }
    else if (Timeunit == "y") {
        Timeunit <- "year"
    }
    bestcut <- surv_cutpoint(data, time = colnames(data)[1], event = colnames(data)[2], variables = GeneName, minprop = minprop)
    cutoff <- bestcut$cutpoint[1, 1]
    data$group <- factor(ifelse(data[, GeneName] > cutoff, "High", "Low"), levels = c("Low", "High"))
    data$group2 <- factor(ifelse(data[, GeneName] > median(data[, GeneName]), "High", "Low"), levels = c("Low", 
        "High"))
    SurvivalType <- ifelse(colnames(data)[2] == "OS", "Overall Survival", ifelse(colnames(data)[2] == "DSS", "Disease-specific Survival", 
        ifelse(colnames(data)[2] == "PFS", "Progression-free Survival", ifelse(colnames(data)[2] == "DFS", "Disease-free Survival", 
            ifelse(colnames(data)[2] == "RFS", "Relapse-free Survival", ifelse(colnames(data)[2] == "DFI", "Disease-free Interval", 
                ifelse(colnames(data)[2] == "PFI", "Progression-free Interval", "Survival")))))))
    time_label <- switch(Timeunit, year = "years", month = "months", day = "days", Timeunit)
    survdiff_formula1 <- as.formula(paste0("Surv(", colnames(data)[1], ",", colnames(data)[2], ") ~ group"))
    diff1 <- survdiff(survdiff_formula1, data = data, na.action = na.exclude)
    pValue1 <- 1 - pchisq(diff1$chisq, length(diff1$n) - 1)
    fit1 <- survfit(survdiff_formula1, data = data)
    fit1$call$formula <- survdiff_formula1
    legend_labs1 <- c(paste0("Low ", GeneName, " (n=", sum(data$group == "Low", na.rm = TRUE), ")"), paste0("High ", 
        GeneName, " (n=", sum(data$group == "High", na.rm = TRUE), ")"))
    surv_plot1 <- ggsurvplot(fit = fit1, data = data, palette = c("#00599F", "#D01910"), linetype = 1, size = 1.2, 
        conf.int = TRUE, conf.int.alpha = 0.15, risk.table = FALSE, pval = TRUE, pval.coord = c(0, 0.05), pval.size = 5, 
        xlab = paste0("Time (", time_label, ")"), ylab = paste0(SurvivalType, " Probability"), legend.title = GeneName, 
        legend.labs = legend_labs1, legend.position = c(0.85, 0.85), ggtheme = theme_classic() + theme(text = element_text(size = 14), 
            axis.line = element_line(size = 0.8), axis.ticks = element_line(size = 0.8), axis.title = element_text(face = "bold"), 
            legend.title = element_text(face = "bold"), legend.background = element_rect(fill = "white", color = "black", 
                size = 0.3), panel.grid.major.y = element_line(color = "white", size = 0.3)))
    survdiff_formula2 <- as.formula(paste0("Surv(", colnames(data)[1], ",", colnames(data)[2], ") ~ group2"))
    diff2 <- survdiff(survdiff_formula2, data = data, na.action = na.exclude)
    pValue2 <- 1 - pchisq(diff2$chisq, length(diff2$n) - 1)
    fit2 <- survfit(survdiff_formula2, data = data)
    fit2$call$formula <- survdiff_formula2
    legend_labs2 <- c(paste0("Low ", GeneName, " (n=", sum(data$group2 == "Low", na.rm = TRUE), ")"), paste0("High ", 
        GeneName, " (n=", sum(data$group2 == "High", na.rm = TRUE), ")"))
    surv_plot2 <- ggsurvplot(fit = fit2, data = data, palette = c("#00599F", "#D01910"), linetype = 1, size = 1.2, 
        conf.int = TRUE, conf.int.alpha = 0.15, risk.table = FALSE, pval = TRUE, pval.coord = c(0, 0.05), pval.size = 5, 
        xlab = paste0("Time (", time_label, ")"), ylab = paste0(SurvivalType, " Probability"), legend.title = GeneName, 
        legend.labs = legend_labs2, legend.position = c(0.85, 0.85), ggtheme = theme_classic() + theme(text = element_text(size = 14), 
            axis.line = element_line(size = 0.8), axis.ticks = element_line(size = 0.8), axis.title = element_text(face = "bold"), 
            legend.title = element_text(face = "bold"), legend.background = element_rect(fill = "white", color = "black", 
                size = 0.3), panel.grid.major.y = element_line(color = "white", size = 0.3)))
    return(list(best_cutoff_plot = surv_plot1, median_cutoff_plot = surv_plot2, best_cutoff = cutoff))
}

plot_venn_1 <- function(data.list, set.names = NULL, fill.colors = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF"), 
    fill.alpha = 0.5, border.color = "black", border.alpha = 1, border.width = 1, border.linetype = 1, set.name.colors = "black", 
    set.name.size = 6, text.color = "black", text.size = 5.5) {
    library(ggplot2)
    library(ggvenn)
    library(dplyr)
    library(tidyr)
    library(purrr)
    if (length(data.list) <= 1) {
        stop("The number of elements in data.list must be greater than 1!")
    }
    if (length(data.list) > 4) {
        stop("This function currently only supports venn plot of 2-4 sets!")
    }
    if (is.null(set.names)) {
        if (is.null(names(data.list))) {
            stop("Elements in data.list must have names, or set.names must be prepared!")
        }
        else {
            set.names <- names(data.list)
        }
    }
    data.length <- length(set.names)
    data <- purrr::map2_df(data.list, names(data.list), ~data.frame(value = .x, Class = .y, Type = TRUE)) %>% tidyr::pivot_wider(names_from = "Class", 
        values_from = "Type") %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(., FALSE)))
    venn_plot <- geom_venn(set_names = set.names, fill_color = fill.colors, fill_alpha = fill.alpha, stroke_color = border.color, 
        stroke_alpha = border.alpha, stroke_linetype = border.linetype, stroke_size = border.width, set_name_color = set.name.colors, 
        set_name_size = set.name.size, text_color = text.color, text_size = text.size, show_percentage = FALSE, 
        digits = 1, show_outside = "none")
    if (data.length == 2) {
        colnames(data)[-1] <- LETTERS[1:2]
        p <- ggplot(data, aes(A = A, B = B)) + theme_void() + coord_fixed() + venn_plot
    }
    if (data.length == 3) {
        colnames(data)[-1] <- LETTERS[1:3]
        p <- ggplot(data, aes(A = A, B = B, C = C)) + theme_void() + coord_fixed() + venn_plot
    }
    if (data.length == 4) {
        colnames(data)[-1] <- LETTERS[1:4]
        p <- ggplot(data, aes(A = A, B = B, C = C, D = D)) + theme_void() + coord_fixed() + venn_plot
    }
    print(p)
}



# Additional in-memory utility: score_model
# Data-frame interface; no input/output paths or implicit normalization.
score_model <- function(expression, coefficients) {
    if (!all(c("gene", "coefficient") %in% names(coefficients)))
        stop("Coefficients require gene and coefficient columns")
    if (anyDuplicated(coefficients$gene) || anyNA(coefficients$gene))
        stop("Gene symbols must be present and unique")
    if (!all(coefficients$gene %in% names(expression))) stop("Missing model genes")
    X <- as.matrix(expression[, coefficients$gene, drop=FALSE])
    if (!is.numeric(X) || any(!is.finite(X)) || !is.numeric(coefficients$coefficient) || any(!is.finite(coefficients$coefficient)))
        stop("Inputs must be finite numeric values")
    eta <- drop(X %*% coefficients$coefficient)
    RS <- exp(eta)
    if (any(!is.finite(RS))) stop("Score overflow; verify input scale")
    result <- data.frame(linear_predictor=eta, RS=RS)
    if ("ID" %in% names(expression)) result <- cbind(ID=expression$ID, result)
    result
}


# Additional in-memory utility: audit_models
# Input is an in-memory table of stored values, not external validation.
audit_models <- function(scores) {
    stopifnot(all(c("ID", "Model", "Cindex") %in% names(scores)),
              !anyDuplicated(scores[c("ID", "Model")]), all(is.finite(scores$Cindex)))
    ranked <- aggregate(Cindex ~ Model, scores, mean)
    ranked <- ranked[order(-ranked$Cindex), ]
    ranked$rank <- seq_len(nrow(ranked))
    ranked
}
