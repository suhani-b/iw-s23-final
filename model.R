# Author: Suhani Balachandran
# Date: Spring 2023
# Project: BSE IW - single cell lm
# Advisors: Mona Singh, Sara Geraghty

# Filename: model.py
# Purpose: implement various linear modeling functions on the data
# Parameters: path to mutations matrix (csv), path to GE matrix (csv)
# Input: SNV matrix (samples x driver genes), GE matrix (samples x all genes)
# Output: file of most significant combinations for qvalue correction

library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
snv_path <- args[1] # ex: /Genomics/grid/users/sb5334/single_cell_lm/
ge_path <- args[2]

# read and format csv data
snv <- read.csv(csv_path)
ge <- t(read.csv(ge_path))

# make df to hold results
lm_results <- data.frame(matrix(ncol = 6, nrow = 0))
ridge_results <- 
lasso_results <-

# FUNCTION FOR REGULAR LM
regular_lm <- function(driver, target) {
    # set up df needed for lm: cells by [cnv of driver, ge of target gene]
    data <- as.data.frame(cbind(snv[, driver], ge[, target]))
    colnames(data) <- c("SNV_Driver", "GE_Target")

    # perform lm
    matrix_coef <- as.data.frame(
        summary(lm(GE_Target ~ CNV_Driver, data))$coefficients)
    if (!is.nan(matrix_coef$"Pr(>|t|)"[2])) {
        matrix_coef <- cbind(matrix_coef[2, ],
                    c(drivers_list[driver]), c(targets_list[target]))
        return(matrix_coef[1,])
    } else {
        print("p value NaN")
    }
}

# lm loop
drivers_list <- colnames(snv)
targets_list <- colnames(ge)

for (driver in seq_len(ncol(snv))) {
    for (target in seq_len(ncol(ge))) {
        lm_results <- rbind(lm_results, regular_lm(driver, target))
gi
    }
}

# reformat results and export
results <- setNames(results,
        c("beta", "std error", "t statistic", "p value", "driver", "target"))
rownames(results) <- NULL
write.csv(results, paste0(path_prefix, "model/output/", title, "LM.csv"))