# Table A5. Early-career controls
library(arrow)
library(data.table)
library(fixest)

data_dir <- "model-data"  # Change this to the folder containing the extracted data.
setDTthreads(8)
arrow::set_cpu_count(8)

files <- list.files(file.path(data_dir, "main"), "[.]parquet$", full.names = TRUE)
papers <- rbindlist(lapply(files, read_parquet))
stopifnot(nrow(papers) == 100304404L)
papers[, `:=`(author_id = as.numeric(author_id), work_id = as.numeric(work_id))]
early <- as.data.table(read_parquet(file.path(data_dir, "early-background.parquet")))
early_row <- match(papers$author_id, as.numeric(early$author_id))
papers <- papers[!is.na(early_row)]
early_row <- early_row[!is.na(early_row)]
stopifnot(nrow(papers) == 84205769L)

# All models use paper weights and cluster by author and paper.
adjustment <- ssc(K.adj = FALSE, G.adj = FALSE, t.df = "min")
cluster <- vcov_cluster(~ author_id + work_id, ssc = adjustment, vcov_fix = FALSE)

estimate <- function(formula, data) {
    model <- feols(formula, data = data, weights = ~ paper_weight, vcov = cluster,
                   ssc = adjustment, fixef.rm = "none", fixef.tol = 1e-9,
                   fixef.iter = 10000L, mem.clean = TRUE, nthreads = 8L, lean = TRUE)
    stopifnot(nobs(model) == nrow(data), all(is.finite(se(model))))
    model
}

# Fit the baseline and add each interaction separately on the same sample.
baseline <- outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2 | author_id + paper_year
models <- list(Base = estimate(baseline, papers))
controls <- c("log_coauthors", "log_papers", "log_training_h", "log_training_citations", "early_idr")
for (control in controls) {
    papers[, early_interaction := idr_c * early[[control]][early_row]]
    models[[control]] <- estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2 + early_interaction
                                 | author_id + paper_year, papers)
}
print(models)
