# Table A7. Fixed five- and ten-year citation windows
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

# Each window uses its own outcome labels and eligible authors.
windows <- c(`5` = 2020L, `10` = 2015L)
expected <- c(`5` = 96284291L, `10` = 68453652L)
columns <- setdiff(names(papers), "subfield_id")
models <- list()
for (window in names(windows)) {
    prefix <- file.path(data_dir, paste0("window-", window))
    labels <- read_parquet(paste0(prefix, ".parquet"))
    authors <- read_parquet(paste0(prefix, "-authors.parquet"))
    label_row <- match(papers$work_id, as.numeric(labels$work_id))
    keep <- papers$paper_year <= windows[[window]] & !is.na(label_row) &
        papers$paper_year == labels$paper_year[label_row] &
        papers$author_id %in% as.numeric(authors$author_id)
    sample <- papers[keep, ..columns]
    sample[, outcome := labels$fixed_outcome[label_row[keep]]]

    # Include eligible author–paper records outside the main sample.
    extra <- as.data.table(read_parquet(paste0(prefix, "-extra.parquet")))
    extra[, `:=`(author_id = as.numeric(author_id), work_id = as.numeric(work_id))]
    sample <- rbindlist(list(sample, extra), use.names = TRUE)
    rm(labels, authors, label_row, keep, extra)
    stopifnot(nrow(sample) == expected[[window]])

    models[[paste0("Fixed_", window)]] <- estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2
                                                | author_id + paper_year, sample)
    rm(sample)
}
print(models)
