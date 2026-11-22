# Table A2. Main interaction, pooled and by domain
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
authors <- as.data.table(read_parquet(file.path(data_dir, "authors.parquet")))
author_row <- match(papers$author_id, as.numeric(authors$author_id))
stopifnot(!anyNA(author_row))
papers[, domain_id := authors$domain_id[author_row]]
rm(author_row)

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

domains <- c(Pooled = NA, Physical = 1L, Life = 2L, Health = 3L, Social = 4L)
models <- list()
for (label in names(domains)) {
    id <- domains[[label]]
    sample <- if (is.na(id)) papers else papers[domain_id == id]

    # idr_rao_c is the interaction of paper interdisciplinarity and early diversity.
    models[[label]] <- estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2
                               | author_id + paper_year, sample)
    rm(sample)
}
print(models)
