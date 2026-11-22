# Table A3. Early diversity and the two components of team interdisciplinarity
library(arrow)
library(data.table)
library(fixest)

data_dir <- "model-data"  # Change this to the folder containing the extracted data.
setDTthreads(8)
arrow::set_cpu_count(8)

files <- list.files(file.path(data_dir, "mechanism"), "[.]parquet$", full.names = TRUE)
teams <- rbindlist(lapply(files, read_parquet))
teams[, `:=`(author_id = as.numeric(author_id), work_id = as.numeric(work_id))]

# Use the same weights and paper-level variables as the main regressions.
keys <- c("author_id", "work_id", "paper_year")
columns <- c(keys, "subfield_id", "paper_weight")
files <- list.files(file.path(data_dir, "main"), "[.]parquet$", full.names = TRUE)
papers <- rbindlist(lapply(files, function(file) read_parquet(file, col_select = all_of(columns))))
papers[, `:=`(author_id = as.numeric(author_id), work_id = as.numeric(work_id))]
teams <- papers[teams, on = keys]
stopifnot(nrow(teams) == 81712439L, !anyNA(teams$paper_weight))
rm(papers)
authors <- as.data.table(read_parquet(file.path(data_dir, "authors.parquet")))
author_row <- match(teams$author_id, as.numeric(authors$author_id))
stopifnot(!anyNA(author_row))

# Cut points for the five early-diversity groups, from lowest to highest.
cuts <- c(0.002835837268295671, 0.017090628509381846,
          0.03338824532905888, 0.05587396864968282)
group <- findInterval(authors$rao, cuts, left.open = TRUE) + 1L
teams[, `:=`(q = group[author_row], log_team_size = log1p(team_size))]
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

models <- list(
    Within = estimate(within_idr ~ i(q, ref = 1) + log_history + log_team_size
                      | paper_year + subfield_id, teams),
    Between = estimate(between_idr ~ i(q, ref = 1) + log_history + log_team_size
                       | paper_year + subfield_id, teams)
)
print(models)
