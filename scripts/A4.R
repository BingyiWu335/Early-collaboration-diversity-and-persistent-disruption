# Table A4. Team interdisciplinarity and the main interaction
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
columns <- c(keys, "outcome", "idr_c", "idr_rao_c", "log_refs",
             "log_team_h", "age2", "paper_weight")
files <- list.files(file.path(data_dir, "main"), "[.]parquet$", full.names = TRUE)
papers <- rbindlist(lapply(files, function(file) read_parquet(file, col_select = all_of(columns))))
papers[, `:=`(author_id = as.numeric(author_id), work_id = as.numeric(work_id))]
teams <- papers[teams, on = keys]
stopifnot(nrow(teams) == 81712439L, !anyNA(teams$paper_weight))
rm(papers)
teams[, `:=`(within_c = within_idr - mean(within_idr),
            between_c = between_idr - mean(between_idr))]

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
    Base = estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2
                    | author_id + paper_year, teams),
    Within = estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2 +
                      within_c + idr_c:within_c | author_id + paper_year, teams),
    Between = estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2 +
                       between_c + idr_c:between_c | author_id + paper_year, teams),
    Joint = estimate(outcome ~ idr_c + idr_rao_c + log_refs + log_team_h + age2 +
                     within_c + idr_c:within_c + between_c + idr_c:between_c
                     | author_id + paper_year, teams)
)
print(models)
