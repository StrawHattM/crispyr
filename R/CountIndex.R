

CountIndex <- function(counts, order) {

  if(missing(order)) {

    withr::with_locale(LC_COLLATE = "en_US.UTF-8", sort(colnames(counts)[]))

  }

}



calculate_totals <- function(raw_counts, conditions, replicates) {

  sample_cols <- 3:ncol(raw_counts)

  totals <- data.frame(
    sample = factor(colnames(raw_counts)[sample_cols],
                    levels = colnames(raw_counts)[sample_cols]),
    total = sapply(sample_cols, function(i) sum(raw_counts[[i]])),
    condition = factor(conditions, levels = unique(conditions)),
    rep = factor(replicates, levels = unique(replicates))
  )

  return(totals)
}
