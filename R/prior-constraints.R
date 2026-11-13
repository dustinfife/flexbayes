# A comparison specifies a cut in the ordered CDF:
# P(X < t) precedes P(X <= t), allowing atoms at t.
# Upper-tail events are complements of these two different cuts.
prior_check_probability_constraints = function(records) {
  indices = which(vapply(records, function(r) r$type == "probability", logical(1)))
  if (!length(indices)) return(list(contradictory = FALSE))
  cuts = lapply(indices, function(i) {
    r = records[[i]]
    lower = r$operator %in% c("<", "<=")
    probability = if (lower) r$target else 1 - r$target
    list(index = i, threshold = r$threshold,
      inclusive = r$operator %in% c("<=", ">"),
      probability = probability, lower = lower,
      minimum = max(0, probability - r$tolerance),
      maximum = min(1, probability + r$tolerance))
  })
  for (cut in cuts) {
    r = records[[cut$index]]
    if (!is.finite(r$target) || r$target < 0 || r$target > 1) {
      return(list(contradictory = TRUE, conflicts = cut$index,
        message = paste0("Belief ", cut$index, " is contradictory: probabilities must be between 0 and 1.")))
    }
  }
  if (length(cuts) < 2L) return(list(contradictory = FALSE))
  for (a in seq_len(length(cuts) - 1L)) {
    for (b in seq.int(a + 1L, length(cuts))) {
      left = cuts[[a]]
      right = cuts[[b]]
      equivalent = left$threshold == right$threshold && left$inclusive == right$inclusive
      ordered = left$threshold < right$threshold ||
        (left$threshold == right$threshold && !left$inclusive && right$inclusive)
      if (!equivalent && !ordered) {
        swap = left
        left = right
        right = swap
      }
      conflict = if (equivalent) max(left$minimum, right$minimum) > min(left$maximum, right$maximum) else left$minimum > right$maximum
      if (!conflict) next
      numbers = c(indices[a], indices[b])
      reason = if (equivalent) {
        "equivalent probability statements disagree."
      } else if (left$lower && !right$lower) {
        total = records[[left$index]]$target + records[[right$index]]$target
        sprintf("they assign a total probability of %.2f to mutually exclusive regions.", total)
      } else {
        "CDF values must be nondecreasing (P(value < t) <= P(value <= t))."
      }
      return(list(contradictory = TRUE, conflicts = numbers,
        message = paste0("Beliefs ", numbers[1], " and ", numbers[2], " are contradictory: ", reason)))
    }
  }
  list(contradictory = FALSE)
}
