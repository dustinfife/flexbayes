require_draw_renderer = function() {
  skip_if_not("prepare_flexplot" %in% getNamespaceExports("flexplot"), "New public flexplot renderer required")
}

prior_band_layers = function(p, geom) {
  Filter(function(layer) inherits(layer$geom, geom) && is.data.frame(layer$data) &&
    all(c(".prior_lower", ".prior_upper") %in% names(layer$data)), p$layers)
}

grid_manual_expected = function(x, prep, ids = seq_len(nrow(x$coefficient_draws))) {
  marginal = attr(prep$grid, "prior_marginal")
  X = model.matrix(delete.response(x$metadata$terms), marginal$cases,
    contrasts.arg = x$metadata$contrasts, xlev = x$metadata$factor_levels)
  linear = x$coefficient_draws[ids, , drop = FALSE] %*% t(X)
  values = switch(x$metadata$family$family, gaussian = linear, binomial = plogis(linear), poisson = exp(linear))
  result = vapply(marginal$groups, function(rows) rowMeans(values[, rows, drop = FALSE]), numeric(length(ids)))
  matrix(result, nrow = length(ids))
}
