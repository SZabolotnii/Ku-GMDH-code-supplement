# MIA-GMDH tournament (algorithm-spec.md 3). The outer architecture is the
# classical multilayer iterative algorithm; only INNER_ESTIMATE and the
# (optional) PMM-loss criterion differ from classical GMDH.
#
# Each node is either an input variable (layer 0) or a fitted KG-2 partial
# model (layer >= 1) referencing two parent nodes by id. Node ids increase with
# layer, so a child's parents always have smaller ids -> the cascade can be
# evaluated by a single id-ordered sweep (see predict.gmdh_pmm).

#' Fit a GMDH-PMM model
#'
#' @param X numeric matrix or data.frame of predictors (n x d).
#' @param y numeric response (length n).
#' @param control a \code{\link{gmdh_pmm_control}} list.
#' @return an object of class \code{gmdh_pmm}.
#' @export
gmdh_pmm <- function(X, y, control = gmdh_pmm_control()) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X); d <- ncol(X)
  if (d < 2) stop("GMDH needs at least 2 predictors.")
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(d))

  # train / validation split. Most experiments use a random split; time-aware
  # C5 experiments can supply fixed indices through the control object.
  if (!is.null(control$train_index) || !is.null(control$val_index)) {
    if (is.null(control$train_index) || is.null(control$val_index)) {
      stop("train_index and val_index must be supplied together.")
    }
    train <- sort(unique(as.integer(control$train_index)))
    val <- sort(unique(as.integer(control$val_index)))
    if (length(train) < 6L) stop("Training partition must contain at least 6 rows.")
    if (length(val) < 1L) stop("Validation partition is empty.")
    if (any(train < 1L | train > n) || any(val < 1L | val > n)) {
      stop("train_index/val_index outside data row range.")
    }
    if (length(intersect(train, val))) stop("train_index and val_index must not overlap.")
  } else {
    if (!is.null(control$seed)) set.seed(control$seed)
    n_train <- max(6L, floor(control$split_ratio * n))
    train <- sort(sample.int(n, n_train))
    val <- setdiff(seq_len(n), train)
  }
  if (length(val) < 1) stop("Validation partition is empty; lower split_ratio.")
  y_train <- y[train]; y_val <- y[val]

  # node bookkeeping, indexed by id
  nodes <- vector("list", d)
  out_train <- vector("list", d)   # node output on the training partition
  out_val <- vector("list", d)     # node output on the validation partition
  for (v in seq_len(d)) {
    nodes[[v]] <- list(id = v, layer = 0L, type = "input", var = v,
                       name = colnames(X)[v])
    out_train[[v]] <- X[train, v]
    out_val[[v]] <- X[val, v]
  }

  current_ids <- seq_len(d)
  next_id <- d + 1L
  cr_best_prev <- Inf
  best <- NULL
  layers <- list()

  for (ell in seq_len(control$L_max)) {
    if (length(current_ids) < 2) break
    pairs <- utils::combn(current_ids, 2)
    cand <- vector("list", ncol(pairs))

    for (k in seq_len(ncol(pairs))) {
      a <- pairs[1, k]; b <- pairs[2, k]
      est <- inner_estimate(out_train[[a]], out_train[[b]], y_train, control)
      yhat_val <- kg2_predict(est$theta, out_val[[a]], out_val[[b]])
      cr <- external_criterion(y_val, yhat_val, method = est$method,
                               diag = est$diag, type = control$criterion,
                               reserve_weight = control$reserve_weight,
                               tail_weight = control$tail_weight,
                               overreserve_weight = control$overreserve_weight,
                               tail_prob = control$tail_prob)
      cand[[k]] <- list(id = next_id, layer = ell, type = "model",
                        parents = c(a, b), theta = est$theta,
                        method = est$method, converged = est$converged,
                        diag = est$diag, cr = cr)
      next_id <- next_id + 1L
    }

    crs <- vapply(cand, function(z) z$cr, numeric(1))
    crs[!is.finite(crs)] <- Inf
    ord <- order(crs)
    keep <- ord[seq_len(min(control$F, length(ord)))]
    selected <- cand[keep]

    for (nd in selected) {
      nodes[[nd$id]] <- nd
      out_train[[nd$id]] <- kg2_predict(nd$theta, out_train[[nd$parents[1]]],
                                        out_train[[nd$parents[2]]])
      out_val[[nd$id]] <- kg2_predict(nd$theta, out_val[[nd$parents[1]]],
                                      out_val[[nd$parents[2]]])
    }

    cr_best_curr <- crs[keep[1]]
    if (is.null(best) || cr_best_curr < best$cr) best <- selected[[1]]
    layers[[ell]] <- list(
      best_cr = cr_best_curr,
      best_id = selected[[1]]$id,    # best node of this layer (for depth curves)
      methods = table(vapply(selected, function(z) z$method, character(1)))
    )

    # relative-improvement stopping (control$epsilon of previous best)
    improved <- (cr_best_prev - cr_best_curr) > control$epsilon * abs(cr_best_prev)
    if (ell > 1 && !improved) break
    cr_best_prev <- cr_best_curr
    current_ids <- vapply(selected, function(z) z$id, numeric(1))
  }

  structure(
    list(nodes = nodes, best = best, best_id = best$id,
         d = d, var_names = colnames(X), control = control,
         n_train = length(train), n_val = length(val), layers = layers),
    class = "gmdh_pmm"
  )
}

#' Predict from a fitted GMDH-PMM model
#'
#' @param object a \code{gmdh_pmm} object.
#' @param newdata numeric matrix/data.frame with the same \code{d} columns
#'   (in the same order) as the training predictors.
#' @param node node id whose output to return (default: the global best). Use
#'   \code{object$layers[[L]]$best_id} to read off the best model at depth L.
#' @param ... ignored.
#' @return numeric vector of predictions.
#' @export
predict.gmdh_pmm <- function(object, newdata, node = object$best_id, ...) {
  X <- as.matrix(newdata)
  storage.mode(X) <- "double"
  if (ncol(X) != object$d) stop("newdata must have ", object$d, " columns.")
  out <- vector("list", length(object$nodes))
  for (nd in object$nodes) {
    if (is.null(nd)) next                       # unselected candidate id -> gap
    if (nd$type == "input") {
      out[[nd$id]] <- X[, nd$var]
    } else {
      out[[nd$id]] <- kg2_predict(nd$theta, out[[nd$parents[1]]],
                                  out[[nd$parents[2]]])
    }
  }
  out[[node]]
}

#' @export
print.gmdh_pmm <- function(x, ...) {
  cat("GMDH-PMM model\n")
  cat(sprintf("  predictors: %d | train/val: %d/%d\n",
              x$d, x$n_train, x$n_val))
  cat(sprintf("  layers grown: %d | best layer: %d | best CR: %.5g\n",
              length(x$layers), x$best$layer, x$best$cr))
  methods <- vapply(Filter(function(z) !is.null(z) && z$type == "model", x$nodes),
                    function(z) z$method, character(1))
  if (length(methods)) {
    cat("  partial-model methods selected: ")
    cat(paste(names(table(methods)), table(methods), sep = "=", collapse = ", "))
    cat("\n")
  }
  invisible(x)
}
