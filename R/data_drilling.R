# Drilling-dataset loaders for the Paper 1 ReserveTW4 replication.
#
# Two real petroleum/geothermal benchmarks replace the UCI Gas Turbine CO
# flagship of the arXiv v1 draft:
#   - Volve well 15/9-F-15 (Equinor Open) -- time-based MWD, has a gamma-ray
#     channel for an OBJECTIVE lithology-transition reserve window.
#   - Utah FORGE Pason logs (CC BY 4.0) -- depth-indexed; no GR, so the
#     formation transition is located by a torque/ROP step PROXY.
#
# Both return a common frame: (X, y, regime, depth, transition_index) where
# the rows are ordered along drilling progression (depth) so a temporal /
# along-hole split is meaningful and the reserve window sits past the
# structural break. The regression frame is documented in data-frame-spec.md.

# Smoothed two-sided step statistic: mean of the k points ahead minus the mean
# of the k points behind. Used to locate the dominant formation transition.
.step_statistic <- function(x, k) {
  n <- length(x)
  s <- rep(NA_real_, n)
  if (n <= 2L * k) return(s)
  for (i in (k + 1L):(n - k)) {
    s[i] <- mean(x[(i + 1L):(i + k)]) - mean(x[(i - k):(i - 1L)])
  }
  s
}

#' Locate the dominant structural break in an ordered channel
#'
#' Returns the row index of the largest-magnitude smoothed step in
#' \code{channel}, restricted to the central \code{[min_frac, max_frac]} of the
#' series so the reserve window keeps a usable number of rows on both sides.
#'
#' @param channel numeric vector ordered along drilling progression.
#' @param k half-width of the smoothing window (default 200).
#' @param min_frac,max_frac admissible location range for the break
#'   (defaults 0.55 / 0.9: the break must leave >=55\% for train/cal/val and
#'   >=10\% for the reserve window).
#' @return integer row index of the detected transition.
#' @export
locate_transition <- function(channel, k = 200L, min_frac = 0.55, max_frac = 0.90) {
  n <- length(channel)
  s <- abs(.step_statistic(channel, k))
  lo <- max(k + 1L, floor(min_frac * n))
  hi <- min(n - k, ceiling(max_frac * n))
  if (hi <= lo) return(floor(0.8 * n))
  band <- rep(NA_real_, n)
  band[lo:hi] <- s[lo:hi]
  which.max(band)
}

#' Load the Volve drilling regression frame
#'
#' @param path CSV exported from \code{volve_onbottom.parquet}.
#' @param segment \code{"gr"} (default) restricts to the gamma-ray-bearing LAS
#'   file (well 15/9-F-15, ~1360-2536 m, 27826 on-bottom rows) so the lithology
#'   transition is labelled objectively by ARC_GR_RT. \code{"full"} keeps all
#'   four LAS files ordered by depth and locates the transition by a torque/ROP
#'   proxy (no GR for files 2-4).
#' @param target \code{"log_rop5"} (default) models \eqn{\log} ROP5 -- the raw
#'   ROP5 is artifact-dominated (skew 15, excess kurtosis 227); the log target
#'   keeps the genuine right-skewed, heavy-tailed residual structure
#'   (\eqn{\gamma_3\approx2.8}, \eqn{g_2\approx0.63}) that motivates PMM2 while
#'   staying numerically tractable. \code{"rop5"} keeps the raw response.
#' @param k transition smoothing half-width (default 200).
#' @return list with \code{X} (n x 4: swob, rpm, tqa, dept), \code{y},
#'   \code{regime} (stick_rt, arc_gr_rt), \code{depth}, \code{transition}
#'   (integer row index of the formation transition), \code{n}, and meta fields.
#' @export
load_volve_drilling <- function(path = "data/volve_onbottom.csv",
                                segment = c("gr", "full"),
                                target = c("log_rop5", "rop5"),
                                k = 200L) {
  segment <- match.arg(segment)
  target <- match.arg(target)
  if (!file.exists(path)) stop("Volve CSV not found: ", path)
  d <- data.table::fread(path)

  if (segment == "gr") {
    d <- d[d$source_file == "WL_RAW_BHPR-GR-MECH_TIME_MWD_1.LAS", ]
    transition_channel <- "arc_gr_rt"
  } else {
    transition_channel <- NULL  # GR absent in files 2-4 -> torque/ROP proxy below
  }
  d <- d[order(d$dept), ]

  predictors <- c("swob", "rpm", "tqa", "dept")
  X <- as.matrix(d[, predictors, with = FALSE])
  storage.mode(X) <- "double"
  y_raw <- d$rop5
  y <- if (target == "log_rop5") log(pmax(y_raw, 1e-6)) else y_raw

  if (!is.null(transition_channel) && all(is.finite(d[[transition_channel]]))) {
    trans <- locate_transition(d[[transition_channel]], k = k)
  } else {
    # Proxy: combine standardized |d(torque)/dt| and |d(ROP)/dt| step magnitude.
    dtq <- abs(c(0, diff(d$tqa)))
    drp <- abs(c(0, diff(y_raw)))
    z <- function(v) (v - mean(v)) / stats::sd(v)
    proxy <- z(dtq) + z(drp)
    trans <- locate_transition(proxy, k = k)
  }

  list(
    X = X, y = y,
    regime = as.matrix(d[, c("stick_rt", "arc_gr_rt"), with = FALSE]),
    depth = d$dept, transition = as.integer(trans), n = nrow(d),
    predictors = predictors, target = target, segment = segment,
    dataset = "volve_15_9_F15"
  )
}

#' Load a Utah FORGE drilling regression frame (Pason log schema)
#'
#' FORGE Pason logs are depth-indexed and carry no gamma-ray channel, so the
#' formation transition is located by a torque/ROP step proxy. The regression
#' frame mirrors Volve: target = log ROP, predictors = WOB, RPM, torque, depth;
#' hookload is excluded (controller-coupled, analogous to the Volve block-
#' velocity exclusion).
#'
#' @param path CSV (e.g. \code{Well_58-32_processed_pason_log.csv}).
#' @param target \code{"log_rop"} (default) or \code{"rop"}.
#' @param k transition smoothing half-width (default 60; FORGE processed logs
#'   are ~7300 rows, coarser than the Volve MWD stream).
#' @return list with the same fields as \code{\link{load_volve_drilling}}
#'   (\code{regime} carries the torque/ROP proxy and depth).
#' @export
load_forge_drilling <- function(path = "data/Well_58-32_processed_pason_log.csv",
                                target = c("log_rop", "rop"), k = 60L) {
  target <- match.arg(target)
  if (!file.exists(path)) stop("FORGE CSV not found: ", path)
  d <- data.table::fread(path)

  col <- function(nm) {
    hit <- names(d)[which(names(d) == nm)]
    if (!length(hit)) stop("FORGE column not found: ", nm)
    as.numeric(d[[hit]])
  }
  depth <- col("Depth(ft)")
  rop   <- col("ROP(1 ft)")
  wob   <- col("weight on bit (k-lbs)")
  rpm   <- col("Rotary Speed (rpm)")
  tq    <- col("Surface Torque (psi)")

  frame <- data.frame(wob = wob, rpm = rpm, tq = tq, dept = depth,
                      rop = rop)
  ok <- stats::complete.cases(frame) & is.finite(rop) & rop > 0 &
        is.finite(rpm) & rpm > 0           # drop off-bottom / non-rotating rows
  frame <- frame[ok, , drop = FALSE]
  frame <- frame[order(frame$dept), ]

  X <- as.matrix(frame[, c("wob", "rpm", "tq", "dept")])
  storage.mode(X) <- "double"
  y <- if (target == "log_rop") log(pmax(frame$rop, 1e-6)) else frame$rop

  dtq <- abs(c(0, diff(frame$tq)))
  drp <- abs(c(0, diff(frame$rop)))
  z <- function(v) { s <- stats::sd(v); if (!is.finite(s) || s <= 0) s <- 1; (v - mean(v)) / s }
  proxy <- z(dtq) + z(drp)
  trans <- locate_transition(proxy, k = k)

  list(
    X = X, y = y,
    regime = cbind(torque_step = dtq, rop_step = drp),
    depth = frame$dept, transition = as.integer(trans), n = nrow(frame),
    predictors = c("wob", "rpm", "tq", "dept"), target = target,
    segment = "processed", dataset = "utah_forge_58_32"
  )
}
