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

  # Transition location. Prefer the requested objective channel (GR) when it has
  # adequate finite coverage; otherwise fall back to a finite-safe torque/ROP proxy.
  # The source actually used is recorded, and a degenerate/non-finite transition is a
  # hard error (so the validation split can never silently change on missing data).
  .zfin <- function(v) {                                   # finite-safe standardization
    fin <- is.finite(v); m <- mean(v[fin]); s <- stats::sd(v[fin])
    if (!is.finite(s) || s <= 0) s <- 1
    out <- (v - m) / s; out[!is.finite(out)] <- 0; out
  }
  min_cov <- 0.80
  use_channel <- !is.null(transition_channel) &&
    isTRUE(mean(is.finite(d[[transition_channel]])) >= min_cov)
  if (use_channel) {
    ch <- d[[transition_channel]]
    if (anyNA(ch) || any(!is.finite(ch)))                  # fill sparse gaps with the channel median
      ch[!is.finite(ch)] <- stats::median(ch[is.finite(ch)])
    trans <- locate_transition(ch, k = k)
    transition_source <- transition_channel
  } else {
    # Proxy: combine standardized |d(torque)/dt| and |d(ROP)/dt| step magnitude.
    dtq <- abs(c(0, diff(d$tqa)))
    drp <- abs(c(0, diff(y_raw)))
    proxy <- .zfin(dtq) + .zfin(drp)
    trans <- locate_transition(proxy, k = k)
    transition_source <- if (is.null(transition_channel)) "torque_rop_proxy"
                         else sprintf("torque_rop_proxy (%s coverage %.2f < %.2f)",
                                      transition_channel, mean(is.finite(d[[transition_channel]])), min_cov)
  }
  if (length(trans) != 1L || !is.finite(trans) || trans < 1L || trans > nrow(d)) {
    stop("load_volve_drilling: failed to locate a single finite transition index (source: ",
         transition_source, ").")
  }

  list(
    X = X, y = y,
    regime = as.matrix(d[, c("stick_rt", "arc_gr_rt"), with = FALSE]),
    depth = d$dept, transition = as.integer(trans),
    transition_source = transition_source, n = nrow(d),
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

#' Resolve the local Ku_NG drilling data root
#'
#' Set \code{KU_NG_DATA_ROOT} to rerun the evidence pipeline on a local data
#' checkout without editing scripts. If the environment variable is unset, the
#' function falls back to a relative \code{data} directory in the current
#' working directory.
#'
#' @param root optional explicit data root. If \code{NULL}, uses
#'   \code{Sys.getenv("KU_NG_DATA_ROOT")} and then \code{"data"}.
#' @return normalized path string; the path is not required to exist.
#' @export
ku_ng_data_root <- function(root = NULL) {
  if (is.null(root) || !nzchar(root)) {
    root <- Sys.getenv("KU_NG_DATA_ROOT", unset = "data")
  }
  normalizePath(root, mustWork = FALSE)
}

.drilling_dataset_catalog <- function() {
  data.frame(
    dataset_id = c(
      "forge58_processed",
      "forge58_raw",
      "forge56_10sec",
      "forge56_1sec",
      "forge78b_10sec",
      "forge78b_1sec",
      "volve_onbottom",
      "volve_real",
      "volve_synthetic"
    ),
    relative_path = c(
      "forge/Well_58-32_processed_pason_log.csv",
      "forge/Well_58-32_raw_pason_log.csv",
      "forge/well-56-32/56-32_10sec_standard.csv",
      "forge/well-56-32/56-32_1sec_standard.csv",
      "forge/well-78b-32/78b-32_10sec_standard.csv",
      "forge/well-78b-32/78b-32_1sec_standard.csv",
      "processed/volve_onbottom.parquet",
      "processed/volve_15_9_F15_real.parquet",
      "processed/volve_synthetic.parquet"
    ),
    role = c(
      "primary_compact_depth_log",
      "primary_large_time_log",
      "external_validation",
      "secondary_sparse_1sec",
      "secondary_sparse_stress",
      "secondary_sparse_1sec",
      "primary_oilfield_onbottom",
      "primary_oilfield_compact",
      "synthetic_not_publication_evidence"
    ),
    stringsAsFactors = FALSE
  )
}

.count_file_lines <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    x <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(x)) break
    n <- n + length(x)
  }
  n
}

#' Build a lightweight manifest for the local drilling data catalog
#'
#' This helper intentionally records only the project-relevant datasets, not
#' every PDF/DLIS sidecar in the raw archives. Use the Paper 1B evidence script
#' when a complete recursive manifest with hashes is needed.
#'
#' @param root Ku_NG data root.
#' @param include_hash if \code{TRUE}, add md5 hashes via \code{tools::md5sum}.
#' @param count_rows if \code{TRUE}, count text-file rows. Disabled by default
#'   because the raw Pason logs are large.
#' @return data.frame with dataset id, path, existence, size, and optional row
#'   counts / hashes.
#' @export
drilling_data_manifest <- function(root = ku_ng_data_root(),
                                   include_hash = FALSE,
                                   count_rows = FALSE) {
  catalog <- .drilling_dataset_catalog()
  rows <- lapply(seq_len(nrow(catalog)), function(i) {
    path <- file.path(root, catalog$relative_path[i])
    exists <- file.exists(path)
    info <- if (exists) file.info(path) else NULL
    ext <- tolower(tools::file_ext(path))
    line_count <- NA_integer_
    if (exists && count_rows && ext %in% c("csv", "las", "lis", "asc")) {
      line_count <- .count_file_lines(path)
    }
    hash <- NA_character_
    if (exists && include_hash) hash <- unname(tools::md5sum(path))
    data.frame(
      dataset_id = catalog$dataset_id[i],
      role = catalog$role[i],
      relative_path = catalog$relative_path[i],
      path = path,
      exists = exists,
      format = ext,
      size_bytes = if (exists) as.numeric(info$size) else NA_real_,
      line_count = line_count,
      md5 = hash,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.clean_drilling_numeric <- function(x, sentinel = -9990) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out) | out <= sentinel] <- NA_real_
  out
}

.zsafe <- function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  m <- mean(v, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  out <- (v - m) / s
  out[!is.finite(out)] <- 0
  out
}

.forge_standard_path <- function(root, well, cadence) {
  well <- match.arg(well, c("56-32", "78B-32"))
  cadence <- match.arg(cadence, c("10sec", "1sec"))
  if (well == "56-32") {
    file.path(root, "forge", "well-56-32",
              sprintf("56-32_%s_standard.csv", cadence))
  } else {
    file.path(root, "forge", "well-78b-32",
              sprintf("78b-32_%s_standard.csv", cadence))
  }
}

#' Load a Utah FORGE standard CSV as a drilling regression frame
#'
#' Supports the FORGE 56-32 / 78B-32 standard Pason exports used as the third-
#' well evidence candidates for Paper 1B. The files contain a second units row;
#' it is removed automatically. The returned object mirrors
#' \code{\link{load_volve_drilling}} / \code{\link{load_forge_drilling}}.
#'
#' @param path optional explicit CSV path. If \code{NULL}, constructed from
#'   \code{root}, \code{well}, and \code{cadence}.
#' @param well \code{"56-32"} or \code{"78B-32"}.
#' @param cadence \code{"10sec"} or \code{"1sec"}.
#' @param root Ku_NG data root.
#' @param target \code{"log_rop"} (default) or \code{"rop"}.
#' @param k transition smoothing half-width.
#' @param rop_cap maximum plausible ROP in ft/hour for quality filtering.
#' @param min_optional_frac minimum non-missing fraction required before an
#'   optional channel is admitted as a predictor.
#' @return list with \code{X}, \code{y}, \code{regime}, \code{depth},
#'   \code{transition}, \code{n}, and metadata fields.
#' @export
load_forge_standard_drilling <- function(path = NULL,
                                         well = c("56-32", "78B-32"),
                                         cadence = c("10sec", "1sec"),
                                         root = ku_ng_data_root(),
                                         target = c("log_rop", "rop"),
                                         k = 200L,
                                         rop_cap = 1000,
                                         min_optional_frac = 0.80) {
  well <- match.arg(well)
  cadence <- match.arg(cadence)
  target <- match.arg(target)
  if (is.null(path)) path <- .forge_standard_path(root, well, cadence)
  if (!file.exists(path)) stop("FORGE standard CSV not found: ", path)

  d <- data.table::fread(path, na.strings = c("", "NA", "N/A"))
  if (nrow(d) > 0 && "Date" %in% names(d) &&
      identical(tolower(as.character(d$Date[1])), "date/time")) {
    d <- d[-1L, ]
  }

  need <- c(
    "Bit Depth", "Depth Hole Total Vertical Depth",
    "Rate of Penetration (Depth/Hour)", "Weight on Bit",
    "Rotary Revolutions per Minute"
  )
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("FORGE standard column(s) not found: ", paste(miss, collapse = ", "))

  get <- function(nm) .clean_drilling_numeric(d[[nm]])
  frame <- data.frame(
    dept = get("Bit Depth"),
    tvd = get("Depth Hole Total Vertical Depth"),
    rop = get("Rate of Penetration (Depth/Hour)"),
    wob = get("Weight on Bit"),
    rpm = get("Rotary Revolutions per Minute"),
    stringsAsFactors = FALSE
  )
  optional_map <- c(
    diff_pressure = "Differential Pressure",
    flow = "Flow In",
    top_drive_torque = "Top Drive Torque",
    downhole_torque = "Downhole Torque",
    gamma = "Gamma Measured while Drilling",
    stick = "Stick Slip"
  )
  for (nm in names(optional_map)) {
    src <- optional_map[[nm]]
    if (src %in% names(d)) frame[[nm]] <- get(src)
  }

  ok <- is.finite(frame$rop) & frame$rop > 0 & frame$rop < rop_cap &
    is.finite(frame$dept) & is.finite(frame$tvd) &
    is.finite(frame$wob) & frame$wob >= 0 &
    is.finite(frame$rpm) & frame$rpm > 0
  frame <- frame[ok, , drop = FALSE]
  frame <- frame[order(frame$dept), , drop = FALSE]
  if (!nrow(frame)) stop("No usable FORGE standard rows after quality filtering: ", path)

  predictors <- c("wob", "rpm", "dept", "tvd")
  optional_predictors <- setdiff(names(frame), c("rop", predictors))
  for (nm in optional_predictors) {
    frac <- mean(is.finite(frame[[nm]]))
    if (is.finite(frac) && frac >= min_optional_frac) {
      predictors <- c(predictors, nm)
    }
  }
  keep <- stats::complete.cases(frame[, c("rop", predictors), drop = FALSE])
  frame <- frame[keep, , drop = FALSE]
  if (nrow(frame) < 30) stop("Too few usable rows for FORGE standard frame: ", nrow(frame))

  X <- as.matrix(frame[, predictors, drop = FALSE])
  storage.mode(X) <- "double"
  y <- if (target == "log_rop") log(pmax(frame$rop, 1e-6)) else frame$rop

  proxy <- .zsafe(abs(c(0, diff(frame$rop)))) +
    .zsafe(abs(c(0, diff(frame$wob)))) +
    .zsafe(abs(c(0, diff(frame$rpm))))
  if ("top_drive_torque" %in% names(frame) &&
      sum(is.finite(frame$top_drive_torque)) >= 30) {
    proxy <- proxy + .zsafe(abs(c(0, diff(frame$top_drive_torque))))
  }
  trans <- locate_transition(proxy, k = k)

  regime_cols <- intersect(c("diff_pressure", "flow", "top_drive_torque",
                             "downhole_torque", "gamma", "stick"),
                           names(frame))
  regime <- if (length(regime_cols)) {
    as.matrix(frame[, regime_cols, drop = FALSE])
  } else {
    cbind(rop_step = abs(c(0, diff(frame$rop))))
  }

  list(
    X = X, y = y, frame = frame,
    regime = regime, depth = frame$dept, transition = as.integer(trans),
    n = nrow(frame), predictors = predictors, target = target,
    segment = cadence, dataset = paste0("utah_forge_", gsub("[^A-Za-z0-9]+", "_", tolower(well))),
    source_path = path
  )
}
