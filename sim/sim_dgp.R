# =========================================================================== *
# sim_dgp.R
#
# Data-generating process for the pulseR simulation study.
# One call to simulate_landscape() returns polygons with environmental
# columns plus the truth at all three levels.
#
#   L1  contiguous regions with non-convex (warped) boundaries; true graded
#       membership u* from distance to each region, width tau.
#   L2  local types inside each region. Types are either SHARED (same mean in
#       every region, so they should merge into one meta-type) or CONFINED
#       (unique to one region). share_types sets the proportion.
#   L3  mosaic zones that cut across regions. Each zone arranges the local
#       types by its own rule, at EQUAL composition, so only arrangement
#       (configuration) separates zones.
#
# All length scales are in units of h = mean cell spacing = 1/sqrt(N), so a
# setting means the same thing at every N.
#
# Requires: sf, spdep, igraph, RANN.
# =========================================================================== *

#' Simulate a landscape with three-level truth
#'
#' @param N Number of units (Voronoi cells). Ignored if `geometry` is given.
#' @param geometry Optional sf POLYGON object (e.g. real catchments). Its
#'   geometry is used instead of Voronoi cells; attributes are simulated.
#' @param jitter Voronoi irregularity, 0 = regular grid, 1 = strong.
#' @param R_true Number of true regions.
#' @param contrast Separation of region means. Pairwise distance between any
#'   two region means is contrast * sqrt(2) noise SDs.
#' @param tau Boundary transition width in cell spacings. 0 = crisp truth.
#' @param warp Amplitude of the sinusoidal coordinate warp that makes region
#'   boundaries non-convex (fraction of domain width).
#' @param p_reg Number of regional signal columns.
#' @param K_types Local types per region, scalar or length R_true.
#' @param share_types Proportion of each region's types that are shared.
#'   Confined types have region-specific means, so share_types < 1 puts
#'   REGION signal into the loc columns. For an L1 null use contrast = 0 AND
#'   share_types = 1.
#' @param delta Minimum pairwise distance between local type means (noise SDs).
#' @param p_loc Number of local-type signal columns.
#' @param zones Character vector of zone arrangements, one zone each:
#'   "block" (coarse patches), "fine" (cell-scale patches), "banded"
#'   (stripes), "random" (iid). length(zones) is the true number of zones.
#' @param grain Patch or stripe scale for "block" and "banded", in cell
#'   spacings.
#' @param n_irrelevant Number of pure-noise columns.
#' @param confound SD of two smooth spatial gradients unrelated to regions,
#'   added as columns conf1, conf2. 0 = no confound columns.
#' @param noise_sd Noise SD on every signal column.
#' @param type_dist Noise distribution on local-type columns.
#' @param seed Integer seed. The caller's RNG state is restored on exit.
#' @return List: polygons (sf with ID and env columns), env_cols (list of
#'   column names by block), truth (list), params (list).
simulate_landscape <- function(N            = 2000L,
                               geometry     = NULL,
                               jitter       = 0.6,
                               R_true       = 6L,
                               contrast     = 2.5,
                               tau          = 1,
                               warp         = 0.08,
                               p_reg        = 5L,
                               K_types      = 3L,
                               share_types  = 0.5,
                               delta        = 2,
                               p_loc        = 4L,
                               zones        = c("block", "fine", "banded", "random"),
                               grain        = 6,
                               n_irrelevant = 0L,
                               confound     = 0,
                               noise_sd     = 1,
                               type_dist    = c("gaussian", "skewed", "heavy"),
                               seed         = 1L) {
    type_dist <- match.arg(type_dist)
    zones     <- match.arg(zones, c("block", "fine", "banded", "random"),
                           several.ok = TRUE)
    params <- as.list(environment())
    params$geometry <- if (is.null(geometry)) "voronoi" else "user"

    .with_seed(seed, {
        ## ---- Geometry and adjacency -----------------------------------------*
        polys <- if (is.null(geometry)) .make_voronoi(N, jitter) else
            .prepare_geometry(geometry)
        N  <- nrow(polys)
        xy <- .unit_coords(polys)                 # centroids rescaled to [0, 1]
        h  <- 1 / sqrt(N)
        nb <- spdep::poly2nb(polys, queen = FALSE, snap = 1e-6)
        g  <- .nb_graph(nb)
        if (igraph::components(g)$no > 1L)
            warning("Simulated adjacency graph is disconnected; ",
                    "polygon_to_network() will bridge it.", call. = FALSE)

        ## ---- L1: regions ----------------------------------------------------*
        z     <- .make_regions(xy, R_true, warp, g, nb)
        R_true <- length(unique(z))
        d_reg <- .dist_to_regions(xy, z, R_true) / h          # N x R, cell units
        ustar <- .soft_membership(d_reg, z, tau)
        mu    <- .spread_means(R_true, p_reg, contrast * sqrt(2))
        X_reg <- ustar %*% mu + matrix(rnorm(N * p_reg, sd = noise_sd), N, p_reg)
        colnames(X_reg) <- paste0("reg", seq_len(p_reg))

        ## ---- L3 zones (needed before L2, they set the arrangement) ----------*
        M_true <- length(zones)
        zone   <- .nearest_seed(xy, .farthest_points(xy, M_true))

        ## ---- L2: local types ------------------------------------------------*
        K_r <- rep_len(as.integer(K_types), R_true)
        loc <- .make_types(xy, h, z, zone, zones, K_r, share_types, grain)
        ids   <- sort(unique(loc$type_id))
        m_loc <- .spread_means(length(ids), p_loc, delta)
        rownames(m_loc) <- ids
        E_loc <- switch(type_dist,
                        gaussian = matrix(rnorm(N * p_loc), N, p_loc),
                        skewed   = matrix(rexp(N * p_loc) - 1, N, p_loc),
                        heavy    = matrix(rt(N * p_loc, df = 3) / sqrt(3), N, p_loc))
        X_loc <- m_loc[as.character(loc$type_id), , drop = FALSE] + noise_sd * E_loc
        rownames(X_loc) <- NULL
        colnames(X_loc) <- paste0("loc", seq_len(p_loc))

        ## ---- Nuisance columns -------------------------------------------------*
        X_conf <- NULL
        if (confound > 0) {
            g1 <- .grf(xy, ell = 1 / 3, n_feat = 50L)          # domain-scale field
            g2 <- xy[, 1] * cos(pi / 5) + xy[, 2] * sin(pi / 5)
            X_conf <- cbind(conf1 = confound * as.numeric(scale(g1)),
                            conf2 = confound * as.numeric(scale(g2))) +
                matrix(rnorm(N * 2, sd = noise_sd), N, 2)
        }
        X_irr <- NULL
        if (n_irrelevant > 0) {
            X_irr <- matrix(rnorm(N * n_irrelevant, sd = noise_sd), N, n_irrelevant)
            colnames(X_irr) <- paste0("irr", seq_len(n_irrelevant))
        }

        ## ---- Reference classifiers (non-spatial Gaussian plug-ins) -----------*
        bayes_l1 <- .gauss_posterior(X_reg, mu, noise_sd)
        vstar    <- .local_posterior(X_loc, z, K_r, loc, m_loc, noise_sd)
    })

    env <- as.data.frame(do.call(cbind, Filter(Negate(is.null),
                                       list(X_reg, X_loc, X_conf, X_irr))))
    polys <- sf::st_sf(ID = sprintf("u%05d", seq_len(N)), env,
                       geometry = sf::st_geometry(polys))

    truth <- list(
        z          = z,                        # L1 hard label, 1..R
        ustar      = ustar,                    # N x R true graded membership
        d_boundary = .boundary_distance(d_reg, z),  # cells to nearest other region
        region_size = tabulate(z, R_true)[z],
        local      = loc$local,                # L2 type within region, 1..K_r
        joint      = paste0(z, ":", loc$local),  # L2 truth across regions
        type_id    = loc$type_id,              # meta-type truth (shared ids merge)
        shared     = loc$shared,               # unit's type is shared
        vstar      = vstar$post,               # N x max(K_r), within true region
        zone       = zone,                     # L3 truth, 1..M
        arrangement = zones[zone],
        K_r        = K_r,
        R_true     = R_true,
        M_true     = M_true,
        n_meta_true = length(ids),
        bayes_l1   = bayes_l1,                 # N x R posterior on reg columns
        bayes_l2   = vstar$hard,
        mu         = mu,
        m_loc      = m_loc,
        xy         = xy,
        h          = h)

    list(polygons = polys,
         env_cols = list(reg = colnames(X_reg), loc = colnames(X_loc),
                         conf = colnames(X_conf), irr = colnames(X_irr)),
         truth    = truth,
         params   = params)
}

# --------------------------------------------------------------------------- *
# Geometry
# --------------------------------------------------------------------------- *

# Jittered-grid Voronoi cells on a 100 km square in EPSG:3035, so that
# polygon_to_network() sees projected metre coordinates.
.make_voronoi <- function(N, jitter) {
    nx <- ceiling(sqrt(N)); ny <- ceiling(N / nx)
    gr <- expand.grid(i = seq_len(nx), j = seq_len(ny))
    gr <- gr[sort(sample.int(nrow(gr), N)), ]
    x  <- (gr$i - 0.5 + jitter * (runif(N) - 0.5)) / nx
    y  <- (gr$j - 0.5 + jitter * (runif(N) - 0.5)) / ny
    S  <- 1e5
    pts <- sf::st_as_sf(data.frame(x = x * S, y = y * S), coords = c("x", "y"),
                        crs = 3035)
    box <- sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = S, ymax = S),
                                     crs = sf::st_crs(3035)))
    v <- sf::st_collection_extract(
        sf::st_voronoi(sf::st_union(pts), envelope = box), "POLYGON")
    v <- sf::st_intersection(v, box)
    hit <- sf::st_intersects(pts, v)
    idx <- vapply(hit, function(k) if (length(k)) k[1L] else NA_integer_, integer(1))
    if (anyNA(idx) || anyDuplicated(idx))
        stop("Voronoi cell matching failed; change `seed` or `jitter`.", call. = FALSE)
    sf::st_sf(geometry = v[idx])
}

.prepare_geometry <- function(geometry) {
    stopifnot(inherits(geometry, "sf") || inherits(geometry, "sfc"))
    geo <- sf::st_geometry(geometry)
    if (isTRUE(sf::st_is_longlat(geo)))
        stop("Project `geometry` to a planar CRS first (e.g. EPSG:3035).", call. = FALSE)
    geo <- sf::st_make_valid(geo)
    geo <- geo[!sf::st_is_empty(geo)]
    sf::st_sf(geometry = geo)
}

# Centroids rescaled so the longer side spans [0, 1].
.unit_coords <- function(polys) {
    cc <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(polys))))
    cc <- cc[, 1:2, drop = FALSE]
    rng <- max(apply(cc, 2, function(v) diff(range(v))))
    sweep(cc, 2, apply(cc, 2, min)) / rng
}

.nb_graph <- function(nb) {
    from <- rep.int(seq_along(nb), lengths(nb))
    to   <- unlist(nb, use.names = FALSE)
    keep <- to > 0L & from < to
    igraph::make_graph(rbind(from[keep], to[keep]), n = length(nb), directed = FALSE)
}

# --------------------------------------------------------------------------- *
# L1
# --------------------------------------------------------------------------- *

.farthest_points <- function(xy, k) {
    idx <- sample.int(nrow(xy), 1L)
    d   <- sqrt(colSums((t(xy) - xy[idx, ])^2))
    while (length(idx) < k) {
        nxt <- which.max(d)
        idx <- c(idx, nxt)
        d   <- pmin(d, sqrt(colSums((t(xy) - xy[nxt, ])^2)))
    }
    xy[idx, , drop = FALSE]
}

.nearest_seed <- function(xy, seeds) {
    as.integer(RANN::nn2(seeds, xy, k = 1L)$nn.idx[, 1])
}

# Nearest seed in warped coordinates, then contiguity repair on the graph.
.make_regions <- function(xy, R, warp, g, nb) {
    ph <- runif(4, 0, 2 * pi); fr <- sample(1:2, 2, replace = TRUE)
    xw <- cbind(xy[, 1] + warp * sin(2 * pi * fr[1] * xy[, 2] + ph[1]),
                xy[, 2] + warp * sin(2 * pi * fr[2] * xy[, 1] + ph[2]))
    z <- .nearest_seed(xw, .farthest_points(xw, R))
    z <- .make_contiguous(z, g, nb)
    match(z, sort(unique(z)))                     # relabel to 1..R
}

# Keep the largest connected piece of each label; refill freed units from the
# majority label of their labelled neighbours. Repeat until stable.
.make_contiguous <- function(z, g, nb, max_iter = 20L) {
    for (it in seq_len(max_iter)) {
        changed <- FALSE
        for (r in unique(z[!is.na(z)])) {
            idx  <- which(z == r)
            comp <- igraph::components(igraph::induced_subgraph(g, idx))
            if (comp$no > 1L) {
                z[idx[comp$membership != which.max(comp$csize)]] <- NA
                changed <- TRUE
            }
        }
        while (anyNA(z)) {
            na <- which(is.na(z))
            for (i in na) {
                lab <- z[nb[[i]][nb[[i]] > 0L]]
                lab <- lab[!is.na(lab)]
                if (length(lab)) z[i] <- as.integer(names(which.max(table(lab))))
            }
            if (length(which(is.na(z))) == length(na)) break   # isolated units
        }
        if (!changed) break
    }
    z
}

# Distance from each unit to the nearest unit of each region (domain units);
# 0 for the unit's own region.
.dist_to_regions <- function(xy, z, R) {
    D <- matrix(0, nrow(xy), R)
    for (r in seq_len(R)) {
        inr <- which(z == r)
        D[, r] <- RANN::nn2(xy[inr, , drop = FALSE], xy, k = 1L)$nn.dists[, 1]
        D[inr, r] <- 0
    }
    D
}

# u*_ir proportional to exp(-d_ir / tau). Distances to other regions are
# shifted by half a cell so a boundary-adjacent unit is half a cell from the
# boundary, not one full cell. tau = 0 gives one-hot.
.soft_membership <- function(d, z, tau) {
    N <- nrow(d); R <- ncol(d)
    if (tau <= 0) { U <- matrix(0, N, R); U[cbind(seq_len(N), z)] <- 1; return(U) }
    dd <- pmax(d - 0.5, 0); dd[cbind(seq_len(N), z)] <- 0
    U  <- exp(-dd / tau)
    U / rowSums(U)
}

.boundary_distance <- function(d, z) {
    d[cbind(seq_len(nrow(d)), z)] <- Inf
    apply(d, 1, min)
}

# n points in p dimensions with minimum pairwise distance `min_dist`.
# Exact regular simplex (all distances equal) when p >= n - 1, rotated into a
# random subspace so every column carries signal. Otherwise repulsion on the
# sphere, rescaled so the closest pair sits at min_dist.
.spread_means <- function(n, p, min_dist) {
    if (n == 1L) return(matrix(0, 1L, p))
    if (p >= n - 1L) {
        M  <- diag(n) - 1 / n
        sv <- svd(M)
        V  <- sv$u[, seq_len(n - 1L), drop = FALSE] %*%
            diag(sv$d[seq_len(n - 1L)], n - 1L)
        Q  <- qr.Q(qr(matrix(rnorm(p * (n - 1L)), p, n - 1L)))
        X  <- V %*% t(Q)
    } else {
        X <- matrix(rnorm(n * p), n, p); X <- X / sqrt(rowSums(X^2))
        for (it in seq_len(300L)) for (i in seq_len(n)) {
            dv <- sweep(X[-i, , drop = FALSE], 2, X[i, ])
            X[i, ] <- X[i, ] - 0.01 * colSums(dv / (rowSums(dv^2) + 1e-9))
            X[i, ] <- X[i, ] / sqrt(sum(X[i, ]^2))
        }
    }
    dmin <- min(stats::dist(X))
    X * (min_dist / dmin)
}

# --------------------------------------------------------------------------- *
# L2 and L3
# --------------------------------------------------------------------------- *

# Stationary Gaussian random field by random Fourier features; ell in domain
# units. Returns a standardised vector.
.grf <- function(xy, ell, n_feat = 200L) {
    W <- matrix(rnorm(2 * n_feat, sd = 1 / ell), 2, n_feat)
    b <- runif(n_feat, 0, 2 * pi)
    f <- cos(xy %*% W + rep(b, each = nrow(xy))) %*% rep(1, n_feat)
    as.numeric(scale(f))
}

# One field per zone, by arrangement. Types are quantile bins of the field
# within each zone x region cell, so every cell has equal type composition and
# zones differ only in arrangement.
.make_types <- function(xy, h, z, zone, zones, K_r, share_types, grain) {
    N <- nrow(xy)
    f <- numeric(N)
    for (m in seq_along(zones)) {
        im <- which(zone == m)
        th <- runif(1, 0, pi)
        f[im] <- switch(zones[m],
            block  = .grf(xy, ell = grain * h)[im],
            fine   = .grf(xy, ell = 1 * h)[im],
            banded = sin(2 * pi * (xy[im, 1] * cos(th) + xy[im, 2] * sin(th)) /
                             (grain * h)) + 0.3 * .grf(xy, ell = grain * h)[im],
            random = runif(length(im)))
    }
    local <- integer(N)
    cell  <- interaction(zone, z, drop = TRUE)
    for (cl in levels(cell)) {
        ic <- which(cell == cl)
        K  <- K_r[z[ic[1]]]
        local[ic] <- pmin(K, 1L + floor(K * (rank(f[ic], ties.method = "random") - 1) /
                                            length(ic)))
    }
    n_sh   <- round(share_types * K_r)            # shared = lowest-numbered types
    shared <- local <= n_sh[z]
    type_id <- ifelse(shared, as.character(local), paste0("r", z, "_", local))
    list(local = local, type_id = type_id, shared = shared)
}

# --------------------------------------------------------------------------- *
# Reference classifiers
# --------------------------------------------------------------------------- *

.gauss_posterior <- function(X, means, sd) {
    D <- sapply(seq_len(nrow(means)), function(r) rowSums(sweep(X, 2, means[r, ])^2))
    L <- -D / (2 * sd^2)
    L <- exp(L - apply(L, 1, max))
    L / rowSums(L)
}

# Posterior over the true region's own types, true means, equal type weights.
.local_posterior <- function(X_loc, z, K_r, loc, m_loc, sd) {
    N <- nrow(X_loc); Kmax <- max(K_r)
    post <- matrix(NA_real_, N, Kmax)
    for (r in unique(z)) {
        ir  <- which(z == r)
        ids <- unique(loc$type_id[ir][order(loc$local[ir])])     # id of type 1..K_r
        P   <- .gauss_posterior(X_loc[ir, , drop = FALSE], m_loc[ids, , drop = FALSE], sd)
        post[ir, seq_len(ncol(P))] <- P
    }
    list(post = post, hard = max.col(replace(post, is.na(post), -1), "first"))
}

# --------------------------------------------------------------------------- *
# RNG helper: run `expr` with a seed, restore the caller's stream afterwards.
# --------------------------------------------------------------------------- *
.with_seed <- function(seed, expr) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
            else assign(".Random.seed", old, globalenv()))
    set.seed(seed)
    eval.parent(substitute(expr))
}
