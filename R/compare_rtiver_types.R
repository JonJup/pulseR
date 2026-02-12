compare_river_types <- function(models) {
        
        # identify the global intersection of variables 
        shared <- Reduce(intersect, lapply(models, function(m) rownames(m$parameters$mean)))
        
        if (length(shared) == 0) {
                warning("No shared variables between all models")
                return(NA_real_)
        }
        
        # Build a list of all types (model + component combinations)
        types <- do.call(rbind, lapply(seq_along(models), function(i) {
                m <- models[[i]]
                K <- ncol(m$parameters$mean)
                data.frame(model = i, component = seq_len(K))
        }))

        # Give each type a label
        types$label <- paste0("M", types$model, "_K", types$component)
        
        n <- nrow(types)
        dist_mat <- matrix(NA, n, n, dimnames = list(types$label, types$label))
        
        for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                        m1 <- models[[types$model[i]]]
                        m2 <- models[[types$model[j]]]
                        
                        mu1    <- m1$parameters$mean[all_vars, types$component[i]]
                        sigma1 <- m1$parameters$variance$sigma[all_vars, all_vars, types$component[i]]
                        
                        mu2    <- m2$parameters$mean[all_vars, types$component[j]]
                        sigma2 <- m2$parameters$variance$sigma[all_vars, all_vars, types$component[j]]
                        
                        d <- bhattacharyya_dist(mu1, sigma1, mu2, sigma2)
                        dist_mat[i, j] <- d
                        dist_mat[j, i] <- d
                }
        }
        
        diag(dist_mat) <- 0
}


