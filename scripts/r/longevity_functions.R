#### Functions used to compute longevity and hazard rates ----

fit_functions <- function(longevities, survival, set_func = "all") {
  
  if (any(longevities <= 0))
    stop("Longevities must be strictly positive.")
  
  surv_obj <- Surv(longevities, survival)
  
  func_names <- c(
    "exponential", "weibull", "lognormal",
    "gompertz", "gamma", "llogis", "gengamma"
  )
  
  results <- list()
  
  # Fit all
  if (identical(set_func, "all")) {
    
    for(func in func_names) {
      
      cat('Fitting', func, '...')
      
      fit <- tryCatch(
        flexsurvreg(surv_obj ~ 1, dist = func),
        error = function(e) {
          cat(' ✗\n')
          return(NULL)
        }
      )
      
      if(!is.null(fit)) {
        results[[func]] <- fit
        cat(' ✓\n')
      }
    }
    
    
  } else {
    
    # Validate requested distributions
    if (!all(set_func %in% func_names)) {
      stop(
        "Invalid distribution name(s).\n",
        "Allowed values are: 'all' or any of:\n  ",
        paste(func_names, collapse = ", ")
      )
    }
    
    for(func in set_func) {
      
      cat('Fitting', func, '...')
      
      fit <- tryCatch(
        flexsurvreg(surv_obj ~ 1, dist = func),
        error = function(e) {
          cat(' ✗\n')
          return(NULL)
        }
      )
      
      if(!is.null(fit)) {
        results[[func]] <- fit
        cat(' ✓\n')
      }
    }
    
  }
  
  return(results)
}


compare_functions <- function(longevities, survival, set_func = "all") {
  
  function_list <- fit_functions(
    longevities = longevities,
    survival    = survival,
    set_func    = set_func
  )
  # function_list <- function_list[!sapply(function_list, is.null)]
  # Extract model names
  model_names <- names(function_list)
  
  # Extract AIC and number of parameters dynamically
  AIC_values <- sapply(function_list, function(x) x$AIC)
  N_params   <- sapply(function_list, function(x) length(x$res.t))
  
  comparison <- data.frame(
    Model    = model_names,
    AIC      = AIC_values,
    N_params = N_params,
    stringsAsFactors = FALSE
  )
  
  # Sort by AIC and calculate Delta AIC and weights
  comparison <- comparison |>
    arrange(AIC) |>
    mutate(
      Delta_AIC = AIC - min(AIC),
      Weight    = exp(-0.5 * Delta_AIC) /
        sum(exp(-0.5 * Delta_AIC)),
      Rank      = row_number()
    )
  
  return(comparison)
}

## sample remaining life using a gengamma distribution
sample_remaining <- function(longevity, mu, sigma, Q, n = length(longevity), seed = 1234) {
  
  set.seed(seed = seed)
  
  Fx <- pgengamma(longevity, mu = mu,
                  sigma = sigma, Q = Q)
  
  u <- runif(n)
  
  u_star <- Fx + u * (1 - Fx)
  
  Time <- qgengamma(u_star, mu = best_c[1],
                    sigma = best_c[2], Q = best_c[3]) |>
    floor()
  # t_remaining <- Time - longevity
  
  return(Time)
}

## persistence after the resampling (the MC step was already caried out, so the original persistence function cannot be used)
long_persistence <- function(sites) {
  
  tmp <- sites |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      mean   = mean(persistence, na.rm = TRUE),
      median = median(persistence, na.rm = TRUE),
      q025   = quantile(persistence, 0.025, na.rm = TRUE),
      q975   = quantile(persistence, 0.975, na.rm = TRUE),
      q25    = quantile(persistence, 0.25, na.rm = TRUE),
      q75    = quantile(persistence, 0.75, na.rm = TRUE),
      # normalised
      mean_n   = mean(norm_pers, na.rm = TRUE),
      median_n = median(norm_pers, na.rm = TRUE),
      q025_n   = quantile(norm_pers, 0.025, na.rm = TRUE),
      q975_n   = quantile(norm_pers, 0.975, na.rm = TRUE),
      q25_n    = quantile(norm_pers, 0.25, na.rm = TRUE),
      q75_n    = quantile(norm_pers, 0.75, na.rm = TRUE),
      # years (for plotting)
      start = median(start, na.rm = TRUE),
      end   = median(end, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(tmp)
}

# compute distance between small sites and large sites
dist_func <- function(df, size = 10) {
  
  phases <- sort(unique(df$start.y), decreasing = TRUE) # must be BP
  tmp <- lapply(
    X = phases, 
    FUN = function(phase) {
      sites <- df[df$start.y == phase,]
      large <- sites[sites$Area >= size,]
      small <- sites[sites$Area < size,]
      
      if (nrow(large) > 0) {
        
        small$dist <- as.numeric(apply(st_distance(small, large), 1, min))
        
        dmat <- st_distance(large, large)
        
        diag(dmat) <- NA  # remove self-distances
        
        
        large$dist <- as.numeric(apply(dmat, 1, min, na.rm = TRUE))
        
        # start of the phase with a city - end of life
        small$new_pers <- small$start.y - small$end.x
        # small$new_n_pers <- small$new_pers / (small$start.x - 510)
        
        large$new_pers <- large$persistence
        # large$new_n_pers <- large$norm_pers
        
        sites <- rbind(small, large)
        return(sites)
        
      } else {
        small$dist <- NA
        small$new_pers <- NA
        # small$new_n_pers <- NA
        
        return(small)
      } 
    }
  )
  tmp <- do.call(rbind, tmp)
  return(tmp)
  
}
