### Main function script ----
# load helper functions
source(file.path('scripts', 'helpers_sampling.R'))

## Rule-based sampling of calendar dates
sample_dates <- function(sites, site_id_col = "Sitio", area_col = "Area",
                         start_method = 'unif', end_method = 'unif', peak_method = 'unif',
                         peak_sd_perc = .1, perc = 5, terminus = 1) {
  
  # check for start/end mistakes. Start must happen before end
  invalid_dates <- sites[sites$start <= sites$end, ]
  if(nrow(invalid_dates) > 0) {
    stop(paste0('Error: start >= end. Are you using BC/AD years instead of BP?\n',
                'Sites affected: ', paste(unique(invalid_dates[[site_id_col]]), collapse = ', '), '\n'
    ))
  }
  # check for multiple measurement in the same phase 
  has_duplicates <- sapply(unique(sites[[site_id_col]]), function(site_id) {
    site_subset <- sites[sites[[site_id_col]] == site_id, ]
    nrow(site_subset) != length(unique(site_subset$start)) # getting a TRUE/FALSE
  })
  
  sites_with_duplicates <- unique(sites[[site_id_col]])[has_duplicates] # extract only TRUE rows
  
  if(length(sites_with_duplicates) > 0) {
    stop(paste0('Error: Duplicates found.\n',
                'Sites affected: ', paste(sites_with_duplicates, collapse = ', '), '\n'))
  }
  
  if (!start_method %in% c("unif", "exp")) stop("start_method must be 'unif' or 'exp'.")
  
  if (!end_method %in% c("unif", "exp")) stop("end_method must be 'unif' or 'exp'.")
  
  if (!peak_method %in% c("unif", "norm")) stop("peak_method must be 'unif' or 'norm'.")
  
  output <- list()
  
  for(id in unique(sites[[site_id_col]])) {
    tsite <- sites[sites[[site_id_col]] == id,]
    
    # step 1: sample start
    # cat(paste0('\nSampling site ', tsite[[site_id_col]][1]))
    
    if (start_method == "unif") {
      
      start_1 <- runif(1, min = tsite$end[1], max = tsite$start[1]) |> floor()
      
    } else if (start_method == "exp") {
      
      # define phase start/end
      a <- tsite$start[1]  # older BP
      b <- tsite$end[1]    # younger BP
      
      delay <- exp_sampling(a, b, perc = perc)
      
      start_1 <- a - delay
      
    } 
    
    # step 2&3: end and peaks 
    if(nrow(tsite) == 1) {
      
      # cat(paste0('\nThe site end during the same phase (single-phased).'))
      if(end_method == 'unif') {
        
        end_1 <- runif(1, min = tsite$end[1], max = start_1) |> floor()
        
        # cat('\nEnd year: \n', end_1, '\nIt was sampled between year ', start_1 + 1, ' and year ', tsite$end) 
      } else if (end_method == "exp") {
        
        delay <- exp_sampling(start_1, tsite$end[1], perc = perc)
        
        end_1 <- start_1 - delay
        
      }
      
      if(peak_method == 'unif') {
        peak_1 <- runif(
          n = 1,
          max = start_1,
          min = end_1
        ) |> floor() 
        
        peaks_all <- peak_1 # for consistency
        # cat(paste0('\nThe site is single-phased, only one peak is sampled and it is at year: \n'))
        # print(peak_1) 
        
      } else if(peak_method == 'norm') {
        
        peak_1 <- norm_peak(start = start_1, end = tsite$end[1], peak_sd_perc = peak_sd_perc)
          
        peaks_all <- peak_1 
        
      } 
      
      years_all <- c(start_1, peaks_all, end_1)
      
    } else {
      # sampling end for multi phase sites
      
      if(end_method == 'unif') {
        
        end_1 <- runif(
          n = 1, 
          max = tsite$start[nrow(tsite)], 
          min = tsite$end[nrow(tsite)]
        ) |> floor()
        
        # cat(paste0('Site ', tsite[[site_id_col]][1], ' has ', nrow(tsite), ' phases.'))
        # cat('\nSampling a the first and last peak.')
        
      } else if(end_method == 'exp') {
        
        a <- tsite$start[nrow(tsite)]
        b <- tsite$end[nrow(tsite)]
        
        delay <- exp_sampling(start = a, end = b, perc = perc)
        
        end_1 <- b + delay 
        
        # cat(paste0('Site ', tsite[[site_id_col]][1], ' has ', nrow(tsite), ' phases.'))
        # cat('\nSampling a the first and last peak.')
      } else stop("end_method must be 'unif' or 'exp'.")
      
      if(peak_method == 'unif') {
        # first peak
        peak_1 <- runif(
          n = 1, 
          max = start_1,
          min = tsite$end[1]
        ) |>
          floor()
        
        # cat('First peak: \n')
        # print(peak_1)
        
        # last peak
        peak_last <- runif(
          n = 1, 
          max = tsite$start[nrow(tsite)],
          min = end_1
        ) |> 
          floor()
        
        # cat('Last peak: \n')
        # print(peak_last)
        peaks_all <- c(peak_1, peak_last)
        
        if(nrow(tsite) > 2) {
          
          # cat('\nSampling the other peaks.')
          peaks_mid <- c()
          for(i in 2:(nrow(tsite)-1)) { # 2: and -1 because the first and last one we have already
            peak_mid <- runif(
              n = 1,
              max = tsite$start[i],
              min = tsite$end[i]
            ) |>
              floor() 
            
            # add it to the list
            peaks_mid <- c(peaks_mid, peak_mid)
            
            # cat(paste0('\nPeak ', i, ': \n'))
            # print(peak_mid)
            
          }
          # make a single list for all peaks
          peaks_all <- c(peak_1, peaks_mid, peak_last)
        } 
        
      } else if(peak_method == 'norm') {
        
        peak_1 <- norm_peak(start = start_1, end = tsite$end[1], peak_sd_perc = peak_sd_perc)
        
        # cat('First peak: \n')
        # print(peak_1)
        
        # the last phase lasts between its start and its sampled end (end_1)
        peak_last <- norm_peak(start = tsite$start[nrow(tsite)], end = end_1, peak_sd_perc = peak_sd_perc)

        # cat('Last peak: \n')
        # print(peak_last)
        peaks_all <- c(peak_1, peak_last)
        
        if(nrow(tsite) > 2) {
          
          # cat('\nSampling the other peaks.')
          peaks_mid <- c()
          
          for(i in 2:(nrow(tsite)-1)) { # 2: and -1 because the first and last one we have already
            
            peak_mid <- norm_peak(start = tsite$start[i], end = tsite$end[i], peak_sd_perc = peak_sd_perc)

            # add it to the list
            peaks_mid <- c(peaks_mid, peak_mid)
            
            # cat(paste0('\nPeak ', i, ': \n'))
            # print(peak_mid)
            
          }
          # make a single list for all peaks
          peaks_all <- c(peak_1, peaks_mid, peak_last)
        } 
        
      }
      
      years_all <- c(start_1, peaks_all, end_1)
    }  
    
    # Build area values to couple with years
    # 0.1 as min area (to start a site) and area for each phase. repeat the last for the end date
    area_values <- c(0.1, tsite[[area_col]], tsite[[area_col]][nrow(tsite)])
    # here we assume the site ends after the phase but we could add a NA/0 for sites that are no longer occupied by the last phase
    
    df <- data.frame(
      id    = id, 
      start = c(tsite$start[1], tsite$start, tsite$start[nrow(tsite)]),
      end   = c(tsite$end[1], tsite$end, tsite$end[nrow(tsite)]),
      year = years_all,
      area = area_values, 
      row.names = NULL
    )
    # if the site has its last peak before the a user defined phase, the area at t = end_1 should be 0 (abandonment)
    if(df$year[nrow(df)] > terminus) {
      df$area[nrow(df)] <- 0
    }
    output[[as.character(id)]] <- df
  }
  return(do.call(rbind, output))
}

### run examples 
# terminus indicates the start year of the last phase. If a site did not get to it then its last area would be forced to be 0
# t1 <- sample_dates(test, terminus = 510)
# t2 <- sample_dates(test, start_method = 'exp', end_method = 'exp', perc = 0.1, terminus = 510)
# ts2 <- t2[t2$id == 200,]

#### Growth rate function ---- 
# G = (P2 - P1) / (t2 - t1)
# we need to calculate the slope of the line connecting the two year estimates
growth_rate <- function(sampled_df) {
  
  required <- c("id", "year", "area")
  
  if (!all(required %in% colnames(sampled_df))) {
    stop("sampled_df columns must contain: id, year, area")
  }
  
  sampled_df$growth_cagr <- NA_real_
  sampled_df$growth_linear <- NA_real_
  
  for(id in unique(sampled_df$id)) {
    
    tsite <- sampled_df[sampled_df$id == id,]
    # make sure it is sorted (years BP)
    tsite <- tsite[order(tsite$year, decreasing = TRUE), ]
    
    prev_area <- c(NA, tsite$area[-nrow(tsite)])
    prev_year <- c(NA, tsite$year[-nrow(tsite)])
    
    year_diff <- pmax(1, prev_year - tsite$year) # make sure the min value is 1 year
    
    # sampled_df$prev_year[sampled_df$id == id] <- prev_year
    # sampled_df$year_diff[sampled_df$id == id] <- year_diff
    # sampled_df$prev_area[sampled_df$id == id] <- prev_area
    
    sampled_df$start[sampled_df$id == id]  <- tsite$start
    sampled_df$end[sampled_df$id == id]    <- tsite$end
    # CAGR: compound annual growth rate (percentage)
    growth_cagr <- (tsite$area / prev_area)^(1 / year_diff) - 1
    growth_cagr[is.infinite(growth_cagr)] <- NA
    sampled_df$growth_cagr[sampled_df$id == id] <- growth_cagr

    # Linear: simple annual growth rate (ha/year)
    growth_linear <- (tsite$area - prev_area) / year_diff
    growth_linear[is.infinite(growth_linear)] <- NA
    sampled_df$growth_linear[sampled_df$id == id] <- growth_linear
  }
  
  return(sampled_df)
}

run_mc <- function(sites, iterations = 1000, 
                   site_id_col = "Sitio", area_col = "Area",
                   start_method = 'unif', end_method = 'unif', perc = 5,
                   peak_method = 'unif', peak_sd_perc = .1, terminus = 1) {
  
  results_list <- list()  
  
  for(i in 1:iterations) { 
    
    if(i %% 100 == 0) cat(paste0("Iteration ", i, "/", iterations, "\n"))  # Optional progress
    
    tmp <- sample_dates(sites = sites, site_id_col = site_id_col,
                        start_method = start_method, end_method = end_method,
                        perc = perc, peak_method = peak_method, 
                        peak_sd_perc = peak_sd_perc, terminus = terminus) |> 
      growth_rate()
    
    tmp$iteration <- i
    results_list[[i]] <- tmp
  }
  
  # Combine and return results
  all_results <- do.call(rbind, results_list)
  return(all_results)
}
  
## run examples 
# tmc1 <- run_mc(test[test$Sitio %in% c(160,200),], terminus = 510, iterations = 1000)
# tmc2 <- run_mc(test[test$Sitio %in% c(160,200),], terminus = 510, iterations = 1000,
#                start_method = 'exp', end_method = 'exp', peak_method = 'norm')