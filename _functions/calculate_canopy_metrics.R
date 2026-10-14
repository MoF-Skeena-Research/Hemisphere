img = "./canopy_photos/standardized_hemi/BB_FdMR2_hemi_masked.png"
options <- expand.grid(
  channel = c("2BG"),
  gamma = c(1.6),
  stretch = c(FALSE)
)
calculate_canopy_metrics <- function(
    img,
    options,
    lat = 54.8,
    tau = 0.6,
    uoc = 0.15
){
  
  cat("Processing", basename(img), "\n")
  
  results <- vector("list", nrow(options))
  
  for(i in seq_len(nrow(options))){
    
    opt <- options[i, ]
    
    result <- binarize_image(
      image_file = img,
      method     = "Otsu",
      channel    = opt$channel,
      gamma      = opt$gamma,
      stretch    = opt$stretch
    )
    
    gapfrac <- gapfrac_fisheye(
      result$binary,
      maxVZA  = 90,
      lens    = "equidistant",
      endVZA  = 70,
      display = FALSE,
      message = FALSE
    )
    
    canopy_metrics <- canopy_fisheye(gapfrac)
    
    par_image <- make_par_image(result$binary)
    
    par <- calculate_yearly_PAR(
      image = par_image,
      lat   = lat,
      tau   = tau,
      uoc   = uoc
    )
    
    par <- as.data.frame(par)
    
    par_summary <- tibble(
      PAR_DirectAbove  = sum(par$DirectAbove,  na.rm = TRUE),
      PAR_DiffuseAbove = sum(par$DiffuseAbove, na.rm = TRUE),
      PAR_DirectUnder  = sum(par$DirectUnder,  na.rm = TRUE),
      PAR_DiffuseUnder = sum(par$DiffuseUnder, na.rm = TRUE)
    ) |>
      mutate(
        
        PAR_TotalAbove =
          PAR_DirectAbove + PAR_DiffuseAbove,
        
        PAR_TotalUnder =
          PAR_DirectUnder + PAR_DiffuseUnder,
        
        PAR_Transmit =
          if_else(
            PAR_TotalAbove > 0,
            PAR_TotalUnder / PAR_TotalAbove,
            NA_real_
          ),
        
        Direct_Transmit =
          if_else(
            PAR_DirectAbove > 0,
            PAR_DirectUnder / PAR_DirectAbove,
            NA_real_
          ),
        
        Diffuse_Transmit =
          if_else(
            PAR_DiffuseAbove > 0,
            PAR_DiffuseUnder / PAR_DiffuseAbove,
            NA_real_
          ))
    
    results[[i]] <- bind_cols(
      canopy_metrics,
      par_summary
    ) |>
      mutate(
        Omega = if_else(L > 0, Le / L, NA_real_),
        
        Image    = result$image,
        Method   = "Otsu",
        Channel  = opt$channel,
        Gamma    = opt$gamma,
        Stretch  = opt$stretch,
        OptionID = i
      )
    
  }
  
  results <- bind_rows(results)
  
  results <- results |>
    mutate(
      plot_id = sub("^([^_]+_[^_]+).*", "\\1", id),
      site    = sub("_.*", "", id),
      canopy = sub("^[^_]+_([^_]+).*", "\\1", id)
    )
  
  results
}


# calculate_canopy_metrics <- function(img, options) {
#   
#   cat("Processing", basename(img), "\n")
#   
#   results <- vector("list", nrow(options))
#   
#   for(i in seq_len(nrow(options))) {
#     
#     opt <- options[i, ]
#     
#     result <- binarize_image(
#       image_file = img,
#       method     = "Otsu",
#       channel    = opt$channel,
#       gamma      = opt$gamma,
#       stretch    = opt$stretch
#     )
#     
#     gapfrac <- gapfrac_fisheye(
#       result$binary,
#       maxVZA = 90,
#       lens = "equidistant",
#       endVZA = 70,
#       display = FALSE,
#       message = FALSE
#     )
#     
#     canopy_metrics <- canopy_fisheye(gapfrac)
#     
#     par_image <- make_par_image(result$binary)
#     
#     par <- calculate_yearly_PAR(
#       image = par_image,
#       lat = 54.8,
#       tau = 0.6,
#       uoc = 0.15
#     )
#     
#     par_summary <- tibble(
#       PAR_DirectAbove  = sum(par[, "DirectAbove"]),
#       PAR_DiffuseAbove = sum(par[, "DiffuseAbove"]),
#       PAR_DirectUnder  = sum(par[, "DirectUnder"]),
#       PAR_DiffuseUnder = sum(par[, "DiffuseUnder"])
#     ) |>
#       mutate(
#         PAR_TotalAbove = PAR_DirectAbove + PAR_DiffuseAbove,
#         PAR_TotalUnder = PAR_DirectUnder + PAR_DiffuseUnder,
#         PAR_Transmit   = PAR_TotalUnder / PAR_TotalAbove
#       )
#     
#     results[[i]] <- bind_cols(
#       canopy_metrics,
#       par_summary
#     ) |>
#       mutate(
#         Image = result$image,
#         Method = "Otsu",
#         Channel = opt$channel,
#         Gamma = opt$gamma,
#         Stretch = opt$stretch,
#         OptionID = i
#       )
#   
#   }
#   
#   results <- bind_rows(results)
#   
#   results <- results |>
#     mutate(
#       plot_id = sub("^([^_]+_[^_]+).*", "\\1", id),
#       site    = sub("_.*", "", id),
#       canopy  = sub("^[^_]+_([^_]+).*", "\\1", id)
#     )
#   
#   results
# }