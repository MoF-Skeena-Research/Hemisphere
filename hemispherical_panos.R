library(magick)
library(hemispheR)
library(dplyr)
library("exifr")
require(tictoc)
  

    ### Load necessary libraries:
    library(tidyverse) # For data manipulation
    
    # ImageMagick is a command line tool for image manipulation. We will call ImageMagick from within R using the magick package, but first you'll need to download [https://imagemagick.org/script/download.php] and install ImageMagick on your machine.
    library(magick) # For image manipulation
    # Check to ensure that ImageMagick is installed.
    magick_config()$version 
    
    # ImageR also requires ImageMagick
    library(imager) # For image display
    
    library(exifr) # For extracting metadata
    
    # For binarizing and calculating some canopy metrics, we will use Chiannuci's hemispheR package, which we need to install from the development version.
    library(devtools)
    devtools::install_git("https://gitlab.com/fchianucci/hemispheR")
    library(hemispheR) # For binarization and estimating canopy measures
    
    ### Run the processing loop
    
    # Define the input path to the directory with panos
    #drive = "C:/Users/myarham/OneDrive - Government of BC/"
    drive = "D:/OneDrive - Government of BC/"
    focal_path <- paste0(drive, "OffSite-Trials/Michelle/Lw_Understory/Insta360_Canopy_photos/test folder") 
    list_of_panos <- list.files(focal_path) # Get the list of all panos
    
    # Define the output directory path
    output <- paste0(drive, "OffSite-Trials/Michelle/Lw_Understory/Insta360_Canopy_photos/Test Results") # Instantiate an empty object to receive the results
    
    # Loop through each pano in the list
    for(i in 1:length(list_of_panos)) {
     # if (i == 1) {
        T0 <- Sys.time() # Used for estimating remaining time
     
      

      
      ###
      
      T1 <- Sys.time() # Used for time check
      
      # Construct the full path to the current image
      focal_image_path <- paste0(focal_path, "/", list_of_panos[i])
      
      # Extract the base name (without extension) of the current image
      focal_image_name <- sub("\\.[^.]+$", "", basename(focal_image_path))
      
      # Print progress for debugging
      print(paste("Processing image:", focal_image_name))
      
      # Add additional processing steps here...
    
    
    # Calculate and print total processing time
    T_total <- Sys.time() - T0
    print(paste("Total processing time:", T_total))
    
      # You can choose which variables you'd like to retain
      xmp_data <- 
        read_exif(focal_image_path) %>%
        select(
          SourceFile,
          Make,
          Model,
          Megapixels,
          GPSLatitude,
          GPSLongitude#,
          #GPSAltitude,
        )
      
      # The first step in the process is to convert the equirectangular image from our phone into a hemispherical image.
      
      ### Convert the equirectangular image to hemisphere
      
      pano <- image_read(focal_image_path)
      
      # Store the pano width to use in scaling and cropping the image
      pano_width <- image_info(pano)$width
      
      # Store the pano heading in order to rotate the hermispherical image to standardize true north as the top of the image. This only matters for analyses like global site factor or through-canopy radiation that require plotting a sunpath over the hemisphere.
      image_heading <- read_exif(focal_image_path)$PoseHeadingDegrees
      
      #edit
      image_heading <- 0
      #
      
      
      # To process the image, we need to scale it, reproject it into polar coordinates, reorient it, and rotate it to true north.
      pano_hemisphere <- pano %>%
        # Crop to retain the upper hemisphere
        image_crop(geometry_size_percent(100, 50)) %>%
        # Rescale into a square to keep correct scale when projecting in to polar coordinate space
        image_resize(geometry_size_percent(100, 400)) %>%
        # Remap the pixels into polar projection
        image_distort("Polar",
                      c(0),
                      bestfit = TRUE) %>%
        image_flip() %>%
        # Rotate the image to orient true north to the top of the image
        image_rotate(image_heading) %>%
        # Rotating expands the canvas, so we crop back to the dimensions of the hemisphere's diameter
        image_crop(paste0(pano_width, "x", pano_width, "-", pano_width/2, "-", pano_width/2))
      
      # Plot the hemispherical image. The image looks funny because the outer pixels are extended ny interpolation and we've rotated the image. Most analyses define a bounding perimeter to exclude any pixels outside of the circular hemisphere, so the weird border shouldn't matter. But, we can add a black mask to make the images look better.
      pano_hemisphere
      
      ### Create black mask for the image (this isn't really neccessary, but makes the images look nicer)
      # Get the image mask vector file
      image_mask <- image_read(paste0(drive, "OffSite-Trials/Michelle/Lw_Understory/Insta360_Canopy_photos/HemiPhotoMask.svg")) %>%
        image_transparent("white") %>%
        image_resize(geometry_size_pixels(width = pano_width, height = pano_width)) %>%
        image_convert("png")
      
      masked_hemisphere <- image_mosaic(c(pano_hemisphere, image_mask))
      setwd(paste0(drive, "OffSite-Trials/Michelle/Lw_Understory/Insta360_Canopy_photos"))
      
      getwd()
      # We'll store the masked hemispheres in their own subdirectory.
      if(dir.exists("./masked_hemispheres/") == FALSE){
        dir.create("./masked_hemispheres/")
      } # If the subdirectory doesn't exist, we create it.
      getwd()
      
      masked_hemisphere_path <- paste0("./masked_hemispheres/", focal_image_name, "hemi_masked.jpg") # Set the filepath for the new image
      
      image_write(masked_hemisphere, masked_hemisphere_path) # Save the masked hemispherical image
      
      # At this point, you can process the hemispherical images however you'd like to calculate canopy and light metrics.
      
      ### For this example, I'm going to use Chiannuci's hemispheR package in order to keep this entire pipeline in R.
      # The next step is to import the image. hemispheR allows for lots of fine-tuning. Check out the docs to learn what all of the options are. These settings most closely replicate the processing I used in my 2021 paper. 
      tic()
      fisheye <- import_fisheye(masked_hemisphere_path,
                                channel = '2BG',
                                circ.mask = list(xc = pano_width/2, yc = pano_width/2, rc = pano_width/2),
                                gamma = 2.2,
                                stretch = FALSE,
                                display = TRUE,
                                message = TRUE)
      toc()
      plot(fisheye)
      # Now, we need to binarize the images, converting all sky pizels to white and everything else to black (ideally). Again, there are lots of optionas available in hemispheR. You can decides which settings are right for you. However, I would suggest keeping zonal set to FALSE. Because spherical panoramas are exposing each of the 36 images separately, there is no need to use zonal FIX THIS SENTENCE.
      # I also suggest keeping export set to TRUE so that the binarized images will be saved into a subdirectory named 'results'.
      tic()
      binimage <- binarize_fisheye(fisheye,
                                   method = 'Otsu',
                                   # We do NOT want to use zonal threshold estimation since this is done by the camera
                                   zonal = FALSE,
                                   manual = NULL,
                                   display = TRUE,
                                   export = TRUE)
      toc()
      plot(binimage)
      # Unfortunately, hemispheR does not allow for estimation of understory light metrics. If you need light estimates, you'll have to take the binarized images and follow my instructions for implementing Gap Light Analyzer.
      
      ### Estimate canopy metrics
      # Assuming all you need is canopy metrics, we can continue with hemispheR and finalize the whole pipeline in R.
      
  tic()
      gapfrac <- gapfrac_fisheye(
        binimage,
        maxVZA = 90,
        # Spherical panoramas are equidistant perforce
        lens = "equidistant",

        # startVZA = 0,
         #endVZA = 70,
        # nrings = 5,
        # nseg = 8,
        display = FALSE,
        message = FALSE

      )
    toc() 
      # Note that the 'x' column here is the gap fraction estimate.
      canopy_report <- canopy_fisheye(
        gapfrac
      )
      
      output_report <- 
        xmp_data %>%
        bind_cols(
          canopy_report
        ) %>%
        rename(
          GF = x,
          HemiFile = id
        )
      glimpse(output_report)
      output <- as.data.frame(output)
      output <-
        bind_rows(output,
                  output_report)
      
      T2 <- Sys.time()
      T_instance <- difftime(T2, T1, units = "secs")
      T_total <- difftime(T2, T0, units = "secs")
      T_average <- T_total/i
      
      print(paste0("Completed ", i, " of ", length(list_of_panos), " images in ", round(T_instance, 0), " seconds."))
      print(paste0("Estimated ", round(((length(list_of_panos) - i) * T_average)/60, 1), " minutes remaining."))
      }
write.csv(output_report, "./canopy_output.csv", row.names = FALSE)
