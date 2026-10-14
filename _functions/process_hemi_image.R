#image_file = files_hemi[1]

process_hemi_image <- function(image_file,
                               output_dir,
                               image_mask){
  
  # image_mask <- image_read(mask_path) |>
  #   image_transparent("white") |>
  #   image_convert("png")
  # 
  hemi <- image_read(image_file)
  
  hemi_width <- image_info(hemi)$width
  
  hemi <- image_resize(
    hemi,
    "3000x3000!"
  )
  # resized_mask <- image_resize(
  #   image_mask,
  #   geometry_size_pixels(
  #     hemi_width,
  #     hemi_width
  #   )
  # )
  
  # hemi <- image_mosaic(c(hemi, resized_mask))
  
  outfile <- file.path(
    output_dir,
    paste0(
      tools::file_path_sans_ext(
        basename(image_file)
      ),
      "_masked.png"
    )
  )
  
  image_write(hemi, outfile, format = "png")
  
  cat("Writing to:", outfile, "\n")
  outfile
}
