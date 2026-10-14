
# image_file = files_360,
# 
# output_dir = "./canopy_photos/converted_360",
# mask_path = mask

process_360_image <- function(
    image_file,
    output_dir,
    mask_path
) {
  
  image_mask <- image_read(mask_path) |>
    image_transparent("white") |>
    image_convert("png")
  
  pano <- image_read(image_file)
  
  pano_width <- image_info(pano)$width
  
  hemi <-
    pano |>
    image_crop(geometry_size_percent(100,50)) |>
    image_resize(geometry_size_percent(100,400)) |>
    image_distort("Polar", c(0), bestfit = TRUE) |>
    image_flip() |>
    image_crop(
      paste0(
        pano_width,
        "x",
        pano_width,
        "-",
        pano_width/2,
        "-",
        pano_width/2
      )
    )
  
  resized_mask <-
    image_resize(
      image_mask,
      geometry_size_pixels(
        pano_width,
        pano_width
      )
    )
  
  hemi <- image_mosaic(c(hemi, resized_mask))
  
  outfile <- file.path(
    output_dir,
    paste0(
      tools::file_path_sans_ext(
        basename(image_file)
      ),
      "_hemi.jpg"
    )
  )
  
  image_write(hemi, outfile)
  
  return(outfile)
}

# process_360_image <- function(image_file,
#                               output_dir,
#                               image_mask) {
#   
#   pano <- image_read(image_file)
#   
#   pano_width <- image_info(pano)$width
#   
#   hemi <-
#     pano |>
#     image_crop(geometry_size_percent(100,50)) |>
#     image_resize(geometry_size_percent(100,400)) |>
#     image_distort("Polar", c(0), bestfit = TRUE) |>
#     image_flip() |>
#     image_crop(
#       paste0(
#         pano_width,
#         "x",
#         pano_width,
#         "-",
#         pano_width/2,
#         "-",
#         pano_width/2
#       )
#     )
#   
#   resized_mask <-
#     image_resize(
#       image_mask,
#       geometry_size_pixels(
#         pano_width,
#         pano_width
#       )
#     )
#   
#   hemi <- image_mosaic(c(hemi,resized_mask))
#   
#   outfile <- file.path(
#     output_dir,
#     paste0(
#       tools::file_path_sans_ext(
#         basename(image_file)
#       ),
#       "_hemi.jpg"
#     )
#   )
#   
#   image_write(hemi,outfile)
#   
#   outfile
# }
