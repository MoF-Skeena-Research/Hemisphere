binarize_image <- function(image_file,
                           method = "Otsu",
                           channel = "2BG",
                           gamma = 2.2,
                           stretch = FALSE){
  image_name <-
    tools::file_path_sans_ext(
      basename(image_file)
    )
  
  outdir <- file.path(
    base_dir,
    "binarized",
    method
  )
  
  dir.create(outdir,
             recursive=TRUE,
             showWarnings=FALSE)
  
  fisheye <- import_fisheye(
    image_file,
    channel = channel,
    gamma = gamma,
    stretch = stretch,
    display = FALSE
  )
  
  binimage <- binarize_fisheye(
    fisheye,
    method=method,
    display=FALSE,
    export=FALSE
  )
  
  outfile <- file.path(
    outdir,
    paste0(
      image_name,
      "_",
      channel,
      "_g", gamma,
      "_s", stretch,
      ".png"
    )
  )
  
  png(
    outfile,
    width=1200,
    height=1200
  )
  
  plot(binimage)
  
  dev.off()
  
  list(
    image=image_name,
    method=method,
    fisheye=fisheye,
    binary=binimage,
    outfile=outfile
  )
  
}
