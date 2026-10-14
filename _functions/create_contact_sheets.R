library(magick)

create_contact_sheets <- function(
    original_dir,
    binary_dir,
    output_dir,
    thumb_size = 700){
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  originals <- list.files(
    original_dir,
    pattern = "\\.png$",
    full.names = TRUE
  )
  
  for(original_file in originals){
    
    image_name <- tools::file_path_sans_ext(
      basename(original_file)
    )
    
    cat("Creating sheet for",
        image_name,
        "\n")
    
    binary_files <- list.files(
      binary_dir,
      pattern = paste0("^", image_name),
      recursive = TRUE,
      full.names = TRUE
    )
    
    if(length(binary_files) == 0)
      next
    
    # Original image
    original <- image_read(original_file)
    
    original <- image_resize(
      original,
      paste0(thumb_size, "x", thumb_size, "!")
    )
    
    original <- image_annotate(
      original,
      "ORIGINAL",
      gravity = "north",
      size = 40,
      color = "white",
      boxcolor = "black"
    )
    
    panels <- list(original)
    
    # Binary images
    for(f in sort(binary_files)){
      
      img <- image_read(f)
      
      img <- image_resize(
        img,
        paste0(thumb_size, "x", thumb_size, "!")
      )
      
      label <- basename(f)
      
      label <- tools::file_path_sans_ext(label)
      
      label <- gsub(
        paste0("^", image_name, "_"),
        "",
        label
      )
      
      label <- gsub("_", "\n", label)
      
      img <- image_annotate(
        img,
        label,
        gravity = "north",
        size = 25,
        color = "white",
        boxcolor = "black"
      )
      
      panels[[length(panels) + 1]] <- img
      
    }
    
    n_panels <- length(panels)
    
    cols <- ceiling(sqrt(n_panels))
    rows <- ceiling(n_panels / cols)
    
    sheet <- image_montage(
      do.call(c, panels),
      tile = paste0(cols, "x", rows),
      geometry = paste0(
        thumb_size,
        "x",
        thumb_size,
        "+10+10"
      )
    )
    
    outfile <- file.path(
      output_dir,
      paste0(
        image_name,
        "_contact_sheet.png"
      )
    )
    
    image_write(
      sheet,
      path = outfile,
      format = "png"
    )
    
  }
}
