library(magick)

size <- 3500

# black background
mask <- image_blank(
  width = size,
  height = size,
  color = "black"
)

# draw white circle
mask <- image_draw(mask)

symbols(
  x = size/2,
  y = size/2,
  circles = size/2,
  inches = FALSE,
  bg = "white",
  fg = "white",
  add = TRUE
)

dev.off()

image_write(mask, "./canopy_photos/mask/hemisphere_mask_3500.png")
