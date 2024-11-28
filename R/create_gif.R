#' Create a gif of blok3D
#'
#' @description  Functions creates a gif of rotating polygons displayed in a currently actively rgl device produced with blok3D
#' @param output_file name and path of output file
#' @param frames integer. Number of frames.
#' @param angle_per_frame numeric. Angle to rotate per frame in degrees.
#' @param fps integer. Frames per second.
#' @return NULL
#' @author Liudas Daumantas
#' @family HespDiv visualization options
#' @importFrom rgl par3d rotationMatrix view3d rotate3d rgl.snapshot
#' @importFrom magick image_read image_animate image_write
#' @notes You can adjust the size of rgl device window to control the size of gif.
#' @export
create_gif <- function(output_file = "rotating_polygons.gif",
                                             frames = 90, angle_per_frame = 5,
                       fps = 10) {

  angles <- rep(angle_per_frame * pi / 180, frames)
  # Directory for temporary files
  temp_dir <- tempdir()
  file_list <- character(frames)  # Store file paths for each frame
  rgl::par3d(userMatrix = rgl::rotationMatrix(-45 * pi / 180, 45, 0, 1))
  # Render each frame with rotation
  for (i in seq_len(frames)) {

    rgl::view3d(userMatrix = rgl::rotate3d(rgl::par3d("userMatrix"),
                                 angles[i], 0, 0, 1))

    # Save the frame as a PNG
    file_list[i] <- file.path(temp_dir, sprintf("frame_%03d.png", i))
    rgl::rgl.snapshot(file_list[i])
  }

  # Combine frames into a GIF
  images <- magick::image_read(file_list)  # Read all frames into magick
  animation <- magick::image_animate(images, fps = fps)  # Create animation
  magick::image_write(animation, path = output_file)  # Save as GIF

  # Clean up temporary files
  unlink(file_list)

  message("GIF saved to: ", output_file)
}
