# Load custom palettes from directory
load_custom_palettes <- function(dir) {
  txt_files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)

  if (length(txt_files) == 0) return(list())
  
  palettes <- lapply(txt_files, function(f) {
    colors <- readLines(f, warn = FALSE)
    colors <- trimws(colors[nzchar(colors)])
    palette_name <- tools::file_path_sans_ext(basename(f))
    list(name = palette_name, colors = colors)
  })
  
  # Return as named list of functions 
  named <- setNames(
    lapply(palettes, function(p) {
      force(p)
      function(n) p$colors[seq_len(min(n, length(p$colors)))]
    }),
    paste0(sapply(palettes, `[[`, "name"))
  )

  named
}