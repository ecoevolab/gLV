#' Generate a Unique ID for Output Files
#'
#' This function generates a unique ID by checking the output directory to ensure no duplicate IDs are created.
#'
#' @param wd Character: The working directory path where output files are stored. The function will check the "Outputs" subdirectory to generate a unique ID based on the existing files.
#'
#' @details The function creates a unique ID by ensuring that the newly generated ID does not match any IDs in the output files within the specified working directory.
#'
#' @return A character string containing a unique ID to be used for saving files.
#'
#' @importFrom ids random_id
#'
#' @examples
#' # Example usage:
#' wd <- "~/Documents/LAB_ECO"
#' uniqueID <- forge_ID(wd)
#' print(uniqueID)


forge_ID <- function(wd) {

  # Ensure the ids package is available
  if (!requireNamespace("ids", quietly = TRUE)) {
    stop("The 'ids' package is required but not installed.")
  }

  repeat {
    directory = file.path(wd,"Outputs")
    ID <- ids::random_id(1, 3)
    all_files <- list.files(directory, full.names = TRUE) # List all files in the directory
    matching_files <- grep(ID, basename(all_files), value = TRUE) # Check if any file name contains the ID as a substring

    # Check if any file matches the pattern
    if (length(matching_files) == 0) {
      return(ID = ID) # Return the unique ID
    }
  }
}
