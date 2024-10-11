#' Read Parameters by Line
#'
#' This function reads the original parameters used for simulation. 
#' These parameters were created with the \link{forge_data} function and saved with the \link{params_line_saver} function.
#'
#' @param uniqueID Character. A unique ID string used for saving files. This ID should be generated previously using the \code{generate_uniqueID} function.
#' @param wd Character. The path to the working directory where files will be saved. The file will be saved at: \code{wd/Parameters/Seeds_save.tsv}.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{Interactions}: A matrix representing the interaction values between species.
#'   \item \code{Growths}: A numeric vector representing the growth rates of the simulated species.
#'   \item \code{Population}: A numeric vector representing the initial abundances of the simulated species.
#'   \item \code{Seeds}: A numeric vector representing the seeds used to generate the populations, interaction matrix, and growth rates.
#' }
#'
#' @import readr 
#' 
#' @examples
#' wd <- "~/Documents/LAB_ECO/Simulations"
#' 
#' # Generate parameters
#' # forge_seeds(n = 200, min = 2, max = 1000, wd)
#' 
#' seeds_path <- file.path(wd, "Seeds.tsv")
#' params <- init_data(N_species = 10, seeds_path, C0 = 0.45, CN = 0.2, Diag_val = -0.5)
#' 
#' # Generate unique ID
#' uniqueID <- forge_id(wd)
#' 
#' # Save parameters by line
#' params_line_saver(params, uniqueID, wd)
#' 
#' # Read original parameters with the function 
#' res <- params_line_reader(uniqueID, wd)
#' 
#' @export


params_line_reader <- function(uniqueID, wd){

  # Ensure the readr package is available
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The 'readr' package is required but not installed.")
  }
  
  
  # Generate path
  path <- file.path(wd, paste0("Parameters/P_",uniqueID,".tsv",sep = "") )
  
  # Leer tsv como un vector de caracteres
  tsv_lines <- readLines(path)
  
  # Obtener lineas donde comienzan mis datos
  Interactions_line <- grep("Interactions", tsv_lines)
  Grows_line <- grep("Grow rates", tsv_lines)
  Population_line <- grep("Initial population", tsv_lines)
  Seed_line <- grep("Seeds", tsv_lines)
  
  # Obtener tabla de interacciones
  tmp <- tsv_lines[(Interactions_line + 1):Grows_line - 2 ]
  Interacs <- read.table(text = tmp, sep = "\t", skip = 2, fill = TRUE)
  
  # Obtener tabla de grows rates
  tmp <- tsv_lines[(Grows_line + 1):Population_line - 2]
  Grows <- read.table(text = tmp, sep = "\t", skip = 2, fill = TRUE)
  
  # Obtener tablas de poblaciones iniciales
  tmp <- tsv_lines[(Population_line + 1):Seed_line - 2]
  Pobl <- read.table(text = tmp, sep = "\t", skip = 2, fill = TRUE)
  
  # Obtener tablas de semillas
  tmp <- tsv_lines[Seed_line: (Seed_line + 4)]
  Seed <- read.table(text = tmp, sep = "\t", skip = 1, fill = TRUE)
  
  # Return parameters as a list
  result <- list(Population = Pobl,
                 Interactions = as.matrix(Interacs),
                 Growths = Grows,
                 Seeds = Seed)
  return(result)
  
}



