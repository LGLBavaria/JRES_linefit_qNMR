#Modified rnmrdata function, supplemented by the slot ‘DATA’ which contains the real numbers of the spectrum intensity, adjusted by the scaling factor nc.proc
library(rnmrfit)
#------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#'
#' Although it's likely to change in the near future, this class currently
#' serves to combine processed 1D data from the pdata folder of an experiment
#' with the acqus and procs parameter files.
#'
#' @slot processed 1r/1i data stored as a data.frame with direct.shift and
#'                 intensity as columns.
#' @slot acqus A list of acqus parameters.
#' @slot procs A list of procs parameters.
#' @slot DATA ADDED: contains the real numbers of the spectrum intensity, 
#'             adjusted by the scaling factor nc.proc
#' @slot parameters ADDED: contains additional aquisition parameters
#' @name NMRData1D-class
#' @export
NMRData1D <- setClass("NMRData1D", 
                      contains = "NMRData",
                      slots = c(DATA = 'numeric', 
                                parameters = 'list'))

# Validity testing consists of simply checking the processed data.frame columns
validNMRData1D <- function(object) {
  
  valid <- TRUE
  msg <- c()
  
  processed <- object@processed
  data <- object@DATA
  parameters <- object@parameters
  
  # Checking column names
  valid.columns <- c('direct.shift', 'intensity')
  new.msg <- sprintf('Processed data must have two columns: %s',
                     paste(valid.columns, collapse = ', '))
  
  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    msg <- c(msg, new.msg)
  }
  
  # Checking if DATA slot is a numeric vector
  if (!is.numeric(data) || length(data) == 0) {
    valid <- FALSE
    msg <- c(msg, "DATA slot must be a non-empty numeric vector")
  }
  
  # Checking if parameters slot is a numeric vector
  if (!is.list(parameters) || length(parameters) == 0) {
    valid <- FALSE
    msg <- c(msg, "parameters slot must be a non-empty numeric vector")
  }
  
  if ( valid ) TRUE
  else msg
}

setValidity("NMRData1D", validNMRData1D)

#------------------------------------------------------------------------
#' Constructor for generating an NMRData1D object
#'
#' Loads data in JCAMP file or Bruker directory and uses it to initialize an 
#' NMRData1D object.
#'
#' @param path Path to a Bruker scan directory or JCAMP file.
#' @param procs.number Specifies processing file number to use when loading
#'                     from a Bruker scan directory. Defaults to the smallest 
#'                     number in the pdata directory. Ignored if loading from
#'                     a JCAMP file.
#' @param blocks.number Specifies block number to use when loading from a
#'                      JCAMP file. Defaults to the first block encountered.
#'                      Ignored if loading from Bruker directory.  
#' @param ntuples.number Specifies ntuple entry number to use when loading 
#'                       from a JCAMP file. Defaults to the first ntuple 
#'                       entry encountered. Ignored if loading from Bruker 
#'                       directory.  
#' 
#' @return An NMRData1D object containing the 1r/1i processed data as well
#'         as the procs and acqus parameters.
#'
#' @export
read_1d_rnmrdata <- function(path, procs.number = NA, 
                       blocks.number = NA, ntuples.number = NA) {
  
  # If file exists at path, treating import as JCAMP, otherwise, Bruker scan
  if ( file.exists(path) && !dir.exists(path) ) {
    # Load jcamp
    jcamp <- read_jcamp(path, process.tags = TRUE, process.entries = TRUE)
    
    # If block number specified, check if it's valid
    if (! is.na(blocks.number) ) {
      # Double check that specified block exists
      msg <- sprintf("Specified block number not found in %s", path)
      if ( blocks.number > length(jcamp$blocks) ) stop(msg)
    }
    else {
      blocks.number <- 1
    }
    
    jcamp$blocks <- jcamp$blocks[blocks.number]
    
    # Check to make sure that ntuples exist
    msg <- 'Import from JCAMP file currently limited to NTUPLES entries'
    if (! 'ntuples' %in% names(jcamp$blocks[[1]]) ) stop(msg)
    
    # If ntuple number specified, check if it's valid
    if (! is.na(ntuples.number) ) {
      # Double check that specified block exists
      msg <- sprintf('Specified NTUPLES number not found in block %i of %s', 
                     blocks.number, path)
      if ( ntuples.number > length(jcamp$blocks[[1]]) ) stop(msg)
    }
    else {
      ntuples.number <- 1
    }
    
    jcamp$blocks[[1]]$ntuples <- jcamp$blocks[[1]]$ntuples[ntuples.number]
    
    # Flattening file
    jcamp.flat <- flatten_jcamp(jcamp)
    
    # Extracting processed data from ntuples
    descriptors <- jcamp.flat$ntuples$descriptors
    pages <- jcamp.flat$ntuples$pages
    
    # Checking variable names
    variables <- tolower(descriptors$var.name)
    real.index <- which(str_detect(variables, 'spectrum.*real'))
    imag.index <- which(str_detect(variables, 'spectrum.*imag'))
    
    # If both real and imaginary data isn't there, abort
    msg <- 'Import from JCAMP file currently limited to real/imaginary spectra'
    if ( length(c(real.index, imag.index)) < 2 ) stop(msg)
    
    # Double check that the first entry is frequency
    msg <- 'Import from JCAMP file currently limited to frequency abscissa'
    if (! str_detect(variables[1], 'freq') ) stop(msg)
    
    # Proceed to extract data
    real.data <- pages[[real.index - 1]]
    imag.data <- pages[[imag.index - 1]]
    
    # Checking that frequency is the same
    real.frequency <- real.data[, 1]
    imag.frequency <- imag.data[, 1]
    msg <- 'Mismatch in real and imaginary frequency, likely parsing error'
    if ( any(abs(real.frequency - imag.frequency) > 1e-6) ) stop(msg)
    
    # Starting with raw values
    frequency <- real.frequency
    real.data <- real.data[, 2]
    imag.data <- imag.data[, 2]
    
    # Scaling if factors is non zero
    scale <- descriptors$factor[1]
    if (scale > 1e-6) frequency <- frequency*scale
    
    scale <- descriptors$factor[real.index]
    if (scale > 1e-6) real.data <- real.data*scale
    
    scale <- descriptors$factor[imag.index]
    if (scale > 1e-6) imag.data <- imag.data*scale
    
    # Offsetting if max-min difference is non zero
    d.max <- descriptors$max[1]
    d.min <- descriptors$min[1]
    if ( (d.max - d.min) > 1e-6 ) {
      frequency <- frequency - max(frequency) + d.max
    }
    
    d.max <- descriptors$max[real.index]
    d.min <- descriptors$min[real.index]
    if ( (d.max - d.min) > 1e-6 ) {
      real.data <- real.data - max(real.data) + d.max
    }
    
    d.max <- descriptors$max[imag.index]
    d.min <- descriptors$min[imag.index]
    if ( (d.max - d.min) > 1e-6 ) {
      imag.data <- imag.data - max(imag.data) + d.max
    }
    
    # Doing one final check on the direct shift to check on offset
    direct.shift <- frequency/jcamp.flat$sf
    
    delta <- jcamp.flat$offset - max(direct.shift)
    if ( abs(delta) > 1e-6 ) direct.shift <- direct.shift + delta
    
    # Finally, combine the data
    intensity <- complex(real = real.data, imaginary = imag.data)
    processed <- data.frame(direct.shift = direct.shift,
                            intensity = intensity)
  
    # Returning class object
    new("NMRData1D", processed = processed, parameters = jcamp.flat,
        procs = list(), acqus = list(), DATA = NULL)       
  }
  else {
    # First, loading procs parameters
    procs <- read_procs_1d(path, procs.number)
    
    # Using the procs file to load the processed data
    processed <- read_processed_1d(path, procs, procs.number)
    
    # Finally, loading the general acquisition parameters
    acqus <- read_acqus_1d(path)
    
    #HINZUGEFUEGT:Einbezug des Bruker Skalierungsparameters
    DATA <- Re(processed$intensity*(2^procs$nc.proc)) 
    
    #HINZUGEFUEGT: parameters Liste a la bbio_spec aus MATLAB
    parameters <- list(TITLE=sub(".*/(W[0-9]+/\\d+).*", "\\1", acqus[["comments"]][["8"]]),
                       ACQUS=acqus,
                       PROCS=procs,
                       Date=sub("^(\\d{4}-\\d{2}-\\d{2}).*", "\\1", acqus[["comments"]][["7"]]),
                       NC.PROC=acqus$nc.proc,
                       BYTORDP=procs$bytordp,
                       LB=procs$lb,
                       PHC0=procs$phc0,
                       PHC1=procs$phc1,
                       SF=procs$sf,
                       file=path,
                       Data=Re(processed$intensity*(2^procs$nc.proc)),
                       Imag=Im(processed$intensity*(2^procs$nc.proc)),
                       maxppm=procs$offset,
                       minppm=procs$offset-(procs$sw.p/procs$sf),
                       sw=procs$sw.p/procs$sf,
                       size=procs$si)
    
    # Returning class object
    new("NMRData1D", processed = processed, parameters = parameters,
        procs = procs, acqus = acqus, DATA = DATA)
  }
}
