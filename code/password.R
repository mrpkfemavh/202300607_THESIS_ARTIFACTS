#' Password Strength
#'
#' Test the strength of a password
#'
#' @param x a character string of the password to test
#' @return a list of two elements
#'  - ENTROPY; a numeric value for the password entropy
#'  - STRENGTH; a strength rating of the password (Very Weak, Weak, Moderate, Strong, or Very Strong)
#'
#' @export

password_strength <- function(x)
{
  CHARSET <- as.numeric(charset_size(x))
  
  ENTROPY <- log2(CHARSET ^ nchar(x))
  
  if (ENTROPY < 28) {
    STRENGTH <- 'Very Weak'
  }
  if (ENTROPY >= 28 & ENTROPY <= 35) {
    STRENGTH <- 'Weak'
  }
  if (ENTROPY >= 36 & ENTROPY <= 59) {
    STRENGTH <- 'Moderate'
  }
  if (ENTROPY >= 60 & ENTROPY <= 127) {
    STRENGTH <- 'Strong'
  }
  if (ENTROPY >= 128) {
    STRENGTH <- 'Very Strong'
  }
  
  return(list(ENTROPY = ENTROPY, STRENGTH = STRENGTH))
}

#' Character Set Size
#'
#' Determine the Charset size based on what combination of characters has been used in the password
#'
#' @param x a character string of the password to test
#' @return a numeric value of the charset size
#'
#' @export

charset_size <- function(x)
{
  CHARSET <- c(0, 0, 0, 0)
  
  if (isTRUE(stringr::str_detect(x , '[:punct:]'))) {
    CHARSET[1] <- 22
  }
  
  if (isTRUE(stringr::str_detect(x , '[:upper:]'))) {
    CHARSET[2] <- 4
  }
  
  if (isTRUE(stringr::str_detect(x , '[:lower:]'))) {
    CHARSET[3] <- 26
  }
  
  if (isTRUE(stringr::str_detect(x , '[:digit:]'))) {
    CHARSET[4] <- 10
  }
  
  CHARSET_SIZE <- sum(CHARSET)
  
  return(CHARSET_SIZE)
}

password_strength("AAATTTTAATT")[[1]]



neighbourhood <- "ACTGGTTTCGGAT"
neighbourhood_complexity <- function(neighbourhood) {
  symbols_possible <- "ACTG"
  symbol_set <- sum(c(0, 4, 0, 0)) # ACTG  add others if more
  symbol_set_size <- symbol_set
  log2(symbol_set_size ^ nchar(neighbourhood))
}
neighbourhood_complexity(neighbourhood)



neighbourhood <- "AATTGGTTCCGGTTAA"



