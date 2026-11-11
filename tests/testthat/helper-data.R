skip_no_data <- function() testthat::skip_if_not_installed("ss3diags")
load_ss3 <- function(name){
  e <- new.env(); utils::data(list = name, package = "ss3diags", envir = e); e
}
