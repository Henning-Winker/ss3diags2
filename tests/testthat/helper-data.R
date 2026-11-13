skip_no_data <- function() invisible(TRUE)
# cookbook models ship with ss3diags2; other ss3diags examples (gob.her) need ss3diags
load_ss3 <- function(name){
  e <- new.env()
  pk <- c(natl.sma = "ss3sma", pac.hke = "ss3phk")
  if(name %in% names(pk)){
    nm <- pk[[name]]
    utils::data(list = c(nm, sub("ss3", "retro.", nm), sub("ss3", "aspm.", nm)), package = "ss3diags2", envir = e)
  } else {
    testthat::skip_if_not_installed("ss3diags")
    utils::data(list = name, package = "ss3diags", envir = e)
  }
  e
}
