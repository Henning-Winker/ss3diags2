# Render the user guide, write the commented R script and update the GitHub pages copy.
# Run from the package root: source("inst/guide/updguide.R")

# 1. R script of the guide: prose as comments, code as is (chunks with eval = FALSE inside if(FALSE))
guide2script <- function(dir = "inst/guide/chapters", out = "inst/ss3diags2_user_guide.R", width = 100){
  wrap <- function(x, prefix = "# "){
    w <- strwrap(x, width = width - nchar(prefix) - 2)
    paste0(prefix, c(w[1], if(length(w) > 1) paste0(if(grepl("^\\* ", w[1])) "  " else "", w[-1])))
  }
  res <- c("# ss3diags2 user guide: R code with comments",
           "# Generated from inst/guide/chapters by inst/guide/updguide.R; edit the chapters, not this file.",
           "# Rendered guide: https://henning-winker.github.io/ss3diags2/", "",
           "library(ss3diags2)", "library(ggplot2)", "")
  for(f in sort(list.files(dir, pattern = "\\.Rmd$", full.names = TRUE))){
    if(grepl("references", f)) next
    l <- readLines(f, warn = FALSE)
    i <- 1; para <- character(0)
    flush <- function(){ if(length(para)) res <<- c(res, wrap(paste(para, collapse = " ")), ""); para <<- character(0) }
    while(i <= length(l)){
      x <- l[i]
      if(grepl("^```\\{r", x)){                          # code chunk
        flush()
        noeval <- grepl("eval *= *FALSE", x)
        j <- i + 1; code <- character(0)
        while(j <= length(l) && !grepl("^```\\s*$", l[j])){ code <- c(code, l[j]); j <- j + 1 }
        res <- c(res, if(noeval) c("if(FALSE){", paste0("  ", code), "}") else code, "")
        i <- j + 1; next
      }
      if(grepl("^#+ ", x)){                              # headings
        flush(); lev <- nchar(sub(" .*", "", x)); h <- sub("^#+ +", "", sub(" *\\{#.*\\}$", "", x))
        res <- c(res, if(lev == 1) c(strrep("#", width), paste("#", toupper(h)), strrep("#", width)) else
          paste0("# ---- ", h, " ", strrep("-", max(0, width - 12 - nchar(h)))), "")
      } else if(grepl("^\\|", x)){                       # tables
        flush(); if(!grepl("^\\|[-| ]+\\|$", x)) res <- c(res, paste("#", gsub("`", "", x)))
        if(i == length(l) || !grepl("^\\|", l[i + 1])) res <- c(res, "")
      } else if(!nzchar(trimws(x))) flush()
      else {
        if(grepl("^(- |> \\*\\*|[0-9]+\\. )", x)) flush()  # new bullet, numbered item or note
        para <- c(para, gsub("\\*\\*|`", "", sub("^> ?", "", sub("^- ", "* ", x))))
      }
      i <- i + 1
    }
    flush()
  }
  writeLines(res, out)
  invisible(out)
}
guide2script()

# 2. HTML guide
rmarkdown::render("inst/guide/ss3diags2_user_guide.Rmd")

# 3. GitHub pages (main branch: index.html at the root and docs/index.html)
dir.create("docs", showWarnings = FALSE)
file.copy("inst/guide/ss3diags2_user_guide.html", "index.html", overwrite = TRUE)
file.copy("inst/guide/ss3diags2_user_guide.html", "docs/index.html", overwrite = TRUE)
