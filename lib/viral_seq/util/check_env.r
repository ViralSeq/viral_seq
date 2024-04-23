dir.create(
  path = Sys.getenv("R_LIBS_USER"),
  showWarnings = FALSE,
  recursive = TRUE
)

packages <- c("ggplot2", "phangorn", "ape", "scales", "ggforce", "cowplot", "magrittr", "gridExtra")

install.packages(setdiff(packages, rownames(installed.packages())), lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
