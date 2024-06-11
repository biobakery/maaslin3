#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

template <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    install.packages("devtools"); library("devtools");
  BiocManager::install(version = "3.10")
  BiocManager::valid()              # checks for out of date packages
  print("BiocManager has been succesfully installed.")
}

template()




