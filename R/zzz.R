
.onLoad <- function(libname, pkgname){
  library.dynam("ssKmeans", pkgname, libname)
  set.seed(NULL)
} # End of .onLoad()

.onUnload <- function(libpath){
  library.dynam.unload("ssKmeans", libpath)
} # End of .onUnload()

