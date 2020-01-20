pkgname <- "metaboCombineR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('metaboCombineR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("metaboCombineR-package")
### * metaboCombineR-package

flush(stderr()); flush(stdout())

### Name: metaboCombineR-package
### Title: A short title line describing what the package does
### Aliases: metaboCombineR-package metaboCombineR
### Keywords: package

### ** Examples

  ## Not run: 
##D      ## Optional simple examples of the most important functions
##D      ## These can be in \dontrun{} and \donttest{} blocks.   
##D   
## End(Not run)



cleanEx()
nameEx("runMetaboCombiner")
### * runMetaboCombiner

flush(stderr()); flush(stdout())

### Name: runMetaboCombiner
### Title: Run metaboCombineR algorithm
### Aliases: runMetaboCombiner

### ** Examples

myExample <- runMetaboCombiner()



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
