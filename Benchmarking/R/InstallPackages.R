local({r <- getOption('repos')                                                                                                                                
r['CRAN'] <- 'http://www.stats.bris.ac.uk/R/'                                                                                                          
options(repos=r)                    
}) 

if(!require(tufte)) {
  install.packages("tufte", dependencies=T)
}