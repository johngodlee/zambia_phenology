# Package citations used in analysis 
# John Godlee (johngodlee@gmail.com)
# 2020-09-04

# List packages 
bibentries <- list(citation(), 
  citation("ade4"),
  citation("cluster"))

# Write to file
fileConn <- file("out/packages.bib")
writeLines(
  unlist(lapply(bibentries, toBibtex)), 
  fileConn)
close(fileConn)



