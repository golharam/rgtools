library(RCurl)
library(XML)

ef   <- xmlTreeParse(getURL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=7157&retmode=xml"), useInternalNodes = T)
ns   <- getNodeSet(ef, "//Gene-commentary_accession")
accn <- sapply(ns, function(x) { xmlValue(x) } )
# get the NM_
accn[grep("NM_", unique(accn))]

