library(org.Hs.eg.db)
library(GO.db)
library(plyr)

mapGOToAccessions <- function() {
  frna2gostr <- read.table("/home/gilesc/Desktop/gammaseq/data/frna2go.tsv",sep="\t")
  v <- unlist(lapply(as.list(GO.db::GOTERM), Term))
  gostr2go <- names(v)
  names(gostr2go) <- v
  v <- NULL

  frna2gostr$GOID <- gostr2go[as.vector(frna2gostr$V2)]
  frna2gostr <- frna2gostr[which(!is.na(frna2gostr$GOID)),]
  go2frna <- dlply(frna2gostr, .(GOID), function(x) {as.vector(x$V1)})

  ## Concatenate frna-GO with gene-GO
  goterms <- union(names(as.list(org.Hs.egGO2EG)), names(go2frna))
  sapply(goterms, function(x) {c(unlist(org.Hs.egGO2EG[[x]]), unlist(go2frna[[x]]))})
}
