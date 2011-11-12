library(ROCR)
data(coexpressed.Hs)
data(ontology.Hs)

noncoding.transcripts <- names(coexpressed.Hs)[substr(names(coexpressed.Hs), 1, 2)=="FR"]
coding.transcripts <- setdiff(names(coexpressed.Hs), noncoding.transcripts)

coexpression.cpt <- function(ontology.id, laplaceSmoothing=T) {
  "Probability of annotation conditioned on the number of coexpressed genes also annotated."
  cox <- coexpressed.Hs
  ids <- unlist(ontology.Hs[ontology.id])
  counts <- lapply(cox,function(v) length(intersect(ids,v[2:length(v)])))
  y <- as.numeric(names(cox) %in% ids)
  result <- table(unlist(counts),y)
  if (laplaceSmoothing)
    result[result==0] <- 1
  result[,"1"] / apply(result,1,sum)
}

cpt.lmfit <- function(ontology.id) {
  y <- coexpression.cpt(ontology.id)
  x <- as.numeric(names(cpt))
  fit <- lm(y ~ x)
  structure(coef(fit), names=c("b","m"))
}

cpt.lmfits <- function() {
  lengths <- lapply(ontology.Hs, length)
  ids <- names(lengths[lengths >= 5 & lengths <= 500])
  result <- as.data.frame(t(as.data.frame(mclapply(ids, cpt.lmfit))))
  result$n <- unlist(lapply(ontology.Hs[ids], length))
  rownames(result) <- ids
  result
}

predict.functions <- function(transcript.id) {
  
}

predict.ontology.performance <- function(ontology.id,
                                         model=coexpression.cpt(ontology.id),
                                         transcripts=names(coexpressed.Hs)) {
  labeled.transcripts <- intersect(transcripts, unlist(ontology.Hs[ontology.id]))
  if (length(labeled.transcripts)==0)
    stop("Gold standard has no positive training examples")
  x <- lapply(coexpressed.Hs[transcripts], function(v) {
    length(intersect(labeled.transcripts, v[2:length(v)]))})
  prediction(model[as.character(x)],
             as.numeric(transcripts %in% labeled.transcripts))
}

contrast.coding.noncoding.ontology.performance <- function(ontology.id) {
  predictions <- mclapply(list(Noncoding=noncoding.transcripts,
                Coding=coding.transcripts,
                All=names(coexpressed.Hs)),
                          function(ts) {
                            predict.ontology.performance(ontology.id, transcripts=ts)})


  layout(matrix(c(1,2,3,3), 2, 2, byrow=T))
  for (i in 1:length(predictions)) {
    plot(performance(predictions[[i]], "ppv"),
         lwd=2, main=names(predictions)[i],
         col=rainbow(3)[i])
  }
}

# Potentially useful packages:
## DOSE
##



