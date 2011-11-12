library(SRAdb)

syncSRA <- function(target.dir, sradb.file="SRAmetadb.sqlite") {
  sradb <- dbConnect(SQLite(), sradb.file)
  studies <- dbGetQuery(sradb,
                      "select distinct(study_accession) from study
                       where study_type LIKE \"RNAseq\"")$study_accession
  samples <- dbGetQuery(sradb,
                      "select distinct(sample_accession)
                       from sample where taxon_id=9606")$sample_accession
  runs <- intersect(sraConvert(studies, sra_con=sradb)$run,
                    sraConvert(studies, sra_con=sradb)$run)
  urls <- tempfile()
  writeLines(listSRAfile(runs, sra_con=sradb)$sra, urls)
  system(paste("wget -nc -nd -P", target.dir, "-i", urls))
}

