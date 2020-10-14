suppressMessages(library(httr))
suppressMessages(library(jsonlite))

combine.htseq <- function(cntDir, pat, outFile){
  # adapted from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4#Combine_individual_HTSeq_files_from_the_.27all.27_mapping_series
  
  tophat.all <- list.files(path = cntDir,
                           pattern = pat,
                           all.files = TRUE,
                           recursive = FALSE,
                           ignore.case = FALSE,
                           include.dirs = FALSE)
  
  # we choose the 'all' series
  myfiles <- tophat.all
  DT <- list()
  
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  for (i in 1:length(myfiles) ) {
    infile = paste(cntDir, myfiles[i], sep = "/")
    DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
    cnts <- gsub("(.*).htseq-count.txt", "\\1", myfiles[i])
    colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
  }
  
  # merge all elements based on first ID columns
  data <- DT[[myfiles[1]]]
  
  # we now add each other table with the ID column as key
  for (i in 2:length(myfiles)) {
    y <- DT[[myfiles[i]]]
    z <- merge(data, y, by = c("ID"))
    data <- z
  }
  
  # ID column becomes rownames
  rownames(data) <- data$ID
  data <- data[,-1]
  
  ## add total counts per sample
  data <- rbind(data, tot.counts=colSums(data))
  
  # take all data rows to a new table
  data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]
  
  # write data to file
  write.csv(data.all, file = paste0(cntDir, outFile))
  
  # cleanup intermediate objects
  rm(y, z, i, DT)
}

get.translation <- function(res, col){
  idSymbolTable <- read.csv(file = "../data/hugo.csv", sep = '\t')
  hugoDict <- as.vector(idSymbolTable$Approved.symbol)
  names(hugoDict) <- idSymbolTable$Ensembl.gene.ID
  message(paste("Total gene:", nrow(res)))
  
  for (i in c(1:nrow(res))){
    if ((i+1) %% 5000 == 0){
      message(paste(i+1, "gene translated"))
    }
    
    id <- sub("(ENS[A-Z]+[0-9]{11}).*", "\\1", res[i, 1])

    if (id %in% names(hugoDict)){
      res[i, col] <- hugoDict[[id]]
    }else{
      url = paste("https://biotools.fr/human/ensembl_symbol_converter/?api=1&id=", id, sep="")
      r <- GET(url)
      output = fromJSON(content(r, "text", encoding = "UTF-8"), flatten=TRUE)

      res[i, col] <- ifelse((length(output) != 0 & !is.null(output[[id]])), output[[id]], res[i, 1])
    }
  }
  
  return(res)
}

print("Translating Ensembl to Symbol...(This may take a while)")

args <- commandArgs(trailingOnly = TRUE)

cntDir <- args[[1]]

combine.htseq(cntDir, "*.htseq-count.txt", "/htseq-counts-all.csv")
counts <- paste0(cntDir, "/htseq-counts-all.csv")
ctsTable <- read.table(counts, check.names=FALSE, header=TRUE, sep = ",", stringsAsFactors = FALSE)

df <- data.frame(as.matrix(ctsTable[,1]), stringsAsFactors = FALSE)
colnames(df)[1] <- "Ensembl.ID"
df[ , "Gene.name"] <- NA
df <- get.translation(df, 2)

message("Dumping translation table")
write.csv(df, "../data/translation_table.csv", row.names = FALSE)