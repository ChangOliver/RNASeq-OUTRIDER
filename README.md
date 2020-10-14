# RNASeq-OUTRIDER

## Feature
This tool supports  
+ OUTRIDER analysis and plotting (Quantile-Quantile plot, Expression Rank plot)

## Usage
```
Rscript OUTRIDER.R [input] [output]

Arguments:
        
    [input]:
        If a directory: path to a directory contaning HTSeq data
        If a file:      existing OUTRIDER rds file
        
    [output]:
        Directory to store program output
        The tool should produce:
          1. OUTRIDER data set object (OUTRIDER_dataset.rds)
          2. OUTRIDER_result.csv
          3. imgs folder containing plots and graphs of analysis's result   
```
The script will first merge all htseq-count files in the input directory into one *htseq-counts-all.csv*, then it will perform the analysis.

## Note
Run genTranslationTable.R beforehand if desired genes are not found in the translation table. This process may take a long time as it needs to query [BioTools.fr](https://biotools.fr/human/ensembl_symbol_converter) if the ensembl isn't present in the downloaded HUGO database (hugo.txt).
```
Rscript genTranslationTable.R [path to directory contaning HTSeq data]
```
