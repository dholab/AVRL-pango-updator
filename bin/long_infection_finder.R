#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# bringing in sequencing run metadata
experiment_number <- args[1]
experiment_date <- as.Date(args[2])
parentdir <- args[3]

# reading in pango lineages and designation dates
lineage_csv <- read.csv(args[4])
dates <- read.csv(args[5])

# ensuring dates are properly formatted
dates$designation_date <- as.Date(dates$designation_date)

# preparing a data frame to hold long infection data, if detected
long_infections <- data.frame(sample = NA,
                              lineage = NA,
                              lineage_designation_date = as.Date(NA),
                              experiment_date = as.Date(NA),
                              experiment_number = NA,
                              pango_version = NA)

for (i in 1:nrow(lineage_csv)){
  
  lineage <- lineage_csv$lineage[i]
  designation_date <- dates[dates$designation_date==lineage, 1]
  
  if ( (experiment_date - designation_date) >= 90 ){
    
    new_row <- c(lineage_csv$taxon[i],
                 lineage,
                 designation_date,
                 experiment_date,
                 experiment_number,
                 lineage_csv$pangolin_version[i])
    long_infections <- rbind(long_infections, new_row)
    
  } else {
    next
  }
  
  
}
long_infections <- long_infections[2:nrow(long_infections),]
rownames(long_infections) <- NULL

if (nrow(long_infections)==0){
  
  null_row <- rep(NA, times = ncol(long_infections))
  null_row[1] <- paste("No putative long infections were identified in experiment",
                       experiment_number, "on", Sys.Date(), sep = " ")
  long_infections <- rbind(long_inf_table, null_row)
  
}

write.csv(long_infections,
          paste(experiment_number,
                "_putative_long_infections_",
                Sys.Date(), 
                ".csv", sep = ""),
          row.names = F, quote = F, na = "")
