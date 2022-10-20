#!/usr/bin/env Rscript

all_seqs <- read.delim("pango_consensus_sequences.fasta", header = F)
deflines <- data.frame(line = all_seqs[grepl(">", all_seqs$V1),],
                       line_number = which(grepl(">", all_seqs$V1)))
ba_2_row <- which(deflines$line==">BA.2")
next_seq <- deflines[ba_2_row+1, "line_number"]
ba_2 <- data.frame(V1 = all_seqs[deflines[ba_2_row,"line_number"]:(next_seq-1),])

write.table(ba_2, "ba_2_ref.fasta", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
