library(dplyr)
library(ggplot2)
library(data.table)

# in the transcripts_count directory
# prepration : convert fasta file to table 
#                 awk -f 'FaToTb' *RSEM.transcripts.fa > transcripts_seq.tb

prefix <- "proj_db"
db_name <- paste0(prefix, ".SQLITE") 

database <- src_sqlite(db_name, create = TRUE)



map_filename <- "Trinity.fasta.gene_trans_map"  #paste0(prefix, "_gene_trans_map")
map_tb <- fread(map_filename, header = FALSE)
colnames(map_tb) <- c("gene_name", "transcript_name")
str(map_tb)


TPM_filename <- "rsem-gene.isoform.TPM.not_cross_norm"  #paste0(prefix, "_TPM.matrix")
TPM_tb <- fread(TPM_filename, header = TRUE)
colnames(TPM_tb)[1] <- "transcript_name"
str(TPM_tb)


# seq_filename <- paste0(db_name, "_trans.fa")
# command <- paste0("awk -f 'FaToTb' ./DS", db_num,
#                   "_trans.fa > ./DS", db_num, "_trans.tb")
# system(command, show.output.on.console = TRUE)
# shell command: convert fasta file to table
#    awk -f 'FaToTb' DS3*fa > DS3_trans.tb


seq_filename <- "Trinity.fasta.RSEM.transcripts.fa"  #paste0(prefix, "_trans.tb")
seq_tb <- fread(seq_filename, header = FALSE)
colnames(seq_tb) <- c("transcript_name", "seq")
str(seq_tb)



blast_filenames <- list.files(pattern="*.outfmt6")
blast_tb <- rbindlist(lapply(blast_filenames, fread, data.table=FALSE, header = FALSE))
colnames(blast_tb) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

    # qseqid: query (e.g., gene) sequence id
    # sseqid: subject (e.g., reference genome) sequence id
    # pident: percentage of identical matches
    # length: alignment length
    # mismatch: number of mismatches
    # gapopen: number of gap openings
    # qstart: start of alignment in query
    # qend: end of alignment in query
    # sstart: start of alignment in subject
    # send: end of alignment in subject
    # evalue: expect value
    # bitscore: bit score

str(blast_tb)



copy_to(database, map_tb, temporary = FALSE, overwrite = TRUE)
copy_to(database, TPM_tb, temporary = FALSE, overwrite = TRUE)
copy_to(database, seq_tb, temporary = FALSE, overwrite = TRUE)
copy_to(database, blast_tb, temporary = FALSE, overwrite = TRUE)
