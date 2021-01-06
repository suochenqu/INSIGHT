library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if ((length(args) <=0) || (length(args)>=4)){
  stop('usage: Rscript extract_barcode_reads.R input_file [output_file barcode_length]')
}
input_file <- switch(tail(strsplit(args[1], '\\.')[[1]],1),
                     'gz' = gzfile(args[1]),
                     'bz' = bzfile(args[1]),
                     'xz' = xzfile(args[1]),
                     args[1]
                     )  # input file containing sequencing result

output_file <- ifelse(length(args)>=2, args[2], 'pool_reads.csv')  # output file storing barcode reads
amplicon <- 'ACACCTGTGCCTGTTAAACCATTGAAGTTGAAATTGACACATTTGTTTTTAACCAAATTAGTAGACTTTTTAGGTCCACAAACAGTTGCTGG'
barcode_length <- ifelse(length(args)==3, as.integer(args[3]), 5L)  # length of each barcode in a pair
blocksize <- 1000000L  # block of lines by which to process the input file

cat('Extracting barcode reads from', args[1], ':\n')
block_num <- 0L
repeat {
  skip_len <- blocksize * block_num
  lines <- scan(input_file, what=character(), sep="", n=blocksize, skip=skip_len, quiet=TRUE)
  if (length(lines) == 0) break
  cat('Processing lines', skip_len + 1L, 'to', skip_len + length(lines), '...\n')

  # extract valid reads
  reads = agrep(amplicon, lines, costs=list(i=4, d=4, s=1), max.distance=2, value=TRUE)

  # extract barcode sequence from reads
  first_5_start <- 1
  first_5_end <- barcode_length
  last_5_start <- nchar(amplicon) + barcode_length + 1
  last_5_end <- nchar(amplicon) + barcode_length * 2

  df <- data.frame(first=substring(reads, first_5_start, first_5_end),
                   last=substring(reads, last_5_start, last_5_end))

  # save to csv file
  if (block_num==0L){
    write.table(df, file=output_file, append=FALSE, sep=',', row.names=FALSE, col.names=TRUE, quote=FALSE)
  } else {
    write.table(df, file=output_file, append=TRUE, sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
  }

  if (length(lines) != blocksize) break # last block
  block_num <- block_num + 1L
}


