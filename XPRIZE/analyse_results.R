##### library and auxiliary functions #####
library(stringr)
library(stringdist)
library(MASS)
pairing <- setNames(c('T', 'G', 'C', 'A'), c('A','C','G','T'))
burst <- function(seq) strsplit(seq, '')[[1]]
compl <- function(seq) paste0(pairing[rev(burst(seq))], collapse='')
hamming <- function(u, v){stringdist(u, v, method='hamming')}



##### preprocess #####
# INPUT: two files barcode_info_file and pool_reads_file
#   1. the barcode_info_file should be a csv file where each row contains
#      information of a specific barcode pair used in stage 1 NASBA reactions,
#      incuding the following fields
#      - ID: sample ID for which the barcode pair is used
#      - right_sequence: 5nt right barcode sequence
#      - left_sequence: 5nt left barcode sequence
#      - stage_1_result: positive or negative
#   2. the pool_reads_file should be a csv file containing sequencing results
#      in the pool, including the following fields
#      - first: first 5 barcode nucleotides in the sequence, this should
#        correspond to the right sequence of some used barcode pair
#      - last: last 5 barcode nucleotides in the sequence, this should
#        correspond to the reverse complement of the left sequence of used
#        barcode pairs
#
# OUTPUT: two dataframes barcode_used and calling_table
#   1. barcode_used contains all info from barcode_info_file, as well as the
#      number of exact matches for each used barcode pair in the pool
#   2. calling_table contains all 'called sequences' of barcode sequences in
#      the pool. Here, we 'call' the first_5 and the reverse complement of
#      last_5 of pool_read_file to one of the right barcodes and left barcodes
#      respectively. For each called sequence, we record its read_count, weather
#      it is directly from a used barcode pair (as opposed from barcode hopping),
#      total counts of other used barcode pairs sharing the same right barcode
#      (right_count), total counts of other used barcode pairs sharing the same
#      left barcode (left_count), and total counts of used barcode pairs at
#      Hamming distance 1 to this called sequence.
#
# Plot:
# a figure reads_barplot.pdf is generated summarising read counts of all used
# barcode pairs.

preprocess <- function(barcode_info_file='barcode_used.csv',
                       pool_reads_file='all_pool_reads.csv'){
  ## load used barcode information ##
  barcode_used <- read.csv(barcode_info_file, stringsAsFactors=F)
  barcode_used$left_sequence <- sapply(toupper(barcode_used$left_sequence), compl)
  barcode_used$right_sequence <- toupper(barcode_used$right_sequence)
  barcode_used$sequence <- paste0(barcode_used$right_sequence,
                                  barcode_used$left_sequence)

  ## load sequencing results ##
  pool_sequence <- read.csv(pool_reads_file, stringsAsFactors=F)
  pool_sequence$combined <- paste0(pool_sequence$first, pool_sequence$last)
  barcode_length <- nchar(pool_sequence$first[1])  # barcode length for each barcode in pair, default = 5

  ## match sequences to barcodes ##
  # generate hamming distance matrix with rows representing sequenced barcodes
  # and columns representing used barcodes
  dist <- outer(pool_sequence$combined, barcode_used$sequence, hamming)
  barcode_used$exact_matches <- colSums(dist == 0)

  ## match to left and right barcodes separately ##
  # first matched to right barcodes and last matched to left barcodes
  right_dist <- outer(pool_sequence$first, barcode_used$right_sequence, hamming)
  left_dist <- outer(pool_sequence$last, barcode_used$left_sequence, hamming)

  pool_sequence$right_hamming_err <- apply(right_dist, 1, min)
  pool_sequence$right_called <- barcode_used$right_sequence[
    apply(right_dist, 1, which.min)]

  pool_sequence$left_hamming_err <- apply(left_dist, 1, min)
  pool_sequence$left_called <- barcode_used$left_sequence[
    apply(left_dist, 1, which.min)]

  pool_sequence$hamming_err <- with(pool_sequence,
                                    right_hamming_err + left_hamming_err)
  pool_sequence$called_seq <- paste0(pool_sequence$right_called,
                                     pool_sequence$left_called)

  ## create a table for all called barcodes, including used and hopped ##
  # get relevant info from pool_sequence
  calling_table <- with(pool_sequence,
                        sort(table(called_seq[hamming_err == 0]),
                             decreasing=T))
  calling_table <- as.data.frame(calling_table, stringsAsFactors=F)
  colnames(calling_table) <- c('sequence', 'read_count')

  calling_table$used <- calling_table$sequence %in% barcode_used$sequence
  calling_table$right_seq <- substring(calling_table$sequence, 1, barcode_length)
  calling_table$left_seq <- substring(calling_table$sequence, barcode_length+1, barcode_length*2)

  # total used barcode readcounts involving called right/left barcodes
  calling_table$left_count <- calling_table$right_count <- rep(0, nrow(calling_table))
  for (i in 1:nrow(calling_table)){
    calling_table$right_count[i] <- with(calling_table,
                        sum(read_count[(right_seq == right_seq[i]) & used]))
    calling_table$left_count[i] <- with(calling_table,
                        sum(read_count[(left_seq == left_seq[i]) & used]))
    if (calling_table$used[i]) {
      calling_table$left_count[i] <- calling_table$left_count[i] - calling_table$read_count[i]
      calling_table$right_count[i] <- calling_table$right_count[i] - calling_table$read_count[i]
    }
  }

  # count all reads from used barcodes within hamming distance 1 to called barcode
  dist_to_used <- outer(calling_table$sequence, barcode_used$sequence, hamming)
  calling_table$hamming1read <- apply(dist_to_used == 1, 1,
                                      function(v) sum(barcode_used$exact_matches[v]))


  ## generate read count barplot for all used barcode pairs ##
  ord <- order(barcode_used$exact_matches, decreasing=T)
  plotfile <- paste0(strsplit(pool_reads_file,'\\.')[[1]][1], '_count_barplot.pdf')
  pdf(plotfile, width = 14, height=6)
  barplot(barcode_used$exact_matches[ord],
          las=2, cex.names=0.5, log='y', offset=1,
          col=ifelse(barcode_used$stage_1_result[ord]=='positive',2,1),
          names.arg=barcode_used$ID[ord], ylab='Read Count', axes=F)
  axis(2,at=10^seq(0, ceiling(log10(barcode_used$exact_matches[ord[1]]))))
  legend('topright', col=c(1,2), legend=c('Stage 1 negative', 'Stage 1 positive'), pch=15)
  dev.off()

  ## return output ##
  return(list(barcode_used=barcode_used, calling_table=calling_table))
}


##### barcode hopping modelling #####
# INPUT: barcode_used and calling_table dataframe from the preprocessing function
#
# OUTPUT: modified barcode_used dataframe, with computed Z_score and stage 2
#         test results. Here, the Z scores represents devation from predicted
#         reads due to barcode hopping alone. High Z scores mean that the
#         barcode pair is more likely to exist in the pool.
#
# Note: if number of hopped sequences is at most 10, or if the max count of
#       hopped sequences is at most 10, then, the max count of hopped sequences
#       is used as cutoff for stage 2 positivity
#
# Plot:
#    1. Z score barplot for each barcode pair used
#    2. Predicted vs actual read counts for all called barcodes. If barcode
#       hopping are well-accounted, all grey dots should lie within confidence
#       bands
barcode_hopping_modelling <- function(barcode_used, calling_table, Zscore_thresh=10){
  ## check number of hopped sequences
  is_hopped <- !calling_table$used
  max_hop <- max(calling_table$read_count[is_hopped])
  if (sum(is_hopped) <= 10 || max_hop <= 10) {
    barcode_used$stage_2_result <- ifelse(barcode_used$exact_matches > max_hop,
                                          'positive', 'negative')
    return(barcode_used)
  }

  ## fit negative binomial model ##
  # fit hopping model using hopped barcode pairs
  mtry <- try(model <- glm.nb(read_count ~ left_count*right_count + hamming1read,
                    link='identity', data=calling_table[!calling_table$used, ]),
              silent=TRUE)
  if (inherits(mtry, "try-error")){
    print('Negative binomial model fitting failed...')
    print('using max hopped count as cutoff instead.')
    barcode_used$stage_2_result <- ifelse(barcode_used$exact_matches > max_hop,
                                          'positive', 'negative')
    return(barcode_used)
  }

  beta <- model$coefficients  # regression coefficients
  theta <- summary(model)$theta  # overdispersion parameter
  calling_table$pred <- predict(model, calling_table) # fitted hopped count
  calling_table$Z_score <- with(calling_table,
                                (read_count - pred) / (qnbinom(0.975, size=theta, mu=pred) - pred) * 2)
  barcode_used$pred <- calling_table$pred[match(barcode_used$sequence, calling_table$sequence)]

  ## Z score calculation ##
  barcode_used$Z_score <- with(barcode_used,
                               (exact_matches - pred) / (qnbinom(0.975, size=theta, mu=pred) - pred) * 2)
  barcode_used$Z_score[is.na(barcode_used$Z_score)] <- 0
  barcode_used$stage_2_result <- ifelse(barcode_used$Z_score >= Zscore_thresh,
                                        'positive',
                                        'negative')

  ## generate plots ##
  # Z score barplot for each barcode pair used #
  ord <- order(barcode_used$Z_score, decreasing=T)
  pdf('Zscore_barplot.pdf', width = 14, height=6)
  barplot(barcode_used$Z_score[ord],
          las=2, cex.names=0.5, log='y', offset=2,
          col=ifelse(barcode_used$stage_1_result[ord]=='positive',2,1),
          names.arg=barcode_used$ID[ord], ylab='Z score', axes=F)
  axis(2,at=10^seq(0, ceiling(log10(barcode_used$Z_score[ord[1]]))))
  abline(h=Zscore_thresh, lty=2, lwd=2)
  legend('topright', col=c(1,2, 1), legend=c('Stage 1 negative', 'Stage 1 positive', 'Zscore threshold'), pch=c(15,15,NA), lty=c(NA,NA,2))
  dev.off()

  # predicted vs actual plot #
  calling_table$colour <- 'grey'
  calling_table$colour[calling_table$used] <- 'black'
  calling_table$colour[match(barcode_used$sequence[barcode_used$stage_1_result=='positive'], calling_table$sequence)] <- 'red'

  pdf('predict_vs_actual_reads.pdf')
  ymax <- max(calling_table$read_count)
  plot(y=calling_table$read_count, x=calling_table$pred,
       xlim=c(1,ymax), ylim=c(1,ymax), pch=20, cex=0.5,
       xlab='predicted reads', ylab='actual reads', log='xy',
       col=calling_table$colour);
  abline(a=0, b=1, col='blue')
  xs <- exp(seq(0,log(ymax),length=1000)); sd <- sqrt(xs + xs^2/theta)
  alpha <- 0.05
  y_upper <- qnbinom(1-alpha/2, size=theta, mu=xs)
  y_lower <- qnbinom(alpha/2, size=theta, mu=xs)
  points(xs, y_upper, type='l', lty=2, col='blue')
  points(xs, y_lower, type='l', lty=2, col='blue')

  alpha <- 0.0001
  y_upper <- qnbinom(1-alpha/2, size=theta, mu=xs)
  y_lower <- qnbinom(alpha/2, size=theta, mu=xs)
  points(xs, y_upper, type='l', lty=2, col='red')
  points(xs, y_lower, type='l', lty=2, col='red')

  points(calling_table$pred[calling_table$used],
         calling_table$read_count[calling_table$used], pch=20, cex=0.5, col=calling_table$colour[calling_table$used])
  legend('bottomright', pch=20, col=c('red','black','grey'), legend=c('Stage 1 positive','Stage 1 negative','hopped barcodes'))
  dev.off()

  # return output
  return(barcode_used)
}




##### main #####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop('usage: Rscript analyse_results.R barcode_info_file Zscore_thresh all_pool_reads_file [neg_pool_reads_file]')
}
Zscore_thresh <- as.numeric(args[2])
output_file <- 'results.csv'

# analyse all pool reads
cat('Analysing all pool reads ...\n')
tmp <- preprocess(barcode_info_file=args[1],
                  pool_reads_file=args[3])
ret <- barcode_hopping_modelling(tmp$barcode_used, tmp$calling_table, Zscore_thresh)
stage_2_pos <- ret$stage_2_result == 'positive'

# analyse negative pool reads
if (length(args) == 4){
  cat('Analysing negative pool reads ...\n')
  tmp2 <- preprocess(barcode_info_file=args[1],
                     pool_reads_file=args[4])
  ret2 <- barcode_hopping_modelling(tmp2$barcode_used, tmp2$calling_table, Zscore_thresh)

  # combine pools for stage 2 results
  stage_2_pos <- stage_2_pos | (ret2$stage_2_result=='positive')
}

# combine stage 1 and stage 2 results
stage_1_pos <- ret$stage_1_result=='positive'
result <- rep('negative', length(stage_1_pos))
result[stage_1_pos | stage_2_pos] <- 'suspected'
result[stage_1_pos & stage_2_pos] <- 'positive'

# write output to csv file
df <- tmp$barcode_used[, 1:4]
df$stage_2_result <- ifelse(stage_2_pos, 'positive', 'negative')
df$combined_result <- result
write.csv(df, file=output_file, row.names=F)
cat('Done. Output saved in', output_file, '\n')

