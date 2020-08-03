##### library and auxiliary functions #####
library(stringr)
library(stringdist)
library(MASS)
pairing <- setNames(c('T', 'G', 'C', 'A'), c('A','C','G','T'))
burst <- function(seq) strsplit(seq, '')[[1]]
compl <- function(seq) paste0(pairing[rev(burst(seq))], collapse='')
hamming <- function(u, v){stringdist(u, v, method='hamming')}

##### load used barcode combinations and corresponding Stage 1 results #####
# barcode_used.csv is a comma separated file such that each row contains
# information of a specific barcode pair used in NASBA reaction, including the
# following fields:
#     right barcode sequence,
#     right barcode name,
#     left barcode sequence,
#     left barcode name,
#     input copies of viral RNA in NASBA reaction,
#     time to detection using fluorescence readout in Stage 1.
barcode_used <- read.csv('data/barcode_hopping/barcode_used.csv',
                         stringsAsFactors=F)
barcode_used$sequence <- paste0(barcode_used$right_sequence,
                                barcode_used$left_sequence)
barcode_used$name <- paste0(barcode_used$right_name, barcode_used$left_name)
barcode_used$pool <- ifelse(is.na(barcode_used$time_to_detection), 'neg', 'pos')

##### load sequencing results #####
all_seq <- readLines('data/barcode_hopping/all_linear_trimmed.seq.gz') # 102nt seq
pool_sequence <- data.frame(
  'first_5' = substr(all_seq, 1, 5),
  'last_5' = substr(all_seq, 6, 10)
)
pool_sequence$combined <- paste0(pool_sequence$first_5, pool_sequence$last_5)
##### match sequences to barcodes #####
# generate hamming distance matrix with rows representing sequenced barcodes
# and columns representing used barcodes
dist <- outer(pool_sequence$combined, barcode_used$sequence, hamming)
barcode_used$exact_matches <- colSums(dist == 0)
barplot(sort(barcode_used$exact_matches, decreasing=T),
        las=2, cex.names=0.8, log='y', offset=1,
        names.arg=barcode_used$name, ylab='Read Count', axes=F)
axis(2,at=10^(0:5))

## match to left and right barcodes separately ##
# first_5 matched to right barcodes and last_5 matched to left barcodes
right_dist <- outer(pool_sequence$first_5, barcode_used$right_sequence, hamming)
left_dist <- outer(pool_sequence$last_5, barcode_used$left_sequence, hamming)

pool_sequence$right_hamming_err <- apply(right_dist, 1, min)
pool_sequence$right_called_name <- barcode_used$right_name[
  apply(right_dist, 1, which.min)]

pool_sequence$left_hamming_err <- apply(left_dist, 1, min)
pool_sequence$left_called_name <- barcode_used$left_name[
  apply(left_dist, 1, which.min)]

pool_sequence$hamming_err <- with(pool_sequence,
                                  right_hamming_err + left_hamming_err)
pool_sequence$called_name <- paste0(pool_sequence$right_called_name,
                                    pool_sequence$left_called_name)

# create a table for all called barcodes, including used and hopped
calling_table <- with(pool_sequence,
                      sort(table(called_name[hamming_err == 0]), decreasing=T))
calling_table <- as.data.frame(calling_table, stringsAsFactors=F)
colnames(calling_table) <- c('name', 'read_count')

calling_table$used <- calling_table$name %in% barcode_used$name
calling_table$right_name <- with(pool_sequence,
                                 right_called_name[match(calling_table$name, called_name)])
calling_table$left_name <- with(pool_sequence,
                                left_called_name[match(calling_table$name, called_name)])

# total used barcode readcounts involving called right/left barcodes
calling_table$left_count <- calling_table$right_count <- rep(0, nrow(calling_table))
for (i in 1:nrow(calling_table)){
  if (!calling_table$used[i]) {
    calling_table$right_count[i] <- with(calling_table,
                                         sum(read_count[(right_name == right_name[i]) & used]))
    calling_table$left_count[i] <- with(calling_table,
                                        sum(read_count[(left_name == left_name[i]) & used]))
  }
}

# copy right/left sequences and pool name from barcode_used to calling_table
calling_table$right_sequence <- barcode_used$right_sequence[
  match(calling_table$right_name, barcode_used$right_name)]
calling_table$left_sequence <- barcode_used$left_sequence[
  match(calling_table$left_name, barcode_used$left_name)]
calling_table$sequence <- paste0(calling_table$right_sequence,
                                 calling_table$left_sequence)

calling_table$right_pool <- barcode_used$pool[match(calling_table$right_name,
                                                    barcode_used$right_name)]
calling_table$left_pool <- barcode_used$pool[match(calling_table$left_name,
                                                   barcode_used$left_name)]

# count all reads from used barcodes within hamming distance 1 to called barcode
dist_to_used <- outer(calling_table$sequence, barcode_used$sequence, hamming)
calling_table$hamming1read <- apply(dist_to_used == 1, 1,
                                    function(v) sum(barcode_used$exact_matches[v]))



##### fit negative binomial model #####
# fit hopping model using hopped barcode pairs
model <- glm.nb(read_count ~ left_count*right_count + hamming1read,
                link='identity', data=calling_table[!calling_table$used,])
AIC(model)
summary(model)
beta <- model$coefficients  # regression coefficients
theta <- summary(model)$theta  # overdispersion parameter

calling_table$pred <- predict(model, calling_table) # fitted hopped count

# Z score calculation
calling_table$Z_score <- (calling_table$read_count - calling_table$pred) /
  (qnbinom(0.975, size=theta, mu=calling_table$pred) - calling_table$pred) * 2

plot(log(calling_table$Z_score+2), col=ifelse(calling_table$used, 'red', 'black'),
     pch=20, cex=0.5)


##### plots #####
ymax <- max(calling_table$read_count)
plot(y=calling_table$read_count[!calling_table$used], x=fitted(model),
     xlim=c(1,ymax), ylim=c(1,ymax), pch=20, cex=0.5,
     xlab='predicted reads', ylab='actual reads', log='xy');
abline(a=0, b=1, col='blue')
xs <- exp(seq(0,log(ymax),length=1000)); sd <- sqrt(xs + xs^2/theta)
alpha=0.05
y_upper <- qnbinom(1-alpha/2, size=theta, mu=xs)
y_lower <- qnbinom(alpha/2, size=theta, mu=xs)
points(xs, y_upper, type='l', lty=2, col='blue')
points(xs, y_lower, type='l', lty=2, col='blue')

alpha=0.0001
y_upper <- qnbinom(1-alpha/2, size=theta, mu=xs)
y_lower <- qnbinom(alpha/2, size=theta, mu=xs)
points(xs, y_upper, type='l', lty=2, col='red')
points(xs, y_lower, type='l', lty=2, col='red')

points(calling_table$pred[calling_table$used],
       calling_table$read_count[calling_table$used], pch=20, cex = 0.5, col = 'green')

