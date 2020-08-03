matplotlib_palette <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
plotting_symbols <- c(1, 2, 6, 5, 0, 3, 4, 8, 7, 9)
##### 2eL (28 Apr) ######
tmp <- read.csv('data/28_Apr/beacon405060.csv', header=TRUE)
timepoints <- c(0, 3.5, 91.3, 103.8, 116.4, 128.9, 141.6, 154.4, 166.95, 179.4, 193.4, 205.7, 218.05)
timepoints <- timepoints[-c(1,2)] - 91.3
n_cycles <- length(timepoints)
rows <- 1:8 # row A, B, C are reaction wells
reaction_rows <- 1:3
reaction_cols <- 3:9
wells <- as.vector(t(outer(LETTERS[rows], 1:12, paste0)))
n_wells <- length(wells)
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp
dat <- as.matrix(tmp[, -(1:4)])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# visualise control well vs reaction wells
boxplot(dat[ctrl_wells, ], ylim=c(50000, 250000), border='blue')
boxplot(dat[reaction_wells, ], add=TRUE, border='orange')

# normalisation
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[,1]
plot(normalisation, type='b')

# plotting by rows (row A, B, C which are reaction rows)
#palet <- matplotlib_palette[seq_along(reaction_cols)]
palet <- matplotlib_palette[c(7,6,4,3,2,1,5)]
#palet <- rep('black', 7)
pch <- plotting_symbols[c(7,6,4,3,2,1,5)]
pch <- rep(1, 7)
title <- c('Beacon 1:40', 'Beacon 1:50', 'Beacon 1:60')
names(title) <- c('A','B','C')
legend <- c(expression(paste("10"^"6", " copies")),
            expression(paste("10"^"5", " copies")),
            expression(paste("10"^"4", " copies")),
            expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')

r='B'
pdf(paste0('fig/paper figures/2eL.pdf'), 8, 6)
plot(range(timepoints), range(dat_norm[paste0(r, reaction_cols), ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[r])
for (c in reaction_cols){
  points(timepoints, dat_norm[paste0(r, c), ], type='b', col=palet[match(c,reaction_cols)], lwd=2, pch=pch[match(c,reaction_cols)])
}
legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
dev.off()

##### 2eR (30 Apr) ######
tmp <- read.csv('data/30_Apr/beacon_replicates_life_Science_comparasion.csv', header=TRUE)
timepoints <- c(0, 12.8, 25.9, 38.55, 51.05, 65.5, 79.4, 91.8, 105.5, 107.4, 119.1)
n_cycles <- length(timepoints)
rows <- 1:8 # row A, B, C are reaction wells
wells <- as.vector(t(outer(LETTERS[rows], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 2:7; reaction_cols <- c(1:9, 12)
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# discard timepoint 9 and redo boxplots
dat <- dat[, -9]
n_cycles <- n_cycles - 1
timepoints <- timepoints[-9]

# normalisation
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[,1]


# plotting by columns
palet <- matplotlib_palette[c(6,4,3,2,1,5)]
#palet <- rep('black', 6)
pch <- plotting_symbols[c(6,4,3,2,1,5)]
pch <- rep(1, 6)

title <- c('Beacon 1/50', 'Beacon 1/60', 'Beacon 1/75', 'Life Sciences')
legend <- c(expression(paste("10"^"5", " copies")),
            expression(paste("10"^"4", " copies")),
            expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')

cols <- list(c(1,2,3),c(4,5,6),c(7,8,9),12)
fig <- 4
pdf(paste0('fig/paper figures/2eR.pdf'), 8, 6)
fig_wells <- as.vector(outer(LETTERS[reaction_rows], cols[[fig]], paste0))
plot(range(timepoints), range(dat_norm[fig_wells, ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[fig])
for (c in cols[[fig]]) for (r in reaction_rows) {
  points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[match(r,reaction_rows)],pch=pch[match(r,reaction_rows)], lwd=2)
}
legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
dev.off()

##### 3cL and 3cR (12 May) #####
tmp <- read.csv('data/12_May/lifescience_barcode_primers_12052020.csv', header=TRUE)
timepoints <- c(0,12.83, 25.35, 41.17, 53.72, 66.42, 78.85, 91.4, 104.13, 119.72)
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:4; reaction_cols <- 1:6
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]

# plotting by columns
palet <- matplotlib_palette[c(3,2,1,5)]
#palet <- rep('black', 4)
pch <- plotting_symbols[c(3,2,1,5)]
pch <- rep(1, 4)

title <- c('Control (P8)', 'Left illumina', 'R A1', 'R A2', 'R A1 + L B1', 'R A2 + L B2')
legend <- c(expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')
cols <- 1:6
fig <- 5
pdf(paste0('fig/paper figures/3cL.pdf'), 8, 6)
fig_wells <- as.vector(outer(LETTERS[reaction_rows], cols[[fig]], paste0))
plot(range(timepoints), range(dat_norm[fig_wells, ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[fig])
for (c in cols[[fig]]) for (r in reaction_rows) {
  points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[r], pch=pch[r], lwd=2)
}
legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
dev.off()

fig <- 6
pdf(paste0('fig/paper figures/3cR.pdf'), 8, 6)
fig_wells <- as.vector(outer(LETTERS[reaction_rows], cols[[fig]], paste0))
plot(range(timepoints), range(dat_norm[fig_wells, ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[fig])
for (c in cols[[fig]]) for (r in reaction_rows) {
  points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[r], pch=pch[r], lwd=2)
}
legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
dev.off()



##### 4b ######
tmp <- read.csv('data/11_May/TRno262_condition_optimization.csv', header=TRUE)
timepoints <- c(0, 12.28, 25.32, 38.95, 52.02, 67.53, 81.17, 94.37, 111.27, 125.67)
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:6; reaction_cols <- 1:6
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]


# plotting by row
palet <- matplotlib_palette[c(6,4,3,2,1,5)]
#palet <- rep('black', 6)
pch <- plotting_symbols[c(6,4,3,2,1,5)]
pch <- rep(1, 6)

title <- c('In-house mixture', 'Trehalose 8 mM', 'low DMSO (5%)', 'hi-T7', 'T7 high yield', 'T7 flash')
legend <- c(expression(paste("10"^"5", " copies")),
            expression(paste("10"^"4", " copies")),
            expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')
r=1
pdf('fig/paper figures/4b.pdf', 8,6)
plot(range(timepoints), range(dat_norm[paste0(LETTERS[r], reaction_cols), ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[match(r,reaction_rows)])
for (c in reaction_cols){
  points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[match(c,reaction_cols)], pch=pch[match(c,reaction_cols)], lwd=2)
}
legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
dev.off()







##### 4cL, 4cR, 4eL, 4eR (5 May) #####
tmp <- read.csv('data/5_May/lifesci_saliva_isothermal_05052020.csv', header=TRUE)
timepoints <- c(3.75, 18.5, 31.75, 45.15, 58.5, 71.7, 84.75, 103.5)
n_cycles <- length(timepoints)
rows <- 1:8 # row A, B, C, D, E are reaction wells
wells <- as.vector(t(outer(LETTERS[rows], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:5; reaction_cols <- 1:5
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -c(1,2)])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]


# plotting by columns
palet <- matplotlib_palette[c(4,3,2,1,5)]
#palet <- rep('black', 5)
pch <- plotting_symbols[c(4,3,2,1,5)]
pch <- rep(1, 5)

title <- c('Saliva 65°C (repeat 1)', 'Saliva 65°C (repeat 2)', 'control 65°C', 'Saliva 41°C (repeat 1)', 'Saliva 41°C (repeat 2)')
legend <- c(expression(paste("10"^"4", " copies")),
            expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')

filenames <- c('4cL', '4cR', 'control', '4eL', '4eR')
cols <- 1:5
for (fig in c(1,2,4,5)){
  pdf(paste0('fig/paper figures/', filenames[fig],'.pdf'), 8, 6)
  fig_wells <- as.vector(outer(LETTERS[reaction_rows], cols[[fig]], paste0))
  plot(range(timepoints), range(dat_norm[fig_wells, ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[fig])
  for (c in cols[[fig]]) for (r in reaction_rows) {
    points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[r], pch=pch[r], lwd=2)
  }
  legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
  dev.off()
}


##### S1f (25 Apr and 28 Apr) #####
path <- 'data/25_Apr/'
n_cycles <- 10
timepoints <- c(10.5, 23.8, 37.7, 41.15, 52.35, 65.8, 77.3, 90.4, 103.4, 117.1)
reaction_rows <- c(2:5,7)
reaction_cols <- 1:5
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], 1:6, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

dat <- array(0, dim=c(n_wells, n_cycles), dimnames=list(wells, NULL))
for (cycle in 1:n_cycles){
  tmp <- read.csv(paste0(path, cycle+1, '.csv'), header=TRUE, skip=8, row.names = 1)
  dat[,cycle] <- t(as.matrix(tmp[1:8,]))
}

# normalisation
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[,1]


# plotting by columns
palet <- matplotlib_palette[1:5]
#palet <- rep('black', 5)
pch <- plotting_symbols[c(1:5)]
pch <- rep(1, 5)

title <- c('Beacon undiluted', 'Beacon 1:10', 'Beacon 1:20', 'Beacon 1:30', 'Beacon 1:40')
# rows: patient, twist 1e5, spike 1e9, spike 1e8, spike 1e7, spike 1e6, no template, no template/no primer
legend <- c(expression(paste('Twist', "10"^"5", " copies")),
            expression(paste('IVT', "10"^"9", " copies")),
            expression(paste('IVT', "10"^"8", " copies")),
            expression(paste('IVT', "10"^"7", " copies")),
            'no template')
filename <- c('', '1', '2', '3', '4')
for (c in 2:5){
  pdf(paste0('fig/paper figures/S1f', filename[c],'.pdf'), 8, 6)

  plot(range(timepoints), range(dat_norm[paste0(LETTERS[reaction_rows], c), ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[c])
  for (r in reaction_rows) {
    points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[match(r, reaction_rows)], pch=pch[match(r, reaction_rows)], lwd=2)
  }
  legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
  dev.off()
}


tmp <- read.csv('data/28_Apr/beacon405060.csv', header=TRUE)
timepoints <- c(0, 3.5, 91.3, 103.8, 116.4, 128.9, 141.6, 154.4, 166.95, 179.4, 193.4, 205.7, 218.05)
timepoints <- timepoints[-c(1,2)] - 91.3
n_cycles <- length(timepoints)
rows <- 1:8 # row A, B, C are reaction wells
reaction_rows <- 1:3
reaction_cols <- 3:9
wells <- as.vector(t(outer(LETTERS[rows], 1:12, paste0)))
n_wells <- length(wells)
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp
dat <- as.matrix(tmp[, -(1:4)])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[,1]
plot(normalisation, type='b')


# plotting by rows (row A, B, C which are reaction rows)
#palet <- matplotlib_palette[seq_along(reaction_cols)]
palet <- matplotlib_palette[c(7,6,4,3,2,1,5)]
#palet <- rep('black', 7)
pch <- plotting_symbols[c(7,6,4,3,2,1,5)]
pch <- rep(1, 7)

title <- c('Beacon 1:40', 'Beacon 1:50', 'Beacon 1:60')
names(title) <- c('A','B','C')
legend <- c(expression(paste("10"^"6", " copies")),
            expression(paste("10"^"5", " copies")),
            expression(paste("10"^"4", " copies")),
            expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')
filename <- c('4','5','6')
for (fig in 1:3){
  r <- LETTERS[fig]
  pdf(paste0('fig/paper figures/S1f', filename[fig], '.pdf'), 8, 6)
  plot(range(timepoints), range(dat_norm[paste0(r, reaction_cols), ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[r])
  for (c in reaction_cols){
    points(timepoints, dat_norm[paste0(r, c), ], type='b', col=palet[match(c,reaction_cols)], pch=pch[match(c,reaction_cols)], lwd=2)
  }
  legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
  dev.off()
}


######## new1abc (15 Jun) ########

tmp <- read.csv('data/15_Jun/959865_15062020.csv', header=TRUE)
timepoints <- c(0, 15.33, 28.22, 40.85, 53.38, 69.42, 77.28, 90.82, 107.32, 179.02)
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- c(1,2,4); reaction_cols <- 1:4
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]

# plotting by rows
palet <- matplotlib_palette[c(3,2,1,5)]
#palet <- rep('black', 7)
pch <- plotting_symbols[c(3,2,1,5)]
pch <- rep(1, 7)

title <- c('5 min at 95 °C', '2 min at 65 °C', ' ', '2 min at 98 °C')
legend <- c(expression(paste("10"^"3", " copies")),
            expression(paste("10"^"2", " copies")),
            expression(paste("10"^"1", " copies")),
            'no template')

filename <- c('a', 'b', ' ', 'c')
for (r in reaction_rows){
  pdf(paste0('fig/paper figures/new1', filename[r], '.pdf'), 8, 6)
  fig_wells <- paste0(LETTERS[r], reaction_cols)
  plot(range(timepoints), range(dat_norm[fig_wells, ]), pch=' ', xlab='minutes', ylab='normlalised fluorescence', main=title[r])

  for (c in reaction_cols) {
    points(timepoints, dat_norm[paste0(LETTERS[r], c), ], type='b', col=palet[c], lwd=2)
  }
  legend('topleft', col=palet, pch=pch, lty=1, legend=legend, cex=0.8, lwd=2)
  dev.off()
}

##### new2 (20 Jul) #####
date <- '20_Jul'; infile <- 'data/20_Jul/saliva_input_20072020.csv'
reaction_rows <- 1:4; reaction_cols <- 1:6
title <- paste0('Column ', reaction_cols)
legend <- paste0('Row ', reaction_rows)

pch=20
palet <- matplotlib_palette[c(3,2,1,5)]
lines <- readLines(infile)
lines <- lines[lines!='']; lines <- lines[-(1:6)]
writeLines(lines, paste0(infile, '.processed'))
tmp <- read.csv(paste0(infile, '.processed'), header=TRUE); tmp[,1] <- paste0(tmp[,1],tmp[,2]); tmp <- tmp[,-c(2,3)]; colnames(tmp)[1]<-'Well'

library(stringr)
timepoints <- as.numeric(matrix(as.numeric(unlist(strsplit(str_replace(str_replace(paste0(strsplit(lines[1],',')[[1]][-c(1:3)], '0'), ' s0', ''), ' min ', ','), ','))),ncol=2, byrow=T)%*%c(1,1/60))[-4] # compute time
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -c(1,5)])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]
plot(normalisation, type='b')

# all wells normalised plot
ylim <- range(dat_norm)
pdf(paste0('fig/paper figures/new2.pdf'), 10, 5)
par(mfrow=c(length(reaction_rows),length(reaction_cols)), mar=c(2,2,2,2))
for (row in reaction_rows){
  for (col in reaction_cols){
    well <- paste0(LETTERS[row], col)
    plot(timepoints, dat_norm[well,], xlab='minutes', ylab='normalised fluorescence', col=palet[row], pch=pch, type='b', lwd=2, ylim=ylim)
  }
}
dev.off()

##### new3 (16 Jul) #####
date <- '16_Jul'; infile <- 'data/16_Jul/160720_seq5.csv'
palet <- putils::matplotlib_palette(10)
lines <- readLines(infile)
lines <- lines[lines!='']; lines <- lines[-(1:6)]
writeLines(lines, paste0(infile, '.processed'))
tmp <- read.csv(paste0(infile, '.processed'), header=TRUE); tmp[,1] <- paste0(tmp[,1],tmp[,2]); tmp <- tmp[,-c(2,3)]; colnames(tmp)[1]<-'Well'

library(stringr)
timepoints <- as.numeric(matrix(as.numeric(unlist(strsplit(str_replace(str_replace(paste0(strsplit(lines[1],',')[[1]][-c(1:3)], '0'), ' s0', ''), ' min ', ','), ','))),ncol=2, byrow=T)%*%c(1,1/60)) # compute time
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:8; reaction_cols <- c(1,3,5,7,9,11)
reaction_wells <- c(as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0))), 'H2','H4','H6','H8','H10','H12')
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]

# all wells normalised plot
palet <- matplotlib_palette[c(7,6,4,3,2,1,5,5)]
ylim <- range(dat_norm)
pdf(paste0('fig/paper figures/new3.pdf'), 10, 10)
par(mfrow=c(8,6), mar=c(2,2,2,2))
for (row in 1:8){
  for (col in reaction_cols){
    well <- paste0(LETTERS[row], col)
    plot(timepoints, dat_norm[well,], xlab='minutes', ylab='normalised fluorescence', col=palet[row], pch=20, type='b', lwd=2, ylim=ylim)
  }
}
dev.off()

# compute time to detection
ttd <- matrix(0, length(reaction_rows), length(reaction_cols), dimnames=list(LETTERS[reaction_rows], reaction_cols))

for (i in seq_along(reaction_rows)){
  for (j in seq_along(reaction_cols)){
    well <- paste0(LETTERS[reaction_rows[i]], reaction_cols[j])
    ttd[i, j] <- time_to_detection(timepoints, dat_norm[well,], 0.2)
  }
}

print(round(ttd, 1))

##### new4 (14 Jul and 15 Jul) #####
# data from 14 Jul
tmp <- read.csv('data/14_Jul/9541_14072020.csv', header=TRUE)
timepoints <- c(0, 12.82, 29.65, 41.32, 53.08, 71.33, 90.6)
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:4; reaction_cols <- 1:3
reaction_wells <- as.vector(t(outer(LETTERS[reaction_rows], reaction_cols, paste0)))
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]

plotting_wells_14 <- as.vector(outer(LETTERS[1:4], 1:3, paste0))[1:11]
dat_14 <- dat_norm
timepoints_14 <- timepoints

# data from 15 Jul
date <- '15_Jul'; infile <- 'data/15_Jul/150720_small_copy_beacon.csv'
lines <- readLines(infile)
lines <- lines[lines!='']; lines <- lines[-(1:6)]
writeLines(lines, paste0(infile, '.processed'))
tmp <- read.csv(paste0(infile, '.processed'), header=TRUE); tmp[,1] <- paste0(tmp[,1],tmp[,2]); tmp <- tmp[,-c(2,3)]; colnames(tmp)[1]<-'Well'

library(stringr)
timepoints <- as.numeric(matrix(as.numeric(unlist(strsplit(str_replace(str_replace(paste0(strsplit(lines[1],',')[[1]][-c(1:3)], '0'), ' s0', ''), ' min ', ','), ','))),ncol=2, byrow=T)%*%c(1,1/60)) # compute time
n_cycles <- length(timepoints)
wells <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))
n_wells <- length(wells)
reaction_rows <- 1:4; reaction_cols <- 1:3
ctrl_wells <- wells[!(wells%in%reaction_wells)]

# read off data matrix from tmp (remove timepoint 0)
dat <- as.matrix(tmp[, -1])
colnames(dat) <- paste0('t', seq_len(n_cycles))
rownames(dat) <- tmp$Well

# normalisation (to median of blank well)
normalisation <- rep(0, n_cycles)
for (cycle in 1:n_cycles) {
  normalisation[cycle] <- median(dat[ctrl_wells, cycle])
}
dat_norm <- sweep(dat, 2, normalisation, '/')
dat_norm <- dat_norm - dat_norm[, 1]
plot(normalisation, type='b')

plotting_wells_15 <- c('D4', as.vector(outer(LETTERS[1:4], 1:3, paste0)))
dat_15 <- dat_norm
timepoints_15 <- timepoints

# all wells normalised plot
palet <- matplotlib_palette[c(2,1,3,5)]
ylim <- range(dat_norm)
pdf(paste0('fig/paper figures/new4_fluorescence_sensitivity.pdf'), 10, 5)
par(mfcol=c(4,6), mar=c(2,2,2,2))
counter <- 0
ylim <- range(dat_14[plotting_wells_14,], dat_15[plotting_wells_15,])
for (well in plotting_wells_14){
  counter <- counter + 1
  plot(timepoints_14, dat_14[well,], xlab='minutes', ylab='normalised fluorescence', col=palet[(counter-1)%%4+1], pch=20, type='b', lwd=2, ylim=ylim)
}
for (well in plotting_wells_15){
  counter <- counter + 1
  plot(timepoints_15, dat_15[well,], xlab='minutes', ylab='normalised fluorescence', col=palet[(counter-1)%%4+1], pch=20, type='b', lwd=2, ylim=ylim)
}

dev.off()
