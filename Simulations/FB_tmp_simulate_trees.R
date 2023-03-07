rm(list = ls())

library(ape)
library(plyr)
library(adephylo)

# not so easy to infer rates in simple way
# symmetrical/asym model for rates?
# state of the root well specified/inferred?
# rate not too fast, not too slow
# dataset size?

# WHY DO WE NEED TO CORRECT THE RATE ESTIMATES BY FACTOR 2?
# this is because in the function definition of rTraitDisc
# there is a step where the rate matrix is multiplied by freq, Q <- Q * freq
# and freq is by default rep(1/k, k). Not sure why given the definition of freq
# thus the simulations are not with a factor 1/k to the rates
# and this should be corrected by setting freq = c(1,1)

#################################### FIRST SET OF SIMULATIONS WITH TREE OF 1000 TIPS  #################################### 

nrep <- 50
tr <- rtree(n = 1000, rooted = T)
median_length <- median(diag(vcv.phylo(tr))) # root to tip distances
cols <- c("blue", "red"); names(cols) <- c("A", "B")
all_rates <- c()
for(rate_multiplier in c(0.5, 0.75, 1, 1.5, 2)){
  
  print(rate_multiplier)
  myrates <- rate_multiplier * c(1/median_length, 2/median_length)
  
  cat('My rates', myrates,'\nSum: ',sum(myrates), '\n')
  for(i in 1:nrep){
    #print(i)
    z <- rTraitDisc(phy = tr, model = "ARD", k = 2, rate = myrates, root.value = sample(x = 1:2, size = 1)) # phenotype
    #print(table(z))
    # plot.phylo(x = tr, tip.color = cols[z])
    rates <- ace(x = z, phy = tr, type = "discrete", method = "ML", model = "ARD", ip = myrates)  # give true values of initial estimates to get rid of potential optim problems
    all_rates <- rbind(all_rates,  c(rate_multiplier, myrates, 2 * rates$rates)) # ad-hoc factor 2 correction
  }
}
all_rates <- data.frame(all_rates)
names(all_rates) <- c("x", "true1", "true2", "est1", "est2")
medians <- ddply(all_rates, .(x), summarise, mtrue1 = median(log10(true1)), mtrue2 = median(log10(true2)), mest1 = median(log10(est1)), mest2 = median(log10(est2)))
lims <- c(-1.2, 0)
pdf("~/ownCloud/coronavirus/phylodynamique/Estimated_true_rates.pdf", width = 12, height = 6)
par(mfrow = c(1,2), mar = c(4,4,1,1))
plot(log10(all_rates$true1), log10(all_rates$est1), pch = 20, cex = 0.5, type = "p", las = 1,
     xlim = lims, ylim = lims, xlab = "True log10 rate", ylab = "Estimated log10 rate")
points((medians$mtrue1), (medians$mest1), pch = 20, cex = 2, col = "red")
points((medians$mtrue1), (medians$mtrue1), pch = 4, cex = 2, col = "red")
abline(0,1)
plot(log10(all_rates$true2), log10(all_rates$est2), pch = 20, cex = 0.5, type = "p", las = 1,
     xlim = lims, ylim = lims, xlab = "True log10 rate", ylab = "Estimated log10 rate")
points((medians$mtrue2), (medians$mest2), pch = 20, cex = 2, col = "red")
points((medians$mtrue2), (medians$mtrue2), pch = 4, cex = 2, col = "red")
abline(0,1)
dev.off()

#################################### SECOND SET OF SIMULATIONS WITH RANDOM SAMPLES  #################################### 

# now simulate  al large tree
tr10000 <- rtree(n = 10000, rooted = T) # random tree

if(data_tree <- TRUE){
  outgroup_label <- "EPI_ISL_402125"
  tr10000 <- read.tree(file = "~/ownCloud/coronavirus/phylodynamique/tree_gisaid_2020.tree") #... or SARS-CoV-2 tree
  tr10000 <- root(tr10000, outgroup = outgroup_label)
  # make it a tree of size 10,000:
  todrop <- sample(tr10000$tip.label[tr10000$tip.label!=outgroup_label], size = length(tr10000$tip.label)-10000, replace = F)
  tr10000 <- drop.tip(phy = tr10000, tip = todrop, rooted = T)
  # tr10000$edge.length[tr10000$edge.length==0] <-0.001
  tr10000 <- root(tr10000, outgroup = outgroup_label, resolve.root = T)
  tr10000 <- multi2di(tr10000)
}

stopifnot(length(tr10000$tip.label)==10000) # check this is a 10000-tips tree
stopifnot(is.rooted(tr10000))

# get rtt distance:
get_rtt_distance <- function(mytip, mytree = tr10000){
  #print(mytip)
  rtt <- 0; flag <- T
  while(flag){
    idx <- which(mytree$edge[, 2] == mytip)
    if(length(idx) == 0) flag <- F else {rtt <- rtt + mytree$edge.length[idx]; mytip <- mytree$edge[idx, 1]}
  }
  stopifnot(mytip == length(mytree$tip.label) + 1)
  return(rtt)
}
median_length <- median(sapply(1:length(tr10000$tip.label), get_rtt_distance)) # median root to tip distances
rate_multiplier <- 1 # set this to 2
myrates <- rate_multiplier * c(1/median_length, 2/median_length)

rates_boot <- c()

for(i in 1:nrep){
  print(i)
  # simulate trait:
  z <- rTraitDisc(phy = tr10000, model = "ARD", k = 2, rate = myrates, root.value = sample(x = 1:2, size = 1)) # phenotype
  for(k in 1:nrep){ # 50 bootstraps of size 100
    todrop <- sample(tr10000$tip.label, size = 10000-100, replace = F)
    sampled_tr <- drop.tip(phy = tr10000, tip = todrop) # sampled tree of size 100
    sampled_z <- z[sampled_tr$tip.label]
    stopifnot(all(sampled_tr$tip.label==names(sampled_z)))
    rates <- ace(x = sampled_z, phy = sampled_tr, type = "discrete", method = "ML", model = "ARD", ip = myrates)  # give true values of initial estimates to get rid of potential optim problems
    rates_boot <- rbind(rates_boot,  c(i, k, myrates, 2 * rates$rates)) # ad-hoc factor 2 correction
  }
}
  
rates_boot <- data.frame(rates_boot)
names(rates_boot) <- c("i", "k", "true1", "true2", "est1", "est2")
median_boot <- ddply(rates_boot, .(i), summarise, mtrue1 = median(log10(true1)), mtrue2 = median(log10(true2)), mest1 = median(log10(est1)), mest2 = median(log10(est2)))

# redraw previous figures
par(mfrow = c(1,1), mar = c(4,4,1,1))
truerate2_full <- unique(log10(all_rates$true2[all_rates$x==1])) # true rate in first simulations
truerate2_boot <- unique(log10(rates_boot$true2)) # true rate in bootstrap simulations
xmin <- 0
xmax <- nrep+1
pdf(paste0("~/ownCloud/coronavirus/phylodynamique/Estimated_true_rates_10000_draws_datatree_", data_tree, ".pdf"), width = 6, height = 6)
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(NULL, xlim = c(xmin, xmax), ylim = c(-1, 1), pch = 20, cex = 0.5, type = "p", las = 1, xlab = "Index of replicate", ylab = "Estimated log10 rate")
# add inference done on bs
colset <- rep(RColorBrewer::brewer.pal(n = 10, name = "Paired"), nrep/10)
for(i in 1:nrep){
  points(log10(rates_boot$true2[rates_boot$i==i])+1 * i, log10(rates_boot$est2[rates_boot$i==i]), pch = 20, cex = 0.5, type = "p", col = colset[i])
  points((median_boot$mtrue2[median_boot$i==i])+1 * i, (median_boot$mest2[median_boot$i==i]), pch = 20, cex = 2, type = "p", col = colset[i])
}
segments(x0 = xmin, x1 = xmax, y0 = median(median_boot$mest2), y1 = median(median_boot$mest2), lwd = 3, lty = 1, col = "red")
segments(x0 = xmin, x1 = xmax, y0 = truerate2_boot, y1 = truerate2_boot, lwd = 3, lty = 2, col = "gray")
legend("topleft", legend = c("Median of draws for one replicate", "Overall median across replicates", "True rate"), lty = c(NA, 1, 2), pch = c(20, NA, NA), lwd = 2, col  = c("gray", "red", "gray"), bty = "n")
dev.off()
 
# evaluate the behaviour of the error as a function of how many replicates we average

# 1) For the true replicated evolutionary history done on the tree of 1000

truc <- all_rates[all_rates$x==1, ] # work on inference for just one parameter value
truevalue <- truc$true1[1]
err <- c()
for(imax in 1:10){ # N = 1 to 10 replicates considered (among the nrep)
  for(j in 1:50){
    samp <- sample(1:nrep, size = imax, replace = F) # choose imax among the nrep replicates
    median_est <- median(truc$est1[samp])
    err <- rbind(err, c(imax, (median_est - truevalue)^2))
  }
}
err <- data.frame(err); names(err) <- c("imax", "error")
err <- ddply(err, .(imax), summarise, median.error = median(error))

pdf("~/ownCloud/coronavirus/phylodynamique/error_replicates_simuls1.pdf", width = 4, height = 3*3)
ymax <- 0.0008
par(mar = c(4,5,1,1), mfrow = c(3, 1))
plot(err$imax, err$median.error, type = "o", pch = 20, las = 1, xlab = "N", ylab = "Median error", ylim = c(0, ymax), bty = "n")
abline(h = 0, lty = 2)
plot(err$median.error * err$imax, type = "o", pch = 20, las = 1, xlab = "N", ylim = c(0, ymax), ylab = "Median error x N", bty = "n")
abline(h = 0, lty = 2)
plot(log(err$imax), log(err$median.error), type = "o", pch = 20, las = 1, xlab = "log(N)", ylab = "log(Median error)", bty = "n")
abline(log(err$median.error), -1, lty = 2)
dev.off()

# 2) For the bootstrapped trees of 100

truevalue <- rates_boot$true1[1]
err <- c()
for(i in 1:50){ # for each replicate of evolutionary history
  truc <- rates_boot[rates_boot$i==i, ] # select this replicate
  for(imax in 1:10){ # N = 1 to 10 draws considered (among the nrep)
    for(j in 1:50){ # we do several choices of which draws we sample, to average the error
      samp <- sample(1:nrep, size = imax, replace = F) # choose imax among the nrep replicates
      median_est <- median(truc$est1[samp])
      err <- rbind(err, c(i, imax, (median_est - truevalue)^2))
    }
  }
}

err <- data.frame(err); names(err) <- c("i", "imax", "error")
err <- ddply(err, .(i, imax), summarise, median.error = median(error)) # mean error when chosing imax replicates

pdf(paste0("~/ownCloud/coronavirus/phylodynamique/error_replicates_simuls2_datatree_", data_tree, ".pdf"), width = 4, height = 4*1)
par(mar = c(4,4,1,1), mfrow = c(1,1))

plot(NULL, xlab = "log(N)", ylab = "log(Median error)", xlim = c(0, 2.5), ylim = c(-12, -5), las = 1, bty = "n")
for(i in 1:50) points(log(err$imax[err$i==i]), log(err$median.error[err$i==i]), type = "o", pch = 20, las = 1, col = colset[i])

lm0 <- lm(log(median.error) ~ log(imax) + as.factor(i), data = err)
abline(lm0$coefficients["(Intercept)"] + mean(c(0, lm0$coefficients[3:51])), lm0$coefficients["log(imax)"], lwd = 2, lty = 2)

dev.off()


















#################################### AN ATTEMPT WITH THREE CHARACTERS  #################################### 

tr <- rtree(n = 200, rooted = T)
median_length <- median(diag(vcv.phylo(tr))) # root to tip distances
cols <- c("blue", "red", "green"); names(cols) <- c("A", "B", "C")
all_rates <- c()
for(rate_multiplier in c(0.5, 0.75, 1, 1.5, 2)){
  print(rate_multiplier)
  myrates <- rate_multiplier * c(1/median_length, 2/median_length, 3/median_length, 4/median_length, 5/median_length, 6/median_length)
  
  for(i in 1:10){
    #print(i)
    z <- rTraitDisc(phy = tr, model = "ARD", k = 3, rate = myrates, root.value = sample(x = 1:3, size = 1)) # phenotype
    #print(table(z))
    #plot.phylo(x = tr, tip.color = cols[z])
    rates <- ace(x = z, phy = tr, type = "discrete", method = "ML", model = "ARD", ip = myrates)  # give true values of initial estimates to get rid of potential optim problems
    all_rates <- rbind(all_rates,  c(rate_multiplier, myrates, 3 * rates$rates)) # ad-hoc correction
  }
}
all_rates <- data.frame(all_rates)
names(all_rates) <- c("x", paste0("true",1:6), paste0("est", 1:6))
medians <- ddply(all_rates, .(x), summarise, mtrue1 = median(true1), mtrue2 = median(true2), mest1 = median(est1), mest2 = median(est2))
par(mfrow = c(2,1), mar = c(4,4,1,1))
plot(log10(all_rates$true1), log10(all_rates$est1), pch = 20, cex = 0.5, type = "p", las = 1, xlab = "True", ylab = "Estimated")
points(log10(medians$mtrue1), log10(medians$mest1), pch = 20, cex = 2, col = "red")
abline(0,1)
plot(log10(all_rates$true2), log10(all_rates$est2), pch = 20, cex = 0.5, type = "p", las = 1, xlab = "True", ylab = "Estimated")
points(log10(medians$mtrue2), log10(medians$mest2), pch = 20, cex = 2, col = "red")
abline(0,1)
