# Script that computes the optimal threshold (according to Youden's J) for a locus being significant
# For a range of per-locus coverages (e.g. 10 to 200), the script simulates various configurations of
# the 4 bases (ACGT) for 4 cases: homozygous germline, heterozygous germline, homzygous+somatic,
# heterozygous + somatic. We then compute the optimal thresholds, such that as many as possible of the
# germline loci are eliminated and as many as possible of the somatic ones are kept

### Euler number
e <- exp(1)
hetero_prior = 0.0005
mut_prior = 1e-6
homo_prior = 1 - hetero_prior - mut_prior
factorials = NULL

comp_fact <- function() {
  fact = rep(1, 170)
  for(i in 2:170) {
    fact[i] = i * fact[i-1]
  }
  for(i in 1:170) {
    fact[i] = log(fact[i])
  }
  return (fact)
}

### function to compute log factorial using Stirling formula
log_fact <- function(n) {
  if (n < 2) return (0)
  if(n>=170)
  {
    return(0.5*log(2*pi*n)+n*log(n/e))
  }
  else return(factorials[n])
}


is_significant <- function(cov, log, res) {
  # decide homo/hetero
  p_homo <- res[4]*log.oneMinusTheta + (cov-res[4])*log.theta + log.homoPrior + log.priorH0
  multinomCoef <- log_fact(cov) - log_fact(res[1]) - log_fact(res[2]) - log_fact(res[3]) - log_fact(res[4])
  ### normalization coefficient (probability of evidence)
  # 1. All true allelles are c1
  prob_all_c1 = homo_prior*((1-theta)^res[4])*(theta/3)^(cov-res[4])
  # 2. The locus is heterozygous c1 c2
  prob_hetero = hetero_prior * ((0.5-theta/3)^(res[3]+res[4])) * ((theta/3)^(res[1]+res[2]))
  # 3 The locus is homozygous + somatic mutation
  prob_homo_som = homo_prior * mut_prior * ((0.75-0.66*theta)^res[4])*(0.25^res[3])*((theta/3)^(res[1]+res[2]))
  # 4 The locus is heterozygous + somatic mutation
  prob_hetero_som = hetero_prior * mut_prior * ((0.5-theta)^res[4])*((0.25)^(res[2]+res[3])) * ((theta/3)^res[1])
  # 5. Two somatic mutations
  prob_two_somatic = hetero_prior*mut_prior^2*(1-theta)^cov
  normCoeff = log(prob_all_c1 + prob_hetero + prob_homo_som + prob_hetero_som + prob_two_somatic)
  # message("p_homo ", p_homo)
  # message("Evidence: ", normCoeff)
  # message("Probability homozygous: ", exp(p_homo+multinomCoef))
  # message("Prob all ", prob_all_c1, " Prob hetero ", prob_hetero, " Prob homo som: ", prob_homo_som, " Prob hetero som: ", prob_hetero_som, " Probability: ", exp(p_homo-normCoeff))
  return (p_homo - normCoeff)
}

is_significant2 <- function(cov, log, res) {
  # decide homo/hetero
  p_homo <- res[4]*log.oneMinusTheta + (cov-res[4])*log.theta + log.homoPrior + log.priorH0
  multinomCoef <- log_fact(cov) - log_fact(res[1]) - log_fact(res[2]) - log_fact(res[3]) - log_fact(res[4])
  ### normalization coefficient
  normCoeff <- log_fact(cov+3) - log_fact(cov) - log(6)

  return (p_homo + normCoeff + multinomCoef)
}

set.seed(123)
### Always do 100000 simulations.
nsim <- 10000
### theta (error rate)
theta <- 0.001
log.theta <- log(theta/3)
log.oneMinusTheta <- log(1-theta)
### homo prior
log.homoPrior <- log(1/4)
### number of cells
ncells <- 140

### prior on H_0 (no somatic mutation)
log.priorH0 <- log(0.998)
### data frame to store the results
data_all <- data.frame()

factorials = comp_fact()

for (meanCov in seq(10,200,10)) {
  message("\nCoverage: ", meanCov)
  cat("Proportion of tumor cells: ")
  for (prop.tumour in seq(0.1,0.9,0.1)) {
    cat(prop.tumour, " ")
    prop.healthy <- 1-prop.tumour

    ### Homozygous germline
    logps.germ.homo <- rep(0,nsim)
    for (i in 1:nsim) {
      # sample coverage
      cov <- rpois(1, meanCov)

      # generate the resulting base counts
      rand <- runif(cov)
      n1 <- sum(rand<theta/3)
      n2 <- sum(rand<2*theta/3)-n1
      n3 <- sum(rand<theta)-n1-n2
      if (n1 == 0 & n2 == 0 & n3 == 0) {
        n1 <- 1
      }
      n4 <- cov - n1 - n2 - n3
      res <- sort(c(n1,n2,n3,n4))

      logps.germ.homo[i] <- is_significant(cov, log, res)
    }

    ### Heterozygous germline
    logps.germ.hetero <- rep(0,nsim)
    for (i in 1:nsim) {
      # sample coverage
      cov <- rpois(1, meanCov)

      # generate the resulting base counts
      rand <- runif(cov)
      n1 <- sum(rand<theta/3)
      n2 <- sum(rand<2*theta/3)-n1
      if (n1 == 0 & n2 == 0) {
        n1 <- 1
      }
      n3 <- sum(rand<0.5+theta/3)-n1-n2
      n4 <- cov - n1 - n2 - n3
      res <- sort(c(n1,n2,n3,n4))

      logps.germ.hetero[i] <- is_significant(cov, log, res)
    }

    ### Homozygous somatic mutation
    logps.som.homo <- rep(0,nsim)
    for (i in 1:nsim) {
      # sample coverage
      cov <- rpois(1, meanCov)

      # generate the resulting base counts
      # without error there is prop.tumour of one allele and prop.healthy of another allele
      # with error:
      # "healthy" allele: f_h(1-theta) + f_t*theta/3
      # "tumour" allele: f_t(1-theta) + f_h*theta/3
      # the other two allels: each theta/3
      # f_h and f_t are the proportion of reads coming from healthy and tumour cells
      rand <- runif(cov)
      cov.healthy <- cov*prop.healthy
      cov.tumour <- cov*prop.tumour
      breaks <- c(theta/3,
        theta/3,
        cov.healthy*(1-theta)/cov+cov.tumour*theta/cov/3,
        cov.tumour*(1-theta)/cov+cov.healthy*theta/cov/3)
      cum.sum <- cumsum(breaks)
      n1 <- sum(rand<cum.sum[1])
      n2 <- sum(rand<cum.sum[2])-n1
      n3 <- sum(rand<cum.sum[3])-n1-n2
      n4 <- cov - n1 - n2 - n3
      res <- sort(c(n1,n2,n3,n4))

      logps.som.homo[i] <- is_significant(cov, log, res)
    }

    ### Heterozygous somatic mutation
    # there are two options: healthy AA, tumour AB; or healthy AB, tumour AA (loss of heterozygosity)
    # the first possibility
    logps.som.hetero.1 <- rep(0,nsim)
    for (i in 1:nsim) {
      # sample coverage
      cov <- rpois(1, meanCov)

      # generate the resulting base counts
      # without error there is prop.healthy + 1/2prop.tumour of one allele and 1/2prop.tumour of another allele
      # with error:
      # "healthy" allele: (f_h+0.5f_t)(1-theta) + 0.5f_t*theta/3
      # "tumour" allele: 0.5f_t(1-theta) + (f_h+0.5f_t)*theta/3
      # the other two allels: each theta/3
      # f_h and f_t are the proportion of reads coming from healthy and tumour cells
      rand <- runif(cov)
      cov.healthy <- cov*prop.healthy
      cov.tumour <- cov*prop.tumour
      breaks <- c(theta/3,
        theta/3,
        (cov.healthy+0.5*cov.tumour)*(1-theta)/cov+0.5*cov.tumour*theta/cov/3,
        0.5*cov.tumour*(1-theta)/cov+(cov.healthy+0.5*cov.tumour)*theta/cov/3)
      cum.sum <- cumsum(breaks)
      n1 <- sum(rand<cum.sum[1])
      n2 <- sum(rand<cum.sum[2])-n1
      n3 <- sum(rand<cum.sum[3])-n1-n2
      n4 <- cov - n1 - n2 - n3
      res <- sort(c(n1,n2,n3,n4))

      logps.som.hetero.1[i] <- is_significant(cov, log, res)
    }

    # the second possibility
    logps.som.hetero.2 <- rep(0,nsim)
    for (i in 1:nsim) {
      # sample coverage
      cov <- rpois(1, meanCov)

      # generate the resulting base counts
      # without error there is 0.5*prop.healthy + prop.tumour of one allele and 1/2prop.healthy of another allele
      # with error:
      # allele A: (0.5*f_h+f_t)(1-theta) + 0.5f_h*theta/3
      # allele B: 0.5f_h(1-theta) + (0.5*f_h+f_t)*theta/3
      # the other two allels: each theta/3
      # f_h and f_t are the proportion of reads coming from healthy and tumour cells
      rand <- runif(cov)
      cov.healthy <- cov*prop.healthy
      cov.tumour <- cov*prop.tumour
      breaks <- c(theta/3,
        theta/3,
        (0.5*cov.healthy+cov.tumour)*(1-theta)/cov+0.5*cov.healthy*theta/cov/3,
        0.5*cov.healthy*(1-theta)/cov+(0.5*cov.healthy+cov.tumour)*theta/cov/3)
      cum.sum <- cumsum(breaks)
      n1 <- sum(rand<cum.sum[1])
      n2 <- sum(rand<cum.sum[2])-n1
      n3 <- sum(rand<cum.sum[3])-n1-n2
      n4 <- cov - n1 - n2 - n3
      res <- sort(c(n1,n2,n3,n4))

      logps.som.hetero.2[i] <- is_significant(cov, log, res)
    }


    data <- data.frame(logp = c(logps.germ.homo, logps.som.hetero.1, logps.som.hetero.2),
     type = c(rep("homozygous germline", nsim),
      #rep("heterozygous germline", nsim),
                                #rep("homozygous mutation", nsim),
                                rep("heterozygous mutation AA+AB", nsim),
                                rep("heterozygous mutation AB+AA", nsim)),
     meanCov = rep(meanCov, 3*nsim),
     prop.tumour = rep(prop.tumour, 3*nsim))
    data_all <- rbind(data_all, data)
  }
}

### save data_all
message("\nWriting all probabilities to Ks_all.csv...")
write.table(data_all, "Ks_all.csv",quote=FALSE,sep=",",row.names=FALSE)

message("Computing optimal Ks...")

###### find absolutely best K, parameters consistent with bayes optimal
pos.all <- 3e9
heteroPrior <- 0.001
pos.hetero <- pos.all*heteroPrior
somMutPrior <- 10e-5
pos.mut <- pos.all*somMutPrior
pos.homo <- pos.all-pos.hetero-pos.mut

best.Ks <- data.frame(meanCov=rep(seq(10,200,10), each=9),
  prop.tumour=rep(seq(0.1,0.9,0.1),20),
  bestK.J=rep(0,180),
  kept.pos.J=rep(0,180),
  kept.mut.J=rep(0,180),
  J=rep(0,180)
  )
index=1
for (meanCov in seq(10,200,10)) {
  message("\nCoverage: ", meanCov)
  cat("Proportion of tumor cells: ")
  for (prop.tumour in seq(0.1,0.9,0.1)) {
   cat(prop.tumour, " ")

   data.tmp <- data_all[data_all$meanCov==meanCov & abs(data_all$prop.tumour-prop.tumour)<0.09,]
   data.tmp$mut <- ifelse(data.tmp$type %in% c("homozygous mutation", "heterozygous mutation AA+AB", "heterozygous mutation AB+AA"),
     "yes",
     "no")
   data.tmp$hetero <- ifelse(data.tmp$type=="heterozygous germline",
    "yes",
    "no")
   data.tmp <- data.tmp[order(data.tmp$logp),]

   best.K.J <- -Inf
   best.kept.pos.J <- -1
   best.kept.mut.J <- -1
   best.J <- -2

   kept.mut <- 0
   kept.pos.homo <- 0
   kept.pos.hetero <- 0

   cutoffs <- unique(data.tmp$logp)

  # for each cut-off except the last one -- the largest one does nothing anyway (we just keep everything)
  for (c in cutoffs[1:(length(cutoffs)-1)]) {
    i <- data.tmp$logp==c
    kept.mut <- kept.mut+sum(data.tmp$mut[i]=="yes")
    kept.pos.hetero <- kept.pos.hetero + sum(data.tmp$hetero[i]=="yes")
    kept.pos.homo <- kept.pos.homo + sum(data.tmp$mut[i]=="no" & data.tmp$hetero[i]=="no")

    tp <- kept.mut*pos.mut/(2*nsim)
    fp <- kept.pos.homo*pos.homo/nsim + kept.pos.hetero*pos.hetero/nsim
    tn <- (pos.homo + pos.hetero)-fp
    fn <- pos.mut-tp

    mutations <- tp/(tp+fn)
    positions <- ((1-heteroPrior)*(1-somMutPrior)*kept.pos.homo/nsim + heteroPrior*(1-somMutPrior)*kept.pos.hetero/nsim)

    J <- tp/(tp+fn) + tn/(tn+fp) -1
    if(J >= best.J) {
      best.J <- J
      best.kept.mut.J <- mutations
      best.kept.pos.J <- positions
      best.K.J <- data.tmp$logp[i][1]
    }
  }
  best.Ks[index,] <- c(meanCov, prop.tumour, best.K.J, best.kept.pos.J,
   best.kept.mut.J, best.J)
  index <- index+1
  }}

  write.csv(best.Ks, "Ks.csv", quote=FALSE, row.names=FALSE)




