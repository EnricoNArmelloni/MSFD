virtualPop_new3 <- function (tincr = 1/12, K.mu = 0.5, K.cv = 0.1, Linf.mu = 80, 
                           Linf.cv = 0.1, ts = 0.5, C = 0.75, LWa = 0.01, LWb = 3, Lmat = 40, 
                           wmat = 8, rmax = 10000, beta = 1, 
                           repro_wt = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), 
                           M = 0.7, harvest_rate = M, 
                           L50 = 0.25 * Linf.mu,
                           wqs = L50 * 0.2,
                           R50 = 35, al = 0.7 , ar = 0.7, pmax = 1, pmin = 0.3, bin.size = 1, timemin = 0, 
                           timemax = 10, timemin.date = as.Date("1980-01-01"), N0 = 5000, 
                           fished_t = seq(timemin + 5, timemax, tincr), lfqFrac = 1,
                           age = c(0,1,2,3,4,5), mortality = c(1,0.8,0.6,0.4,0.2,0.1),
                           Lopt = 2/3*Linf.mu,
                           L_cut = Linf.mu,
                           selectivity = "logistic",
                           progressBar = TRUE) 
{
  require(dublogistic)
  require(tidyverse)
  require(magrittr)
  require(fishdynr)
  #source("s/inverse_vonBert.R")
  
    timeseq <- seq(from = timemin, to = timemax, by = tincr)
  if (!zapsmall(1/tincr) == length(repro_wt)) 
    stop("length of repro_wt must equal the number of tincr in one year")
  repro_wt <- repro_wt/sum(repro_wt)
  repro_t <- rep(repro_wt, length = length(timeseq))
  lfq <- vector(mode = "list", length(timeseq))
  names(lfq) <- timeseq
  indsSamp <- vector(mode = "list", length(timeseq))
  names(indsSamp) <- timeseq
  tmaxrecr <- (which.max(repro_wt) - 1) * tincr
  phiprime.mu = log10(K.mu) + 2 * log10(Linf.mu)
  
  
  # calculate natural mortality
  age1 <- age
  mortality1 <- mortality
  increase <- c()
  decrease <- c()
  age2 <- c()
  mortality2 <- c()
  aa <- c()
  mm <- c()
  for (i in 1:length(age1)) {
    
    increase[i] <- (age1[i+1] - age1[i])/2
    age2[i] <- age1[i] + increase[i]
    decrease[i] <- (mortality1[i] - mortality1[i+1])/2  
    mortality2[i] <- mortality1[i] - decrease[i]
  }
  
  aa <- seq(max(age1) +1, max(age1 + 15),1)
  mm <- rep(min(mortality1), length(aa))
  
  age3 <- c(age1, aa, age2)
  mortality3 <- c(mortality1, mm, mortality2)
  
  df <- as.data.frame(cbind(age3, mortality3))
  df <- na.omit(df)
  lm <- lm(mortality3 ~ age3, data = df)
  require(segmented)
  control <- seg.control(fix.npsi = TRUE)
  my.seg <- segmented(lm, 
                      seg.Z = ~ age3,
                      psi = age1[2:(length(age1)-1)], control = control)
  mort <- ggplot() +
    geom_point(aes(x = df$age3 , y = my.seg$fitted.values)) +
    geom_line(aes(x = df$age3 , y = my.seg$fitted.values), color = "red")
  pdf(file = "maps/natural_mortality.pdf")
  print(mort)
  dev.off()
  
  
  make.inds <- function(id = NA, A = 0, L = 0, W = NA, mat = 0, 
                        K = NA, Winf = NA, Linf = NA, phiprime = NA, F = NA, 
                        Z = NA, Fd = 0, alive = 1, M = NA) {
    inds <- data.frame(id = id, A = A, L = L, W = W, Lmat = LWa * 
                         L^LWb, mat = mat, K = K, Linf = Linf, Winf = Winf, 
                       phiprime = phiprime, F = F, Z = Z, Fd = Fd, alive = alive, M = NA)
    lastID <<- max(inds$id)
    return(inds)
  }
  express.inds <- function(inds) {
    inds$Linf <- ifelse(is.na(inds$Linf), Linf.mu * rlnorm(nrow(inds), 
                                                           0, Linf.cv), inds$Linf)
    inds$K <- ifelse(is.na(inds$K), K.mu * rlnorm(nrow(inds), 
                                                  0, K.cv), inds$K)
    inds$Winf <- LWa * inds$Linf^LWb
    inds$W <- LWa * inds$L^LWb
    inds$phiprime <- log10(inds$K) + 2 * log10(inds$Linf)
    inds$Lmat <- rnorm(nrow(inds), mean = Lmat, sd = wmat/diff(qnorm(c(0.25, 
                                                                       0.75))))
    return(inds)
  }
  grow.inds <- function(inds) {
    L2 <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, 
                         C = C, L1 = inds$L, t1 = t - tincr, t2 = t)
    inds$L <- L2
    inds$W <- LWa * inds$L^LWb
    inds$A <- inds$A + tincr
    return(inds)
  }
  mature.inds <- function(inds) {
    inds$mat <- ifelse((inds$L > inds$Lmat | inds$mat == 
                          1), 1, 0)
    return(inds)
  }
  death.inds <- function(inds) {
    
    if (length(age) > 1) {
      
      inds$M <- predict(my.seg, newdata = data.frame(age3 = inds$A))
      
    } else {
      
      inds$M <- M
    }
    
    
    # choose net selectivity
    if (selectivity == "logistic") {
      pSel <- logisticSelect(inds$L, L50, wqs)
      inds$F <- pSel * Fmax
      
      # sel_plot <- ggplot() +
      #   geom_point(aes(x = inds$L, y = pSel)) +
      #   geom_line(aes(x = inds$L, y = pSel), color = "red")
      # print(sel_plot)
      
      
    } else if (selectivity == "dublogistic") { 
      pSel <- dublogistic.f(inds$L, L50, R50, al, ar, pmax, pmin)
      inds$F <- pSel$selectivity * Fmax
    }
    
    inds$Z <- inds$M + inds$F
    pDeath <- 1 - exp(-inds$Z * tincr)
    dead <- which(runif(nrow(inds)) < pDeath)
    if (length(dead) > 0) {
      inds$alive[dead] <- 0
      tmp <- cbind(inds$F[dead], inds$Z[dead],inds$M[dead])
      Fd <- apply(tmp, 1, FUN = function(x) {
        sample(c(0, 1), size = 1, prob = c(x[3]/x[2], x[1]/x[2]))
      })
      inds$Fd[dead] <- Fd
      rm(tmp)
    }
    
    return(inds)
  }
  remove.inds <- function(inds) {
    dead <- which(inds$alive == 0)
    if (length(dead) > 0) {
      inds <- inds[-dead, ]
    }
    return(inds)
  }
  reproduce.inds <- function(inds) {
    if (repro > 0 & sum(inds$mat) > 0) {
      SSB <- sum(inds$W * inds$mat)
      n.recruits <- ceiling(srrBH(rmax, beta, SSB) * repro)
      offspring <- make.inds(id = seq(lastID + 1, length.out = n.recruits))
      offspring <- express.inds(offspring)
      inds <- rbind(inds, offspring)
    }
    return(inds)
  }
  record.inds <- function(inds, ids = 1:10, rec = NULL) {
    if (is.null(rec)) {
      rec <- vector(mode = "list", length(ids))
      names(rec) <- ids
      inds <- inds
    }
    else {
      ids <- as.numeric(names(rec))
    }
    if (length(rec) > 0) {
      inds.rows.rec <- which(!is.na(match(inds$id, ids)))
      if (length(inds.rows.rec) > 0) {
        for (ii in inds.rows.rec) {
          match.id <- match(inds$id[ii], ids)
          if (is.null(rec[[match.id]])) {
            rec[[match.id]] <- inds[ii, ]
          }
          else {
            rec[[match.id]] <- rbind(rec[[match.id]], 
                                     inds[ii, ])
          }
        }
      }
    }
    rec
  }
  lastID <- 0
  inds <- make.inds(id = seq(N0))
  inds <- express.inds(inds)
  res <- list()
  res$sea_pop <- list()
  res$pop <- list(dates = yeardec2date(date2yeardec(timemin.date) + 
                                         (timeseq - timemin)), N = NaN * timeseq, B = NaN * timeseq, 
                  SSB = NaN * timeseq, L95 = NaN * timeseq, L95_c = NaN *timeseq,
                  L95_p = NaN * timeseq, L50mat_p = NaN * timeseq, Lopt_p = NaN * timeseq, 
                  L95_pc = NaN * timeseq, L50mat_pc = NaN * timeseq, Lopt_pc = NaN * timeseq)
  if (progressBar) 
    pb <- txtProgressBar(min = 1, max = length(timeseq), 
                         style = 3)

  for (j in seq(timeseq)) { 
    t <- timeseq[j]
    if (length(fished_t) == 0) {
      Fmax <- 0
      lfqSamp <- 0
    }
    else {
      if (min(sqrt((t - fished_t)^2)) < 1e-08) {
        Fmax <- harvest_rate
        lfqSamp <- 1
      } else  if (harvest_rate == 0) {
          Fmax <- 0
          lfqSamp <- 0
      }
      
      else {
        Fmax <- 0
        lfqSamp <- 0
      }
    }
    repro <- repro_t[j]
    inds <- grow.inds(inds)
    inds <- mature.inds(inds)
    inds <- reproduce.inds(inds)
    inds <- death.inds(inds)
    
    
    if (lfqSamp) {
      samp <- try(sample(seq(inds$L), ceiling(sum(inds$Fd) * 
                                                lfqFrac), prob = inds$Fd), silent = TRUE)
      if (class(samp) != "try-error") {
        lfq[[j]] <- inds$L[samp]
        indsSamp[[j]] <- inds[samp, ]
      }
      rm(samp)
    }
    
    inds <- remove.inds(inds)
    res$pop$N[j] <- nrow(inds)
    res$pop$B[j] <- sum(inds$W)
    res$pop$SSB[j] <- sum(inds$W * inds$mat)
    res$pop$L95[j] <- quantile(inds$L, probs = 0.95)
    if (harvest_rate != 0) {
      res$pop$L95_c[j] <- quantile(indsSamp[[j]]$L, probs = 0.95)
    }
    
    
    # calculate proportions for sea pop
    L95 <- inds %>% count(L >= res$pop$L95[j])
    res$pop$L95_p[j] <- (L95[2,2] / res$pop$N[j])*100
    L50mat <- inds %>% count(L >= Lmat)
    res$pop$L50mat_p[j] <- (L50mat[2,2] / res$pop$N[j])*100
    Lopt_p <- inds %>% count(L >= Lopt)
    res$pop$Lopt_p[j] <- (Lopt_p[2,2] / res$pop$N[j])*100
    
    # calculate proportions for catches
    if(lfqSamp == 1 & harvest_rate != 0) {
      catch <- indsSamp[[j]]
      L95_c <- catch %>% count(L >= res$pop$L95[j])
      res$pop$L95_pc[j] <- (L95[2,2] / res$pop$N[j])*100
      L50mat_c <- catch %>% count(L >= Lmat)
      res$pop$L50mat_pc[j] <- (L50mat[2,2] / res$pop$N[j])*100
      Lopt_c <- catch %>% count(L >= Lopt)
      res$pop$Lopt_pc[j] <- (Lopt_c[2,2] / res$pop$N[j])*100
    }
    
    # keep the last 12 months
    if (j > length(timeseq)-12) {
      #aggregate data
      bb <- seq(0,max(inds$L)+5,1)
      res$sea_pop[[j]] <- inds
      res$sea_pop[[j]] %<>% 
        mutate(size_cl = cut(L, breaks = bb, include.lowest = TRUE)) %>% 
        group_by(size_cl) %>% 
        summarise(biomass = sum(W),
                  n = n()) %>%
        ungroup() %>% 
        complete(size_cl) %>% 
        mutate(biomass = ifelse(is.na(biomass), 0, biomass),
               n = ifelse(is.na(n), 0, n),
               length = row_number()) %>% 
        filter(length <= L_cut)
    }
    
    
    if (progressBar) 
      setTxtProgressBar(pb, j)
  }
  if (progressBar) 
    close(pb)
  if (harvest_rate != 0) {
    lfq2 <- lfq[which(sapply(lfq, length) > 0)]
    dates <- yeardec2date(date2yeardec(timemin.date) + (as.numeric(names(lfq2)) - 
                                                          timemin))
    Lran <- range(unlist(lfq2))
    Lran[1] <- floor(Lran[1])
    Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 
                  1) * bin.size
    bin.breaks <- seq(Lran[1], Lran[2], by = bin.size)
    bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
    res$lfqbin <- list(sample.no = seq(bin.mids), midLengths = bin.mids, 
                       dates = dates, catch = sapply(lfq2, FUN = function(x) {
                         hist(x, breaks = bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
                       }))
    indsSamp <- indsSamp[which(sapply(indsSamp, length) > 0)]
    res$inds <- indsSamp
    res$growthpars <- list(K = K.mu, Linf = Linf.mu, C = C, ts = ts, 
                           phiprime = phiprime.mu, tmaxrecr = tmaxrecr)
    
    
    # aggregate data for captures
    for (i in 1:length(res$inds)) {
      
      cc <- seq(0,max(res$inds[[i]]$L)+5,1)
      res$inds[[i]] %<>% mutate(size_cl = cut(L, breaks = cc, include.lowest = TRUE)) %>% 
        group_by(size_cl) %>% 
        summarise(biomass = sum(W),
                  n = n()) %>%
        ungroup() %>% 
        complete(size_cl) %>% 
        mutate(biomass = ifelse(is.na(biomass), 0, biomass),
               n = ifelse(is.na(n), 0, n),
               length = row_number()) %>% 
        filter(length <= L_cut)
    }
  }
  # erase empty years
  res$sea_pop[sapply(res$sea_pop, is.null)] <- NULL
  

  
  return(res)
}



