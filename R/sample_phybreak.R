#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak} object (\code{sample.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
#' @param classic The proportion of tree updates with the classic protocol (reference see below), in which within-host
#'   minitrees are proposed by simulating coalescent times and tree topology. In the current default protocol only 
#'   coalescent times are proposed with the minitree topology kept intact. This is followed by removing and reconnecting
#'   the sampling tips one by one. This results in better mixing of the mcmc-chain if there is much 
#'   genetic information (many SNPs) and/or if there are many possible within-host minitree topologies 
#'   (e.g. many samples per host). The classic protocol is faster in terms of updates/second and can thus be more efficient
#'   with little genetic information.
#' @param keepphylo The proportion of tree updates keeping the phylogenetic tree intact, only possible if there is one 
#'   sample per host and the \code{wh.model = "linear"} with complete bottleneck.
#' @param withinhost_only The proportion of tree updates in which only the within-host minitree is sampled, and 
#'   the transmission tree and infection times are kept unchanged.
#' @param parameter_frequency The relative frequency by which the model parameters are updated relative to updating each host.
#' @param status_interval The number of seconds between each on-screen print of the progress of the mcmc-chain.
#' @return The \code{phybreak} object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample_phybreak <- function(x, nsample, thin = 1, thinswap = 1, classic = 0, keepphylo = 0, withinhost_only = 0, 
                            parameter_frequency = 1, status_interval = 10, 
                            verbose = 1, historydist = 0.5,
                            nchains = 1, heats = NULL, all_chains = FALSE, parallel = FALSE, shiny = FALSE, ...) {
  
  if (parallel)
    return(sample_phybreak_parallel(x, nsample, thin, thinswap, classic, keepphylo, withinhost_only, 
                                    parameter_frequency, status_interval, 
                                    verbose, historydist, nchains, heats, all_chains, ...)) 
      
  ### tests
  if(nsample < 1) stop("nsample should be positive")
  if(thin < 1) stop("thin should be positive")
  if(is.null(x$p$wh.bottleneck)) {
    x$p$wh.bottleneck <- choose_whbottleneck("auto", x$p$wh.model)
  }
  if(classic < 0 | classic > 1) stop("classic should be a fraction")
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(withinhost_only < 0 | withinhost_only > 1) stop("withinhost_only should be a fraction")
  if(withinhost_only + keepphylo + classic > 1) stop("classic + keepphylo + withinhost_only should be a fraction")
  if(keepphylo > 0) {
    if(any(duplicated(x$d$hostnames)) || !(x$p$wh.model %in% c(3, "linear")) || x$p$wh.bottleneck == "wide") {
      keepphylo <- 0
      warning("model incompatible with keepphylo-updates: they will not be used", immediate. = TRUE)
    } 
  }
  
  #le <- environment()
  x <- add_modules_to_phybreak(x, phyb.obj = FALSE, ...)
  
  if(!x$p$mult.intro & historydist > 0) {
    historydist <- 0
  }
  if (is.null(heats))
    heats <- 1/(1+1*(1:nchains-1))
  else if (inherits(heats, "numeric") & length(heats) != nchains)
    stop("length of heats is not the same as number of chains")
  else if (!inherits(heats, "numeric"))
    stop("heats is not a numeric vector")
  
  ### WHY NOT printmessage AND printlog AS INPUT?
  if(verbose == 1){
    printmessage = TRUE
    printlog = TRUE
  } else if (verbose == 0){
    printmessage = FALSE
    printlog = FALSE
  } else if (verbose == 2){
    printmessage = FALSE
    printlog = TRUE
  }
  
  ### add distance model if not present
  # if(is.null(x$p$dist.model)) {
  #   x$p$dist.model <- "none"
  #   x$p$dist.exponent <- 2
  #   x$p$dist.scale <- 1
  #   x$p$dist.mean <- 1
  #   x$s$dist.e <- c()
  #   x$s$dist.s <- c()
  #   x$s$dist.m <- c()
  #   x$h$est.dist.s <- FALSE
  #   x$h$est.dist.e <- FALSE
  #   x$h$est.dist.m <- FALSE
  # }
  
  protocoldistribution <- c(1 - classic - keepphylo - withinhost_only, classic, keepphylo, withinhost_only)
  historydistribution <- c(historydist, 1 - historydist)
  
  ### create room in s to add the new posterior samples
  s.post <- list(inftimes = with(x, cbind(s$inftimes, matrix(NA, nrow = p$obs, ncol = nsample))),
                 infectors = with(x, cbind(s$infectors, matrix(NA, nrow = p$obs, ncol = nsample))),
                 nodetimes = with(x, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                 nodehosts = with(x, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                 nodeparents = with(x, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples - 1, ncol = nsample))),
                 introductions = c(x$s$introductions, rep(NA, nsample)),
                 mu = c(x$s$mu, rep(NA, nsample)), 
                 mG = c(x$s$mG, rep(NA, nsample)), 
                 mS = c(x$s$mS, rep(NA, nsample)), 
                 wh.s = c(x$s$wh.s, rep(NA, nsample)), 
                 wh.e = c(x$s$wh.e, rep(NA, nsample)), 
                 wh.0 = c(x$s$wh.0, rep(NA, nsample)),
                 logLik = c(x$s$logLik, rep(NA, nsample)),
                 heat = c(x$s$heat, rep(NA, nsample)))
  
  ### create room for additional parameters from modules
  s.post.modules <- lapply(setdiff(names(x$s), names(s.post)), function(n){
    if(n == "cnt.types"){
      with(x, cbind(s$cnt.types, matrix(NA, nrow = 2, ncol = nsample - 1)))
    } else{
      return(c(x$s[[n]], rep(NA, nsample)))
    }
  })
  names(s.post.modules) <- setdiff(names(x$s), names(s.post))
  s.post <- c(s.post, s.post.modules)
    
                 # ir = c(x$s$ir, rep(NA, nsample)),
                 # wh.h = c(x$s$wh.h, rep(NA, nsample)), 
                 # dist.e = c(x$s$dist.e, rep(NA, nsample)), 
                 # dist.s = c(x$s$dist.s, rep(NA, nsample)), 
                 # dist.m = c(x$s$dist.m, rep(NA, nsample)), 
                 # R = c(x$s$R, rep(NA, nsample)))
  
  # if (!is.null(userenv$samplers)){
  #   for (n in names(userenv$samplers)){
  #     s.post <- c(s.post, list(c(x$s[[n]], rep(NA, nsample))))
  #     names(s.post)[length(s.post)] <- n
  #   }
  # }
  #     
  s.posts <- lapply(1:nchains, function(i) s.post)
    
  build_pbe(x)
  
  envirs <- list()
  npars <- sum(grepl("est", names(x$h)))
  
  for (n in 1:nchains){
    heat <- heats[n]
    copy2pbe0("heat", environment())
    chain <- n
    copy2pbe0("chain", environment())
    envirs[[n]] <- as.environment(as.list(pbe0, all.names = TRUE))
  }
  
  if (printmessage)   
    message(paste0("  sample      logLik  introductions       mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  if (printlog)
    print_screen_log(length(x$s$mu))
  
  curtime <- Sys.time()
  
  swap_thin <- 0
  shared_heats <- heats
  
  for (sa in tail(1:length(s.posts[[1]]$mu), nsample)) {
    #if(shiny == TRUE) incProgress(1/length(tail(1:length(s.posts[[1]]$mu), nsample)))
    for (rep in 1:thin) {
      if(Sys.time() - curtime > status_interval & printlog == TRUE) {
        #saveRDS(s.posts, "~/Documents/phybreak/test.rds")
        heats <- sapply(envirs, function(e) return(e[["heat"]]))
        for (i in ls(envir=envirs[[which(heats == 1)]])) copy2pbe0(i, envirs[[which(heats == 1)]])
        print_screen_log(sa)
        curtime <- Sys.time()
      }
      
      for (i in 1:nchains){
        envirs[[i]]$heat <- shared_heats[i]
      }
      
      envirs <- lapply (envirs, function(e) {
        for (i in ls(envir=e)) copy2pbe0(i, e)
        
        for(i in  sample(c(rep(-(1:npars), parameter_frequency), 0:x$p$obs))) {
          if(i >= 0) {
            which_protocol <- sample(c("edgewise", "classic", "keepphylo", "withinhost"),
                                     1,
                                     prob = protocoldistribution)
            history <- sample(c(TRUE, FALSE), 1, prob = historydistribution)
            update_host(i, which_protocol, history || i == 0)
          }
          
          if (i == -1 && x$h$est.mu) update_mu()
          if (i == -2 && x$h$est.mG) update_mG()
          if (i == -3 && x$h$est.mS) update_mS()
          if (i == -4 && x$h$est.wh.s)  update_wh_slope()
          if (i == -5 && x$h$est.wh.e)  update_wh_exponent()
          if (i == -6 && x$h$est.wh.0)  update_wh_level()
          if (i < -6)  x$updaters[[-6 - i]]()
          
          # if (i == -7 && x$h$est.wh.h) update_wh_history()
          # if (i == -8 && x$h$est.dist.e)  update_dist_exponent()
          # if (i == -9 && x$h$est.dist.s)  update_dist_scale()
          # if (i == -10 && x$h$est.dist.m)  update_dist_mean()
          # if (i == -11 && x$h$est.ir) update_ir()
          # if (i == -12 && x$h$est.R) update_R()
          # if (i < -12){
          #   if (!is.null(userenv$helpers)) {
          #     if (userenv$helpers[[-12 - i]]) {
          #       userenv$updaters[[-12 - i]]()
          #     }
          #   }
          # }
          
        }
        
        as.environment(as.list(pbe0, all.names = TRUE))
      })
      
      if(nchains > 1 & swap_thin %% thinswap == 0){
        shared_lik <- do.call(cbind, lapply(envirs, function(xx){
          sum(sapply(names(xx)[grepl("logLik", names(xx))], function(n) xx[[n]]))
          # sum(xx$logLikcoal, xx$logLikgen, xx$logLiksam, xx$logLikdist, 
          #     xx$logLikseq)
        }))
        
        shared_heats <- swap_heats(shared_heats, shared_lik)
      } 
      swap_thin <- swap_thin + 1
    }
    
    s.posts <- lapply(envirs, function(e){
      for (i in ls(envir=e)) copy2pbe0(i, e)
      
      chain <- pbe0$chain
      s.post <- s.posts[[chain]]
      
      vars_to_log <- environment2phybreak(pbe0$v)
      s.post$inftimes[, sa] <- vars_to_log$inftimes
      s.post$infectors[, sa] <- vars_to_log$infectors
      s.post$nodetimes[, sa] <- vars_to_log$nodetimes[vars_to_log$nodetypes == "c"]
      s.post$nodehosts[, sa] <- vars_to_log$nodehosts[vars_to_log$nodetypes == "c"]
      s.post$nodeparents[, sa] <- vars_to_log$nodeparents
      s.post$introductions[sa] <- sum(vars_to_log$infectors == 0)
      s.post$mu[sa] <- pbe0$p$mu
      #s.post$hist_dens[sa] <- pbe0$h$dist[1]
      #s.post$hist.mean[sa] <- pbe0$p$hist.mean
      s.post$mG[sa] <- pbe0$p$gen.mean
      s.post$mS[sa] <- pbe0$p$sample.mean
      s.post$wh.s[sa] <- pbe0$p$wh.slope
      s.post$wh.e[sa] <- pbe0$p$wh.exponent
      s.post$wh.0[sa] <- pbe0$p$wh.level
      s.post$logLik[sa] <- sum(sapply(names(pbe0)[grepl("logLik", names(pbe0))], function(n) pbe0[[n]]))
      s.post$heat[sa] <- pbe0$heat
      
      if(length(s.post) > 14){
        for(n in setdiff(names(s.post)[15:length(names(s.post))], "chain")){
          if (n == "cnt.types"){
            s.post[[n]][,sa] <- pbe0$p[[n]]
          } else {
            s.post[[n]][sa]<- pbe0$p[[n]]
          }
        }
      }
      # s.post$ir[sa] <- pbe0$p$intro.rate
      # s.post$wh.h[sa] <- pbe0$p$wh.history
      # s.post$dist.e[sa] <- pbe0$p$dist.exponent
      # s.post$dist.s[sa] <- pbe0$p$dist.scale
      # s.post$dist.m[sa] <- pbe0$p$dist.mean
      # s.post$R[sa] <- pbe0$p$R
      
      # if (!is.null(userenv$samplers)){
      #   for (i in 1:length(userenv$samplers)){
      #     s.post[[names(userenv$samplers)[i]]][sa] <- tail(pbe0$p, length(userenv$samplers))[i][[1]] 
      #   }
      # }
      s.posts[[chain]] <- s.post
      
    })
  }
  
  s.posts <- lapply(1:length(heats), function(nheat){
    s <- s.post
    for (chain in 1:nchains){
      smp <- which(s.posts[[chain]]$heat == heats[nheat])
      for (n in names(s)){
        if(inherits(s[[n]], "matrix"))
          s[[n]][,smp] <- s.posts[[chain]][[n]][,smp]
        else
          s[[n]][smp] <- s.posts[[chain]][[n]][smp]
      }
    }
    s$chain <- s.posts[[nheat]]$heat
    return(s)
  })
    
  s.post <- s.posts[[unlist(lapply(1:length(s.posts), function(i){
    if (all(s.posts[[i]]$heat == 1)) return(i)
  }))]]
  
  if(all_chains){
    chains <- lapply(s.posts, function(x){
      x <- list(d = pbe0$d, v = environment2phybreak(pbe0$v), p = pbe0$p, h = pbe0$h, s = x,
                hist = pbe0$v$inftimes[1])
      class(x) <- c("phybreak", "list")
      return(x)
    })
    rm(list = ls(pbe0), envir = pbe0)
    rm(list = ls(pbe1), envir = pbe1)
    return(chains)
  } else {
    res <- destroy_pbe(s.post)
    return(res)
  }
}



#' @rdname sample_phybreak
#' @export
sample.phybreak <- function(...) {
  .Deprecated("sample_phybreak")
  sample_phybreak(...)
}

