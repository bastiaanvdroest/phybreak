#' Outbreak simulation.
#' 
#' Simulate outbreaks of class \code{phybreakdata}, with the outbreak model of \pkg{phybreak} (\code{sim.phybreak} is deprecated).
#' 
#' @param obsize The outbreak size (number of cases) to obtain. If \code{obsize = NA}, \code{popsize} should be provided.
#' @param popsize The population size in which to simulate. If it is not defined (default), 
#'   an optimal population size will be chosen based on R0 and obsize. Be aware that choosing a \code{popsize} and
#'   an \code{obsize} can severely increase the simulation time, depending on \code{R0}.
#' @param samplesperhost Number of samples to be taken per host, either a vector or a single number.
#' @param R0 The basic reproduction ratio used for simulation. The offspring distribution is Poisson.
#' @param introductions The number of index cases when simulating a multiple introduced outbreak.
#' @param intro.rate The rate in which introductions appear. The \code{intro.rate} is the parameter of the exponential waiting time distribution between the 
#' introductions.
#' @param spatial If \code{TRUE}, the hosts are placed on a square with density 1, and a distance kernel is used to 
#'   model transmission probabilities between the hosts.
#' @param gen.shape The shape parameter of the gamma-distributed generation interval.
#' @param gen.mean The mean generation interval.
#' @param sample.shape The shape parameter of the gamma-distributed sampling interval.
#' @param sample.mean The mean sampling interval (for the first sample of each host).
#' @param additionalsampledelay Sampling intervals since first sampling times of each host. Values in this vector will be 
#'   used first for all additional samples of host 1, then of host 2, etc.
#' @param wh.model The model for within-host pathogen dynamics (effective pathogen population size = 
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Names and numbers are allowed.
#'   Options are:
#'   \enumerate{
#'     \item "single": effective size = 0, so coalescence occurs 'just before' transmission in the infector (complete bottleneck)
#'     \item "infinite": effective size = Inf, with complete bottleneck, so coalescence occurs 'just after' transmission in the infectee
#'     \item "linear": effective size at time t after infection = \code{wh.level + wh.slope * t} (complete or wide bottleneck; if complete, \code{wh.level = 0})
#'     \item "exponential": effective size at time t after infection = \code{wh.level * exp(wh.exponent * t)} (wide bottleneck)
#'     \item "constant": effective size = wh.level (wide bottleneck)
#'   }
#' @param wh.bottleneck Whether the bottleneck should be complete or wide, which is only an option if \code{wh.model = "linear"} 
#'   (in that case, \code{"auto"} defaults to \code{"complete"}).
#' @param wh.slope Within-host slope, used if \code{wh.model = "linear"}.
#' @param wh.exponent Within-host exponent, used if \code{wh.model = "exponential"}
#' @param wh.level Within-host effective pathogen size at transmission, used if \code{wh.bottleneck = "wide"}
#'   (if \code{wh.model = "exponential"} or \code{"constant"}, and optional if \code{wh.model = "linear"})
#' @param wh.history Within-host coalescent rate of the history host. Rate = 1/\code{wh.history}.
#' @param dist.model The distance kernel to use if \code{spatial = TRUE}. Options are:
#'   \enumerate{
#'     \item "power": a power law function pr(dist) ~ 1 / (1 + (dist/dist.scale) ^ dist.exponent)
#'     \item "exponential": an exponential function pr(dist) ~ exp(-dist.exponent * dist)
#'   }
#' @param dist.exponent Distance model exponent.
#' @param dist.scale Distance model scale, only with power law distance model.
#' @param mu Expected number of mutations per nucleotide per unit of time along each lineage. 
#' @param sequence.length Number of available nucleotides for mutations.
#' @param ... If arguments from previous versions of this function are used, they may be interpreted correctly through 
#'   this argument, but it is better to provide the correct argument names.
#' @return The simulation output as an object of class \code{'phybreakdata'} with sequences (class \code{'phyDat'}) and 
#'   sampling times (which would be the observations), and infection times, infectors, and phylogenetic tree 
#'   of class \code{\link[ape]{phylo}}.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' simulation <- sim_phybreak()
#' @export
sim_phybreak <- function(obsize = 50, popsize = NA, samplesperhost = 1,
                         introductions = 1, intro.rate = 1, outbreak.duration = 10,
                         R0 = 1.5, spatial = FALSE, contact = FALSE, contact.symmetric = TRUE, contact.probs = 1,
                         gen.shape = 10, gen.mean = 1, 
                         sample.shape = 10, sample.mean = 1,
                         additionalsampledelay = 0,
                         wh.model = "linear", wh.bottleneck = "auto", wh.slope = 1, wh.exponent = 1, wh.level = 0.1,
                         wh.history = 100,
                         dist.model = "power", dist.exponent = 2, dist.scale = 1,
                         cnt.invest.trans = 1, cnt.invest.nontrans = 0.1, cnt.rep = 0.9, cnt.rep.false = 0,
                         mu = 0.0001, sequence.length = 10000, ...) {
  
  ### parameter name compatibility 
  old_arguments <- list(...)
  if(exists("old_arguments$shape.gen")) gen.shape <- old_arguments$shape.gen
  if(exists("old_arguments$mean.gen")) gen.mean <- old_arguments$mean.gen
  if(exists("old_arguments$shape.sample")) sample.shape <- old_arguments$shape.sample
  if(exists("old_arguments$mean.sample")) sample.mean <- old_arguments$sample.mean

 
  ### tests
  if(all(is.na(c(obsize, popsize)))) {
    stop("give an outbreak size (obsize) and/or a population size (popsize)")
  }
  if(all(!is.na(c(obsize, popsize)))) {
    warning("giving both an outbreak size (obsize) and a population size (popsize) can take a very large simulation time",
            immediate. = TRUE)
  }
  if(all(!is.na(c(obsize, popsize))) && obsize > popsize) {
    stop("outbreak size (obsize) cannot be larger than population size (popsize)")
  }
  if(R0 <= 1) {
    stop("R0 should be larger than 1")
  }
  if(any(c(gen.shape, gen.mean, sample.shape, sample.mean, wh.slope, mu) <= 0)) {
    stop("parameter values should be positive")
  }
  wh.model <- choose_whmodel(wh.model)
  wh.bottleneck <- choose_whbottleneck(wh.bottleneck, wh.model)
  wh.level <- wh.level * (wh.bottleneck == "wide")
  
  ### simulate step by step
  if(is.na(obsize)) {
    if(spatial) {
      res <- sim_outbreak_spatial(popsize, R0, introductions, intro.rate, gen.shape, gen.mean,
                          sample.shape, sample.mean, dist.model, dist.exponent, dist.scale)
    } else {
      res <- sim_outbreak(popsize, R0, introductions, intro.rate, gen.shape, gen.mean,
                          sample.shape, sample.mean)
    }
     obsize <- res$obs
     if(obsize == 1) return(c("Outbreak size = 1"))
  } else {
    if(spatial) {
      res <- sim_outbreak_size_spatial(obsize, popsize, R0, introductions, intro.rate, gen.shape, gen.mean,
                               sample.shape, sample.mean, dist.model, dist.exponent, dist.scale)
    } else {
      res <- sim_outbreak_size(obsize, popsize, introductions, intro.rate, outbreak.duration, R0, gen.shape, gen.mean,
                               sample.shape, sample.mean)
    }
  }
  if(any(samplesperhost < 1)) stop("samplesperhost should be positive")
  if(any(additionalsampledelay < 0)) stop("additionalsampledelay cannot be negative")
  res <- sim_additionalsamples(res, samplesperhost, additionalsampledelay)
  res <- sim_phylotree(res, wh.model, wh.bottleneck, wh.slope, wh.exponent, wh.level, wh.history, sample.mean)
  res <- sim_sequences(res, mu, sequence.length)
  if(contact)
    res <- sim_contact_matrix(res, contact.symmetric, contact.probs, 
                              cnt.invest.trans, cnt.invest.nontrans, cnt.rep, cnt.rep.false)

  hostnames <- paste0("host.", 1:obsize)
  samplenames <- paste0("sample.", res$nodehosts[1:res$Nsamples], ".", nthsample(res))

  names(res$sequences) <- samplenames
  if(spatial) {
    rownames(res$locations) <- hostnames
  }
  if(contact) {
    rownames(res$contact.matrix) <- colnames(res$contact.matrix) <- hostnames
  }
  
  ### make a phylo tree
  treeout <- phybreak2phylo(vars = environment2phybreak(res), samplenames = samplenames, simmap = FALSE)
  if(spatial) {
    toreturn <- with(res,
                     phybreakdata(
                       sequences = sequences,
                       sample.times = c(samtimes, addsampletimes),
                       spatial = locations,
                       sample.names = samplenames,
                       host.names = hostnames[c(1:obs, addsamplehosts)],
                       sim.infection.times = inftimes,
                       sim.infectors = infectors,
                       sim.tree = treeout
                ))
  } else {
    toreturn <- with(res,
                     phybreakdata(
      sequences = sequences,
      sample.times = c(samtimes, addsampletimes),
      sample.names = samplenames,
      host.names = hostnames[c(1:obs, addsamplehosts)],
      sim.infection.times = inftimes,
      sim.infectors = infectors,
      sim.tree = treeout
    ))
    #toreturn$sim.hist.time <- histtime
  }
  toreturn$admission.times <- res$admissiontimes
  if (contact){
    toreturn$contact.matrix <- res$contact.matrix
    toreturn$contact.categories <- res$contact.categories
  }
  return(toreturn)
}

#' @rdname sim_phybreak
#' @export
sim.phybreak <- function(...) {
  .Deprecated("sim_phybreak")
  sim_phybreak(...)
}


### simulate an outbreak of a particular size by repeating
### simulations until obsize is obtained
sim_outbreak_size <- function(obsize, Npop, intronr, introrate, duration, R0, aG, mG, aS, mS) {
  if(is.na(Npop)) {
    Npop <- obsize
    while(1 - obsize/Npop < exp(-R0* obsize/Npop)) {Npop <- Npop + 1}
  } 
  
  sim <- sim_outbreak(Npop, intronr, introrate, duration, R0, aG, mG, aS, mS)
  
  while(sim$obs != obsize) {
    sim <- sim_outbreak(Npop, intronr, introrate, duration, R0, aG, mG, aS, mS)
  }
  
  return(sim)
}

sim_outbreak_size_spatial <- function(obsize, Npop, R0, intronr, introrate, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale) {
  if(is.na(Npop)) {
    Npop <- obsize
    while(1 - obsize/Npop < exp(-R0* obsize/Npop)) {Npop <- Npop + 1}
  } 
  
  sim <- sim_outbreak_spatial(Npop, R0, intronr, introrate, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale)
  
  while(sim$obs != obsize) {
    sim <- sim_outbreak_spatial(Npop, R0, intronr, introrate, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale)
  }
  
  return(sim)
}


### simulate an outbreak
sim_outbreak <- function(Npop, intronr, introrate, duration, R0, aG, mG, aS, mS) {
  ### initialize
  introintervals <- rexp(intronr - 1, rate = introrate)
  inftimes <- c(0, cumsum(introintervals), rep(10000, Npop - intronr))
  sources <- rep(0, Npop)
  nrcontacts <- rpois(Npop, R0)
  admission.dates <- sort(round(runif(Npop, max = 10)))
  nth.infection <- 1
  currentID <- 1
  
  ### by order of infection, sample secondary infections
  # currentID is the infected host under consideration
  while(nth.infection <= Npop & inftimes[currentID] != 10000) {
    #when does currentID make infectious contacts?
    #reverse sorting so that with double contacts, the earliest will be used last
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG), decreasing = TRUE)
    
    #who are these contacts made with?
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    
    #are these contact possible, i.e. was the host present at that time?
    #possible <- admission.dates[whocontacted] < whencontacts
    possible <- whencontacts == whencontacts
    
    #are these contacts successful, i.e. earlier than the existing contacts with these hosts?
    successful <- whencontacts < inftimes[whocontacted] & whencontacts < duration & whocontacted > intronr
    
    #change infectors and infection times of successful contactees
    sources[whocontacted[possible & successful]] <- currentID
    inftimes[whocontacted[possible & successful]] <- whencontacts[possible & successful]
    
    #go to next infected host in line
    nth.infection <- nth.infection + 1
    currentID <- order(inftimes)[nth.infection]
  }
  
  ### determine outbreaksize and sampling times
  obs <- sum(inftimes<10000)
  samtimes <- inftimes + rgamma(Npop, aS, aS/mS)
  
  ### order hosts by sampling times, and renumber hostIDs
  ### so that the uninfected hosts can be discarded
  orderbysamtimes <- order(samtimes)
  sources <- sources[orderbysamtimes]
  infectors <- match(sources,orderbysamtimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamtimes][1:obs]
  samtimes <- samtimes[orderbysamtimes][1:obs]
  adtimes <- admission.dates[orderbysamtimes[1:obs]]
  
  
  return(
    list(
      obs = obs,
      samtimes = samtimes,
      inftimes = inftimes,
      infectors = infectors,
      admissiontimes = adtimes
    )
  )
}

### simulate a spatial outbreak
sim_outbreak_spatial <- function(Npop, R0, intronr, introrate, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale) {
  ### initialize spatial population model
  x <- runif(Npop, 0, sqrt(Npop))
  y <- runif(Npop, 0, sqrt(Npop))
  distances <- as.matrix(dist(cbind(x, y)))
  dist_densities <- if(dist.model == "exponential") {
    dist.exponent * exp(-dist.exponent * distances)
  } else {
    dist.exponent * sin(pi/dist.exponent) / (pi * dist.scale * (1 + (distances/dist.scale) ^ dist.exponent))
  }
  dist_densities[cbind(1:Npop, 1:Npop)] <- 0
  # matrixR0 <- max(eigen(dist_densities)$values)
  matrixR0 <- sum(dist_densities)/Npop
  R0_matrix <- R0 * dist_densities / matrixR0

  ### initialize outbreak
  introintervals <- rexp(intronr - 1, rate = introrate)
  inftimes <- c(0, cumsum(introintervals), rep(10000, Npop - intronr))
  sources <- rep(0, Npop)
  nrcontacts <- rpois(Npop, rowSums(R0_matrix))
  nth.infection <- 1
  currentID <- 1
  
  ### by order of infection, sample secondary infections
  # currentID is the infected host under consideration
  while(nth.infection <= Npop & inftimes[currentID] != 10000) {
    #when does currentID make infectious contacts?
    #reverse sorting so that with double contacts, the earliest will be used last
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    
    #who are these contacts made with?
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE, prob = R0_matrix[currentID, ])
    
    #are these contacts successful, i.e. earlier than the existing contacts with these hosts?
    successful <- whencontacts < inftimes[whocontacted] & whocontacted > intronr
    
    #change infectors and infection times of successful contactees
    sources[whocontacted[successful]] <- currentID
    inftimes[whocontacted[successful]] <- whencontacts[successful]
    
    #go to next infected host in line
    nth.infection <- nth.infection + 1
    currentID <- order(inftimes)[nth.infection]
  }
  
  ### determine outbreaksize and sampling times
  obs <- sum(inftimes<10000)
  samtimes <- inftimes + rgamma(Npop, aS, aS/mS)
  
  ### order hosts by sampling times, and renumber hostIDs
  ### so that the uninfected hosts can be discarded
  orderbysamtimes <- order(samtimes)
  sources <- sources[orderbysamtimes]
  infectors <- match(sources,orderbysamtimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamtimes]
  samtimes <- samtimes[orderbysamtimes]
  locations <- cbind(x, y)[orderbysamtimes, ]

  ### return the outbreak
  return(
    list(
      obs = obs,
      samtimes = samtimes[1:obs],
      inftimes = inftimes[1:obs],
      infectors = infectors,
      locations = locations[1:obs, ]
    )
  )
}


### simulate additional samples in a transmission tree
sim_additionalsamples <- function(sim.object, samperh, addsamdelay) {
  # recycle too short arguments
  addsamplesizes <- rep_len(samperh - 1, sim.object$obs)
  addsamplesizes[addsamplesizes < 0] <- 0
  alldelays <- rep_len(addsamdelay, sum(addsamplesizes))
  
  # vectors with additional samplehosts and sample times
  addsamplehosts <- rep(1:sim.object$obs, addsamplesizes)
  addsampletimes <- sim.object$samtimes[addsamplehosts] + alldelays
  addsampletimes <- addsampletimes[order(addsamplehosts, addsampletimes)]
  
  return(within(sim.object, {
    Nsamples <- sim.object$obs + sum(addsamplesizes)
    addsamplehosts <- addsamplehosts
    addsampletimes <- addsampletimes
  }))
}

### simulate a phylogenetic tree given a transmission tree
sim_phylotree <- function (sim.object, wh.model, wh.bottleneck, wh.slope, wh.exponent, wh.level, wh.history, sample.mean) {
  list2env(list(v = sim.object, 
                p = list(wh.model = wh.model, wh.bottleneck = wh.bottleneck, wh.slope = wh.slope, wh.exponent = wh.exponent,
                         wh.level = wh.level, wh.history = wh.history, sample.mean = sample.mean),
                d = list(nsamples = sim.object$Nsamples)), pbe1)
  
  pbe1$v$nodeparents <- rep(-1, 2 * sim.object$Nsamples + sim.object$obs - 1)  #initialize nodes: will contain parent node in phylotree
  pbe1$v$nodetimes <- c(sim.object$samtimes, sim.object$addsampletimes, 
                        rep(0, sim.object$Nsamples - 1), sim.object$inftimes)   #initialize nodes: will contain time of node
  pbe1$v$nodehosts <- c(1:sim.object$obs, sim.object$addsamplehosts, 
                        rep(-1, sim.object$Nsamples - 1), sim.object$infectors)   #initialize nodes: will contain host carrying the node
  pbe1$v$nodetypes <- c(rep("s", sim.object$obs), rep("x", sim.object$Nsamples - sim.object$obs), 
                        rep("c", sim.object$Nsamples - 1), rep("t", sim.object$obs))  #initialize nodes: will contain node type (sampling, additional sampling, coalescent)
  
  
  if(wh.bottleneck == "wide") {
    invisible(sapply(1:sim.object$obs, rewire_pullnodes_wide))
  } else {
    invisible(sapply(0:sim.object$obs, rewire_pullnodes_complete))
  }  
  res <- as.list.environment(pbe1)$v
  return(res)
  
}

### simulate sequences given a phylogenetic tree
sim_sequences <- function (sim.object, mu, sequence.length) {
  with(sim.object,{
    ### simulate the mutations on the phylotree
    #number of mutations
    edgelengths <- nodetimes - c(0, nodetimes)[1 + nodeparents]
    edgelengths[edgelengths < 0] <- 0  #rounding errors close to 0
    nmutations <- rpois(1, mu * sequence.length * sum(edgelengths))
    #place mutations on edges, order by time of edge (end)
    mutedges <- sample(length(edgelengths), size = nmutations, replace = TRUE, prob = edgelengths)
    mutedges <- mutedges[order(nodetimes[mutedges])]
    #sample mutations: which locus, to which nucleotide
    mutsites <- sample(sequence.length, size = nmutations, replace = TRUE)
    mutsites <- match(mutsites, unique(mutsites))  #place mutations at front of sequence 
    mutnucl <- sample(4, size = nmutations, replace = TRUE)
    
    ### construct the strains from the simulation by going backwards
    ### through the phylotree from each tip and placing mutations
    nodestrains <- matrix(data = rep(sample(4, nmutations, replace = TRUE), each = Nsamples), nrow = Nsamples)
    for(i in 1:Nsamples) {
      currentedge <- i
      recentmutations <- rep(FALSE, nmutations) #keep more recent mutations on each locus
      while(nodeparents[currentedge] != 0) {
        nodestrains[i, mutsites[mutedges == currentedge & !recentmutations]] <-
          mutnucl[mutedges == currentedge & !recentmutations]
        recentmutations <- recentmutations | mutedges == currentedge
        currentedge <- nodeparents[currentedge]
      }
    }
    # place single unchanged acgt at front of sequence to force these at first positions in phyDat-object
    nodestrains <- cbind(matrix(data = rep(1:4, each = Nsamples), ncol = 4), nodestrains)
    
    nodestrains[nodestrains == 1] <- "a"
    nodestrains[nodestrains == 2] <- "c"
    nodestrains[nodestrains == 3] <- "g"
    nodestrains[nodestrains == 4] <- "t"
    
    rownames(nodestrains) <- 1:Nsamples
    
    # make phyDat-object and change attributes to get entire sequence, with single acgt at front removed
    nodestrains <- phangorn::as.phyDat(nodestrains)
    mutlocs <- sample(4, max(0, sequence.length - nmutations), replace = TRUE)
    attr(nodestrains, "index") <- c(attr(nodestrains, "index")[-(1:4)], mutlocs)
    attr(nodestrains, "weight")[1:4] <- attr(nodestrains, "weight")[1:4] + tabulate(mutlocs, 4) - 1
    
    return(
      within(sim.object,{
        sequences <- nodestrains
      })
    )
  }
  )
}

sim_contact_matrix <- function(sim.object, contact.symmetric, contact.probs, 
                               cnt.invest.trans, cnt.invest.nontrans, cnt.rep, cnt.rep.false) {
  with(sim.object, {
    n = obs
    
    if (!is.matrix(contact.probs))
      if (contact.probs == 1) 
        contact.probs <- t(contact.probs)
    
    categories <- sample(1:nrow(contact.probs), n, replace = TRUE)
    if (contact.symmetric){
      if (!isSymmetric(contact.probs))
        stop("Contact probabilities should be a symmetric matrix if contact.symmetric = TRUE")
    }
    
    # transmission pair with contact
    tp.wc <- (1-cnt.invest.trans)*cnt.rep.false + cnt.invest.trans*cnt.rep
    # non-transmission pair with contact
    np.wc <- (1-cnt.invest.nontrans)*cnt.rep.false + cnt.invest.nontrans*cnt.rep
    # transmission pair no contact
    tp.nc <- (1-cnt.invest.trans)*(1-cnt.rep.false) + cnt.invest.trans*(1-cnt.rep)
    # non-transmission pair no contact
    np.nc <- (1-cnt.invest.nontrans)*(1-cnt.rep.false) + cnt.invest.nontrans*(1-cnt.rep)
    
    mat <- matrix(0, nrow = n, ncol = n)
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        if (contact.symmetric){
          if (infectors[i] == j | infectors[j] == i){
            mat[i,j] <- sample(c(0,1), 1, prob = c(tp.nc, tp.wc) * contact.probs[categories[i],categories[j]])
          } else {
            mat[i,j] <- sample(c(0,1), 1, prob = c(np.nc,  np.wc) * contact.probs[categories[i],categories[j]])
          }
        } else {
          if (infectors[i] == j) mat[i,j] <- sample(c(0,1), 1, prob = c(tp.nc, tp.wc) * contact.probs[categories[i],categories[j]])
          else if (infectors[j] == i) mat[j,i] <- sample(c(0,1), 1, prob = c(tp.nc, tp.wc) * contact.probs[categories[i],categories[j]])
          else {
            mat[i,j] <- sample(c(0,1), 1, prob = c(np.nc, np.wc) * contact.probs[categories[i],categories[j]])
            mat[j,i] <- sample(c(0,1), 1, prob = c(np.nc, np.wc) * contact.probs[categories[i],categories[j]])
          }
        }
      }
    }
    if (contact.symmetric) contact <- mat + t(mat)
    else contact <- mat
    
    return(within(sim.object,{
      contact.matrix <- contact
      contact.categories <- categories
    }))
  })
}

nthsample <- function(sim.object) {
  with(sim.object, {
    nth <- rep(0, Nsamples)
    sapply(1:obs, function(x) suppressWarnings(nth[which(nodehosts[nodetypes %in% c("s", "x")] == x)] <<- 0:Nsamples))
    return(nth)
  })
}

get_contact_matrix <- function(n, tp.wc, np.wc, tp.nc, np.nc){
  mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (sim$sim.infectors[i] == paste0("host.",j) | sim$sim.infectors[j] == paste0("host.",i)){
        mat[i,j] <- sample(c(0,1), 1, prob = c(tp.nc, tp.wc))
      } else {
        mat[i,j] <- sample(c(0,1), 1, prob = c(np.nc,  np.wc))
      }
    }
  }
  
  #mat[upper.tri(mat)] <- sample(c(0,1), sum(upper.tri(mat)), replace = T)
  mat <- mat + t(mat)
  return(mat)
}

contact_matrix_probability <- function(m, infectors, tp.wc, np.wc, tp.nc, np.nc){
  lik <- c()
  for(i in 1:ncol(m)){
    for(j in 1:ncol(m)){
      if (i != j){
        if (infectors[j] == i | infectors[i] == j){
          if (m[i,j] == 1){
            lik <- c(lik, tp.wc)
          } else{
            lik <- c(lik, tp.nc)
          }
        } else {
          if (m[i,j] == 1){
            lik <- c(lik, np.wc)
          } else {
            lik <- c(lik, np.nc)
          }
        }
      }
    }
  }
  return(prod(lik))
}

