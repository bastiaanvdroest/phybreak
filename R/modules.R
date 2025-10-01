#' Create a phybreak-object from data and prior distributions.
#' @export
add_modules_to_phybreak <- function(le, 
                                    multiple.introductions = TRUE, 
                                    spatial = FALSE, 
                                    contact = FALSE, 
                                    infectivity = FALSE, 
                                    phyb.obj = TRUE, ...){
  
  extras <- list(...)
  extras.introductions <- extras[names(extras) %in% names(as.list(args(introductions_parameters)))]
  extras.spatial <- extras[names(extras) %in% names(as.list(args(spatial_parameters)))]
  extras.contact <- extras[names(extras) %in% names(as.list(args(contact_parameters)))]
  extras.infectivity <- extras[names(extras) %in% names(as.list(args(infectivity_parameters)))]
  
  if(phyb.obj == TRUE){
    if(multiple.introductions) le <- do.call(introductions_parameters, c(le, extras.introductions))
    else {
      le$parameterslot$mult.intro = FALSE
      le$parameterslot$wh.history = 1
    }
    
    if(spatial) le <- do.call(spatial_parameters, c(le, extras.spatial))
    else le$parameterslot$spatial = FALSE
    
    if(contact) le <- do.call(contact_parameters, c(le, extras.contact))
    else le$parameterslot$contact = FALSE
    
    if(infectivity | !is.null(extras$infectivity_file) | !is.null(le$dataset$removal.times)) 
      le <- do.call(infectivity_parameters, c(le, extras.infectivity))
    else le$parameterslot$infectivity = FALSE
    
  } else {
    le[["likelihoods"]] <- list()
    le[["updaters"]] <- list()
    if(le$p$mult.intro) le <- introductions_functions(le)
    if(le$p$spatial) le <- spatial_functions(le)
    if(le$p$contact) le <- contact_functions(le)
    if(le$p$infectivity) le <- infectivity_functions(le)
  }
  
  return(le)
}

#####
### Multiple introductions ###

# Add data and parameters to phybreak.object
#' @export
introductions_parameters <- function(le, introductions = 1, 
    wh.history = 1, intro.rate = 1, reproduction.rate = 1,
    est.intro.rate = TRUE, prior.introductions.mean = 1, prior.intro.rate.shape = 0.1,
    est.wh.history = TRUE, prior.wh.history.shape = 1, prior.wh.history.mean = 100,
    use.NJtree = TRUE){
  
  # Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    introductions = introductions,
    intro.rate = intro.rate,
    mult.intro = TRUE,
    use.NJtree = use.NJtree,
    wh.history = wh.history))
  
  # Helperslot
  le$helperslot <- c(le$helperslot, list(
    si.ir = 2.38*sqrt(trigamma(introductions)),
    est.ir = est.intro.rate,
    est.wh.h = est.wh.history,
    ir.sc = prior.introductions.mean / (as.numeric(max(le$dataslot$sample.times) - min(le$dataslot$sample.times)) * prior.intro.rate.shape),
    ir.sh = prior.intro.rate.shape,
    wh.h.sh = prior.wh.history.shape,
    wh.h.av = prior.wh.history.mean))
  
  le$sampleslot <- c(le$sampleslot, list(
    introductions = c(),
    intro.rate = c(),
    wh.history = c()))
  
  return(le)
}

# Add likelihood functions and parameter functions to object
introductions_functions <- function(le){
  le$updaters[["update_wh_history"]] <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.wh.h) return()
    
    ### change to proposal state
    p$wh.history <- exp(log(p$wh.history) + rnorm(1, 0, h$si.wh))
    #if (p$wh.history > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$wh.history) - log(pbe0$p$wh.history)
    
    ### calculate likelihood
    propose_pbe("wh.history")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
      dgamma(pbe1$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE) - 
      dgamma(pbe0$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("wh.history")
    }
  }

  le$updaters[["update_ir"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.ir) return()
    
    p$intro.rate <- exp(log(p$intro.rate) + rnorm(1, 0, h$si.ir))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$intro.rate) - log(pbe0$p$intro.rate)
    
    ### calculate likelihood
    propose_pbe("ir")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen - pbe0$logLikgen + logproposalratio + 
      dgamma(pbe1$p$intro.rate, shape = h$ir.sh, scale = h$ir.sc, log = TRUE) - 
      dgamma(pbe0$p$intro.rate, shape = h$ir.sh, scale = h$ir.sc, log = TRUE)
    
    ### accept
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("ir")
    }
  }

  return(le)
}

#####
### Spatial distance ###
spatial_parameters <- function(le,
    dist.model = "power", dist.exponent = 2, dist.scale = 1, dist.mean = 1,
    est.dist.exponent = TRUE, prior.dist.exponent.shape = 1, prior.dist.exponent.mean = 1,
    est.dist.scale = TRUE, prior.dist.scale.shape = 1, prior.dist.scale.mean = 1,
    est.dist.mean = TRUE, prior.dist.mean.shape = 1, prior.dist.mean.mean = 1){
  
  dist.model <- choose_distmodel(dist.model, le$dataset$distances)
  
  # Dataslot
  le$dataslot <- c(le$dataslot, list(
   locations = le$dataset$locations,
   distances = le$dataset$distances))
  
  if(is.null(le$dataslot$locations)){
    if(!is.null(le$dataslot$distances)) {
      l <- cmdscale(le$dataslot$distances)
    }
  } else {
    l <- le$dataslot$locations
  }
  if(exists("l")) {
    le$dataslot$area <- (max(l[,1])-min(l[,2])) * (max(l[,2]) - min(l[,2]))
  }
  
  #Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    spatial = TRUE,
    dist.model = dist.model,
    dist.exponent = dist.exponent,
    dist.scale = dist.scale,
    dist.mean = dist.mean))
  
  #Helperslot
  le$helperslot <- c(le$helperslot, list(
    si.dist = 2.38*sqrt(trigamma(le$parameterslot$obs - 1)),
    est.dist.e = est.dist.exponent && dist.model %in% c("power", "exponential"),
    est.dist.s = est.dist.scale && dist.model == "power",
    est.dist.m = est.dist.mean && dist.model == "poisson",
    dist.e.sh = prior.dist.exponent.shape,
    dist.e.av = prior.dist.exponent.mean,
    dist.s.sh = prior.dist.scale.shape,
    dist.s.av = prior.dist.scale.mean,
    dist.m.sh = prior.dist.mean.shape,
    dist.m.av = prior.dist.mean.mean))
  
  # Sampleslot
  le$sampleslot <- c(le$sampleslot, list(
    dist.exponent = c(),
    dist.scale = c(),
    dist.mean = c()))
  
  return(le)
}

spatial_functions <- function(le){
  
  # calculate the log-likelihood of distances
  le$likelihoods[["logLikdist"]] <- function(le){#dist.model, dist.exponent, dist.scale, dist.mean, infectors, distances, area) {
    dist.model = le$p$dist.model
    dist.scale = le$p$dist.scale
    dist.exponent = le$p$dist.exponent
    dist.mean = le$p$dist.mean
    infectors = le$v$infectors
    distances = le$d$distances
    area = le$d$area
    
    dist.model <- ifelse(is.null(dist.model), "none", dist.model)
    if(dist.model == "none") return(0)
    distancevector <- distances[cbind(which(infectors!=0), infectors[infectors!=0])]
    sum((infectors == 0) * dunif(1, min = 0, max = area, log = TRUE)) + 
      switch(dist.model,
             power = sum(log(
               dist.exponent * sin(pi/dist.exponent) / 
                 (dist.scale * pi * (1 + (distancevector/dist.scale)^dist.exponent))
             )),
             exponential = sum(
               log(dist.exponent) - dist.exponent * distancevector
             ),
             poisson = sum(
               -dist.mean + distancevector * log(dist.mean) - lgamma(1 + distancevector)
             )
      )
  }
  
  # update parameters of loglik of distances
  le$updaters[["update_dist_exponent"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.exponent <- 1 + exp(log(p$dist.exponent - 1) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.exponent - 1) - log(pbe0$p$dist.exponent - 1)
    
    ### calculate likelihood
    propose_pbe("dist.exponent")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.exponent")
    }
  }
  le$updaters[["update_dist_scale"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.scale <- exp(log(p$dist.scale) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.scale) - log(pbe0$p$dist.scale)
    
    ### calculate likelihood
    propose_pbe("dist.scale")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.scale")
    }
  }
  le$updaters[["update_dist_mean"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.mean <- exp(log(p$dist.mean) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.mean) - log(pbe0$p$dist.mean)
    
    ### calculate likelihood
    propose_pbe("dist.mean")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.mean")
    }
  }
  
  return(le)
}

#####
### Contact parameters ###
#’ Define and initialize contact parameters for a phybreak object
#’
#’ This function extends a \code{phybreak} object with parameters related
#’ to contact structures, including fractions of transmission occurring via different
#’ contact matrices, contact proportions, estimation flags, and priors. It sets
#’ initial values, translates prior means into appropriate parameters, and prepares
#’ empty containers for posterior samples.
#’
#’ @param le A \code{phybreak} object given by the \code{phybreak} function.
#’ @param contact.fracs Initial fractions of transmission attributed to each
#’   contact type. If a single contact matrix is provided, the default is \code{0.5}.
#’   If multiple contact matrices are provided, the default is equal allocation across
#’   matrices (i.e., \code{1 / (number of matrices + 1)}).
#’ @param est.cnt.fracs Logical; whether to estimate the contact fractions (default = \code{TRUE}).
#’ @param prior.cnt.fracs.means Prior mean values for the contact fractions. If not specified,
#’   defaults to a symmetric prior with equal weights (of 1) across all contact types.
#’ @param prior.cnt.fracs.strength Numeric; strength parameter of the prior (interpreted
#’   as a concentration parameter for the Dirichlet prior). Default = 10.
#’ @param contact.prop Initial contact proportion(s). If not specified, calculated as the
#’   mean connectivity of the contact matrix (or matrices).
#’ @param est.cnt.prop Logical; whether to estimate the contact proportion(s) (default = \code{FALSE}).
#’
#’ @return The input object \code{le}, extended with:
#’ \itemize{
#’   \item \strong{parameterslot}: containing contact fractions, contact proportions, and flags.
#’   \item \strong{helperslot}: containing estimation flags and transformed prior values for
#’   contact fractions and contact proportions.
#’   \item \strong{sampleslot}: empty containers for posterior samples of contact-related parameters.
#’ }
#’
#’ @details
#’ This module allows incorporation of heterogeneous contact structures into the model.
#’ A single contact matrix or a list of contact matrices can be provided in the dataset.
#’ Fractions of transmission are initialized either from user input or as uniform defaults,
#’ and priors on these fractions are specified through Dirichlet parameters derived from
#’ \code{prior.cnt.fracs.means} and \code{prior.cnt.fracs.strength}.
#’
#’ @examples
#’ # Example: add contact parameters to a phybreak object
#’ MCMCstate <- phybreak(dataset, contact = T, est.cnt.fracs = TRUE, est.cnt.prop = TRUE)
#’
#’ @seealso \code{\link{phybreak}}, \code{\link{introductions_parameters}}, \code{\link{spatial_parameters}}
#’ @export
contact_parameters <- function(le,
    contact.fracs = NA, est.cnt.fracs = T, 
    prior.cnt.fracs.means = NA, prior.cnt.fracs.strength = 10,
    contact.prop = NA, est.cnt.prop = F){
  
  # Dataslot
  le$dataslot$contact <- le$dataset$contact.matrix
  
  if (inherits(le$dataset$contact.matrix, "matrix")){
    if (is.na(contact.fracs)) contact.fracs <- 0.5
    else if (length(contact.fracs) != 1) 
      stop("Contact fractions vector should be the same length as the number of contact matrices")

    if (any(is.na(contact.prop))){
      contact.prop <- sum(le$dataset$contact.matrix)/(ncol(le$dataset$contact.matrix)^2-ncol(le$dataset$contact.matrix))
    }
  } else if (inherits(le$dataset$contact.matrix, "list")){
    if (is.na(contact.fracs)) contact.fracs <- rep(1/(length(le$dataset$contact.matrix)+1), length(le$dataset$contact.matrix))
    else if (length(contact.fracs) != length(le$dataset$contact.matrix)){
      stop("Contact fractions vector should be the same length as the number of contact matrices")
    }
    if(any(is.na(contact.prop))){ 
      contact.prop <- do.call(c, lapply(le$dataset$contact.matrix, function(m){
        return(sum(m)/(ncol(m)^2-ncol(m)))
      }))
    }
  } else {
    stop("Provide either a contact matrix or a list of contact matrices")
  }    
  
  print(prior.cnt.fracs.strength)
  ### translate means into right format
  if(all(is.na(prior.cnt.fracs.means))) {
    alpha = rep(1,length(contact.fracs)+1)
  } else if (length(prior.cnt.fracs.means) != length(contact.fracs)+1) {
    alpha <- prior.cnt.fracs.strength * c(1-sum(prior.cnt.fracs.means), prior.cnt.fracs.means)
  } else {
    alpha <- prior.cnt.fracs.strength * prior.cnt.fracs.means
  }

  # Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    contact = TRUE,
    contact.fracs = contact.fracs,
    contact.prop = contact.prop))
  
  # Helperslot
  le$helperslot <- c(le$helperslot, list(
    est.cnt.fracs = est.cnt.fracs,
    cnt.fracs.alpha = alpha,
    est.cnt.prop = est.cnt.prop,
    cnt.prop = contact.prop
  ))
  
  # Sampleslot
  le$sampleslot <- c(le$sampleslot, list(
    contact.fracs = NULL,
    contact.prop = NULL))
  return(le)
}

contact_functions <- function(le){
  
  # calculate the log-likelihood of contacts
  le$likelihoods[["logLikcontact"]] <- function(le){

    # Likelihood for contact fractions
    lik.fracs <- with(le, {
      # # Add fraction and proportion of unknown route
      # p$contact.fracs <- p$contact.fracs 
      # p$contact.prop <- c(p$contact.prop, 1)

      R <- p$R

      #   For each host
      if (any(v$infectors > 0)){
        lik.host <- sapply(which(v$infectors != 0), function(i){
          # For each contact route
          lik.i <- sapply(seq_len(dim(contactarray)[3]), function(r){
            # Compute loglikelihood
            return(R*(p$contact.fracs[r]/p$contact.prop[r]) * contactarray[v$infectors[i],i,r])
          })
          return(log(R*(1-sum(p$contact.fracs)) + sum(lik.i)))
        })
      } else {
        lik.host <- 0
      }

      # Remove fraction and proportion of unknown route before storing
      # p$contact.fracs <- p$contact.fracs[-1]
      # p$contact.prop <- p$contact.prop[-1]
      return(sum(lik.host) - R * p$obs)
    })
  }

  le$updaters[["update_contact_fractions"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe0$p

    ### calculate probabilities from rates
    pbe0.theta <- c(p$contact.fracs, 1-sum(p$contact.fracs))

    ### propose new probabilities
    # Function to apply ALR transformation
    # alr_transform <- function(theta) {
    #   return(log(theta[-length(theta)] / theta[length(theta)]))
    # }

    # # Function to apply inverse ALR transformation
    # alr_inverse <- function(z) {
    #   exp_z <- exp(z)
    #   theta <- c(exp_z, 1) / (1 + sum(exp_z))
    #   return(theta)
    # }

    # # Function to propose a new Dirichlet sample
    # propose_dirichlet_logit <- function(old_theta, sd = 0.05) {
    #   # Transform to unconstrained space
    #   z_old <- alr_transform(old_theta)

    #   # Propose in transformed space
    #   z_new <- z_old + rnorm(length(z_old), mean = 0, sd = sd)

    #   # Transform back to simplex
    #   new_theta <- alr_inverse(z_new)

    #   return(new_theta)
    # }
    # pbe1.theta <- propose_dirichlet_logit(pbe0.theta)
    propose_dirichlet_simplex <- function(old_theta, concentration = 50) {
      new_theta <- rgamma(length(old_theta), shape = old_theta * concentration)  
      new_theta <- new_theta / sum(new_theta)  # Normalize to sum to 1  
      return(new_theta)
    } 
    pbe1.theta <- propose_dirichlet_simplex(pbe0.theta)
    
    # pbe1.theta <- pbe0.theta + rnorm(length(pbe0.theta), mean = 0, sd = 0.01)
    # pbe1.theta <- abs(pbe1.theta) / sum(abs(pbe1.theta))


    ### Compute new contact coefficients directly from updated pbe1.theta
    # p$contact.coeff <- sapply(seq_len(length(p$contact.coeff)), function(n) {
    #   (pbe1.theta[n+1] * p$R) / (p$contact.prop[n] * pbe1.theta[1])
    # })

    ### update proposal environment
    p$contact.fracs <- head(pbe1.theta, length(p$contact.fracs))
    copy2pbe1("p", le)
    ### calculate proposalratio
    log_gamma_density <- function(x, shape) {
      return((shape - 1) * log(x) - lgamma(shape))
    }

    log_proposal_ratio <- function(theta_old, theta_new, concentration) {
      log_q_new_given_old <- sum(log_gamma_density(theta_new, shape = concentration * theta_old))
      log_q_old_given_new <- sum(log_gamma_density(theta_old, shape = concentration * theta_new))
      return(log_q_old_given_new - log_q_new_given_old)
    }
    logproposalratio <- log_proposal_ratio(pbe0.theta, pbe1.theta, concentration = 50)
    
    ### calculate likelihood
    propose_pbe("contact")
    
    ### calculate acceptance probability
    #print(c(pbe1$logLikcontact, pbe0$logLikcontact))
    logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio +
    log(dirichlet_pdf(pbe1.theta, h$cnt.fracs.alpha)) - log(dirichlet_pdf(pbe0.theta, h$cnt.fracs.alpha))

    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("contact")
    }
  }

  le$updaters[["update_contact_prop"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe0$p

    ### check whether to update
    if (!h$est.cnt.prop)
      return()

    ### sample 1 of the coefficients
    n <- sample(length(p$contact.prop), 1)

    ### change to proposal state
    prop.new <- rnorm(1, mean = h$cnt.prop, sd = 0.05)
    if (prop.new < 0 || prop.new >= 1)
      return()

    p$contact.prop[n] <- prop.new

    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$contact.prop[n]) - log(pbe0$p$contact.prop[n])
    
    ### calculate likelihood
    propose_pbe("contact.prop")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("contact.prop")
    }
  }

  return(le)
}

#####
### Generation time
infectivity_parameters <- function(le, admission.times = NULL, removal.times = NULL,
    trans.model = "gamma", trans.growth = 1,
    infectivity_file = NULL){
  
  if(!is.null(infectivity_file)) trans.model = "user"
  # Use the Gamma distribution
  if(trans.model == "gamma"){
    # If no removal times, use standard Gamma distribution
    if (is.null(le$dataset$removal.times)) return(le)
    
    # If removal times are present, use adjusted Gamma distribution
    # Dataslot

      if (!is.null(le$dataset$removal.times)) {
        removal.times = le$dataset$removal.times
      }
    
    le$dataslot <- c(le$dataslot, list(
      admission.times = admission.times,
      removal.times = removal.times
    ))
    
    # Parameterslot
    le$parameterslot <- c(le$parameterslot, list(
      infectivity = TRUE,
      removal.rate = 5
    ))
    
    # Helperslot
    le$helperslot <- c(le$helperslot, list(
      est.removal.rate = FALSE 
    ))
    
    # Sampleslot
    le$sampleslot <- c(le$sampleslot, list(
      removal.rate = c()
    ))
    
    le$parameterslot[["inf_function"]] <- function(time, inftimes, le, nodetimes, 
                                               host, log = FALSE,
                                               test.arguments = FALSE) {
  
    p <- le$p
    v <- le$v
    
   # Convert removal times to numeric relative to reference
    all_remtimes <- as.numeric(le$d$removal.times - le$d$reference.date) - v$inftimes
    admission.times <- if (is.null(le$d$admission.times)) {
      rep(0, length(inftimes))
    } else {
      le$d$admission.times[match(inftimes, v$inftimes)] - inftimes
    }

    hosttimes <- as.numeric(time - inftimes)
    remtimes <- as.numeric(le$d$removal.times[match(inftimes, v$inftimes)] - le$d$reference.date) - inftimes

    # --- Normalization factor ---
    # Expected cumulative infectiousness until removal is pgamma(remtime)
    # Compute expected cumulative infectiousness for all hosts
    all_aucs <- pgamma(all_remtimes, shape = p$gen.shape, scale = p$gen.mean / p$gen.shape)

    # Normalization factor using all hosts
    norm_factor <- 1 / mean(all_aucs, na.rm = TRUE)    

    # --- Compute infectiousness ---
    probs <- numeric(length(hosttimes))
    
    before_admission <- hosttimes < admission.times
    before_removal   <- hosttimes >= admission.times & hosttimes <= remtimes
    after_removal    <- hosttimes > remtimes
    
    # gamma density until removal
    probs[before_removal] <- dgamma(
      hosttimes[before_removal],
      shape = p$gen.shape,
      scale = p$gen.mean / p$gen.shape
    )
    
    # decay after removal
    probs[after_removal] <- dgamma(
      remtimes[after_removal],
      shape = p$gen.shape,
      scale = p$gen.mean / p$gen.shape
    ) * exp(-p$removal.rate * (hosttimes[after_removal] - remtimes[after_removal]))
    
    # admission constraint: already 0 by default
    
    # --- Return ---
    if (log) {
      return(log(probs * norm_factor))
    } else {
      return(probs * norm_factor)
    }
  }
    return(le)
  }
  
  else {
    # Check for removal times
    if (is.null(removal.times)){
      if (!is.null(le$dataset$removal.times)) {
        removal.times = le$dataset$removal.times
      } else {
        stop("Provide removal times in same order as hosts")
      }
    }
    
    # Load user-defined infectivity
    # if(trans.model == "user") {
    #   if(is.null(infectivity_file))
    #     stop("Please provide a R file stating the infectivity function")
    #   else 
    #     source(infectivity_file, local = userenv)
    # } else {
    #   datas <- NULL
    #   parameters <- NULL
    #   helpers <- NULL
    #   samplers <- NULL
    # }

    #Add data for user-defined function
    le$dataslot <- c(le$dataslot, list(
      removal.times = removal.times
    ))
    
    # Parameterslot
    le$parameterslot <- c(le$parameterslot, list(
      trans.init = 1e-4,
      trans.removal = 5,
      trans.growth = trans.growth,
      trans.sample = 1,
      trans.model = "user",
      infectivity = TRUE
    ))
    
    # Helperslot
    le$helperslot <- c(le$helperslot, userenv$helperslot)
    
    # Sampleslot
    le$sampleslot <- c(le$sampleslot, userenv$sampleslot)
    
    # Infectivity function
    le$parameterslot[["inf_function"]] <- function(time, inftimes, le, nodetimes, 
                                                                      host = NULL, log = FALSE,
                                                                      test.arguments = FALSE){
      
      d <- le$d
      p <- le$p
      v <- le$v
      # --- argument checks ---
      if (is.null(d$removal.times)) stop("removal times of hosts must be provided")
      if (is.null(p$trans.init))    stop("initial fraction infected is missing")
      if (is.null(p$trans.growth))  stop("growth factor of infectiousness is missing")
      if (is.null(p$trans.sample))  stop("reduction factor after first positive sample is missing")
      if (is.null(p$trans.removal)) stop("decay factor after removal is missing")
      
      if (test.arguments) return()

      # --- make sure all times are numeric ---
      check_time_class <- function(times){
        if (inherits(times, "Date")) {
          if (is.null(d$reference.date)) {
            stop("reference.date must be provided when using Date() input for times")
          }
          return(as.numeric(difftime(times, d$reference.date, units = "days")))
        } else {
          return(as.numeric(times))
        }
      }
      
      removal.times <- check_time_class(d$removal.times)
      inftimes <- check_time_class(inftimes)
      v$inftimes <- check_time_class(v$inftimes)
      time <- check_time_class(time)

      # --- parameters ---
      a <- (1 - p$trans.init) / p$trans.init
      r <- p$trans.growth
      S <- p$trans.sample
      C <- p$trans.removal 

      # --- normalisation: mean AUC per host ---
      st <- as.numeric(v$nodetimes[1:length(v$inftimes)] - v$inftimes)
      rt <- as.numeric(removal.times - v$inftimes)

      AUCs <- ifelse(r * st < 100,
                    (log(a + exp(r * st)) - log(a + 1)) / r +
                    S * (log(a + exp(r * rt)) - log(a + exp(r * st))) / r +
                    (S / (1 + a * exp(-r * rt))) / C,
                    
                    (r * st - log(a + 1)) / r +
                    S * (r * (rt - st)) / r +
                    S / C)
      norm_factor <- 1 / mean(AUCs)

      # --- times relative to infection ---
      remtimes <- removal.times[match(inftimes, v$inftimes)] - inftimes
      samtimes <- as.numeric(nodetimes - inftimes)
      hosttimes <- as.numeric(time - inftimes)

      # --- infectivity calculation ---
      if (length(hosttimes) == 0) {
        probs <- 1
      } else if (is.null(host)) {
        probs <- numeric(length(hosttimes))
        
        # vectorized conditions
        cond1 <- hosttimes < 0
        cond2 <- hosttimes >= 0 & hosttimes < samtimes
        cond3 <- hosttimes >= samtimes & hosttimes < remtimes
        cond4 <- hosttimes >= remtimes & hosttimes < (remtimes + 5)
        
        probs[hosttimes < 0] <- 0
        probs[hosttimes >= 0 & hosttimes < samtimes] <- 1 / (1 + a * exp(-r * hosttimes[cond2]))
        probs[hosttimes >= samtimes & hosttimes < remtimes] <- S / (1 + a * exp(-r * hosttimes[cond3]))
        probs[hosttimes >= remtimes & hosttimes < (remtimes + 5)] <- S / (1 + a * exp(-r * remtimes[cond4])) * 
                        exp(-C * (hosttimes[cond4] - remtimes[cond4]))
        
      } else {
        ht <- hosttimes
        st <- samtimes[host]
        rt <- remtimes[host]
        
        if (ht < 0) {
          probs <- 0
        } else if (ht < st) {
          probs <- 1 / (1 + a * exp(-r * ht))
        } else if (ht < rt) {
          probs <- S / (1 + a * exp(-r * ht))
        } else if (ht < rt + 5) {
          probs <- S / (1 + a * exp(-r * rt)) * exp(-C * (ht - rt))
        } else {
          probs <- 0
        }
      }
      
      out <- probs * norm_factor
      if (log) return(log(out))
      return(out)
    }

      # d <- le$d
      # p <- le$p
      # v <- le$v
      
      # if (is.null(d$removal.times)) {
      #   stop("removal times of hosts must be provided")
      # } else {
      #   removal.times <- d$removal.times
      #   # if(inherits(removal.times, "Date")){
      #   #   removal.times <- as.numeric(removal.times - d$reference.date)
      #   # }
      # }
      
      # if(test.arguments) return()

      # if(is.null(p$trans.init))
      #   stop("initial fraction infected is missing")
      # if(is.null(p$trans.growth))
      #   stop("growth factor of infectiousness is missing")
      # if(is.null(p$trans.sample))
      #   stop("reduction factor after first positive sample is missing")
      # if(is.null(p$trans.removal))
      #   stop("decay factor after removal is missing")
      
      # a <- (1-p$trans.init)/p$trans.init
      # r <- p$trans.growth
      # S <- p$trans.sample
      # C <- p$trans.removal
      
      # # Calculate normalization factor by calculating mean AUC of infectiousness function
      # AUCs <- unlist(lapply(1:length(v$inftimes), function(i){
      #   samtime = as.numeric(v$nodetimes[i] - v$inftimes[i])
      #   cultime = as.numeric(removal.times[i] - v$inftimes[i])
      #   if (r*samtime < 100){
      #     probs = sum((log(a+exp(r*samtime)) - log(a+1)) / r,
      #                 S * ( log(a+exp(r*cultime)) - log(a+exp(r*samtime)) ) / r,
      #                 (S / (1 + a*exp(-r*cultime))) / C)
      #   } else {
      #     probs = sum((r*samtime - log(a+1)) / r,
      #                 S * ( r*(cultime - samtime) ) / r,
      #                 S / C)
      #   }
      #   return(probs)
      # }))
      # norm_factor <- 1/mean(AUCs)
      
      # # Use removal times of infectors in rest of calculations
      # cultimes <- removal.times[match(inftimes, v$inftimes)]
      # samtimes <- as.numeric(nodetimes - inftimes)
      # cultimes <- as.numeric(cultimes - inftimes)
      # hosttimes <- as.numeric(time - inftimes)
      
      # if (length(hosttimes) == 0){
      #   probs = 1
      # } else {
        
      #   if(is.null(host)){
      #     if(length(hosttimes) != length(nodetimes)){
      #       probs <- 0.1
      #       j <- 1
      #     } else {
      #       probs <- c()
      #       j <- 0
      #     }
      #     for (i in 1:length(samtimes)){
      #       if(hosttimes[i+j] < 0)
      #         probs <- c(probs, 0)
      #       else if(hosttimes[i+j] < samtimes[i])
      #         probs <- c(probs, 1/(1+a*exp(-r*hosttimes[i+j])))
      #       else if(hosttimes[i+j] >= samtimes[i] & hosttimes[i+j] < cultimes[i])
      #         probs <- c(probs, S/(1+a*exp(-r*hosttimes[i+j])))
      #       else if(hosttimes[i+j] >= cultimes[i] & hosttimes[i+j] < cultimes[i] + 5)
      #         probs <- c(probs, S/(1+a*exp(-r*cultimes[i])) * exp(-C*(hosttimes[i+j]-cultimes[i])))
      #       else 
      #         probs <- c(probs, 0)
      #     }
      #   } else {
      #     if(hosttimes < 0)
      #         probs <- 0
      #       else if(hosttimes < samtimes[host])
      #         probs <- 1/(1+a*exp(-r*hosttimes))
      #       else if(hosttimes >= samtimes[host] & hosttimes < cultimes)
      #         probs <- S/(1+a*exp(-r*hosttimes))
      #       else if(hosttimes >= cultimes[host] & hosttimes < cultimes + 5)
      #         probs <- S/(1+a*exp(-r*cultimes)) * exp(-C*(hosttimes-cultimes))
      #       else 
      #         probs <- 0
      #   }
      # }
      # if(log)
      #   return(log(probs*norm_factor))
      # else
      #   return(probs*norm_factor)
      # }
    
    return(le)
  }
}

infectivity_functions <- function(le){
  le$updaters[["removal.rate"]] <- function(){
    NULL
  }
  
  return(le)
}
