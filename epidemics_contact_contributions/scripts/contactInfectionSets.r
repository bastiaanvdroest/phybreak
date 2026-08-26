contactInfectionSets <- function(x, which.hosts = "all", samplesize = Inf, 
output = c("host", "matrix", "proportions", "infections")){
  
  obs <- x$p$obs
  hostnames <- unique(x$d$hostnames)
  nhosts <- length(hostnames)

  ### tests
  if (!(is.character(which.hosts) | is.numeric(which.hosts)) | any(is.na(which.hosts))) {
    stop("'which.hosts' should be numeric or \"all\", or should contain exact host names")
  }
  if (!inherits(x$d$contact, 'list')) x$d$contact <- list("contact" = x$d$contact)
  if(is.numeric(which.hosts)) {
    which.hosts <- as.integer(which.hosts)
    which.hosts <- unique(which.hosts)
    which.hosts <- which.hosts[which.hosts >= 1 & which.hosts <= obs]
  } else if (any(which.hosts == "all")) {
    which.hosts <- 1:obs
  } else {
    which.hosts <- match.arg(which.hosts, hostnames, several.ok = TRUE)
    which.hosts <- unique(which.hosts)
    which.hosts <- match(which.hosts, hostnames)
  }

  output <- match.arg(output)

  ### samplesize
  infectors <- x$s$infectors
  chainlength <- ncol(infectors)
  samplesize <- min(samplesize, chainlength)
  samplerange <- (chainlength - samplesize + 1):chainlength

  ### resize data according to given samplesize
  contact.fracs <- x$s$contact.fracs[,samplerange,drop=FALSE]
  r <- x$s$r[samplerange]
  infectors <- infectors[,samplerange]

  ### for each host and each cycle retrieve relative risk
  transmission.routes <- vector("list", nhosts)
  for (cycle in seq_len(ncol(contact.fracs))){
    fracs <- contact.fracs[,cycle,drop=FALSE]
    fracs <- c(1-colSums(fracs), fracs[,1])
    routes <- lapply(seq_len(nrow(infectors)), function(infectee){
      infector <- infectors[infectee,cycle]
      if (infector != 0){
        contact <- c(1, sapply(seq_along(x$d$contact), function(contact.matrix){
          return(ifelse(x$d$contact[[contact.matrix]][infector,infectee] == 1, 1, 0))
        }))
        infect.prob <- c(0, contact * fracs)
        return(infect.prob / sum(infect.prob))
        #trans.route <- sample(1:4, 1, prob = contact * rel_coeff)
        #return(trans.route - 1)
      } else {
        return(c(1,rep(0, length(x$d$contact)+1)))
        #return(-1)
      }
    })

    transmission.routes <- lapply(seq_along(routes), function(host){
      transmission.routes[[host]] <- cbind(transmission.routes[[host]], t(t(routes[[host]])))
    })
  }

  ### Store mean relative risk over complete chain in a list of hosts
  if (output == 'host') {
    res <- lapply(seq_along(transmission.routes), function(i){
      m <- transmission.routes[[i]]
      df <- data.frame(type = c("Introduction", "Unknown", names(x$d$contact)), proportion = rowSums(m) / ncol(m))
      df <- df[order(df$proportion, decreasing = T),]
      return(df)
    })
    names(res) <- unique(x$d$hostnames)
    return(res[which.hosts])
  }

  ### Store relative risk in matrix of nhosts by ntypes of contact, including introduction and unknown
  if (output == 'matrix') {
    contact.matrix <- do.call(cbind,lapply(seq_len(nrow(transmission.routes[[1]])), function(i){
      m <- do.call(rbind, lapply(transmission.routes, function(m) m[i, , drop = FALSE]))
      m <- rowSums(m) / ncol(m)
    }))
    contact.matrix <- t(apply(contact.matrix, 1, function(row){
        row / sum(row)
    }))
    rownames(contact.matrix) <- hostnames
    colnames(contact.matrix) <- c("introduction", "unknown", names(x$d$contact))

    return(contact.matrix[which.hosts, ])
  }
  
  if(output == "proportions"){
    ### Make realization of transmission route for each cycle by sampling from probabilities
    transmission.routes.freq <- do.call(rbind, lapply(transmission.routes, function(m){
      apply(m, 2, function(column){
        if (any(column[2:length(column)] > 0)){
          return(sample(1:(length(column)-1),size = 1, prob = column[2:length(column)]) - 1)
        } else {
          return(-1)
        } 
      })
    }))

    ### Calculate for each cycle the number of transmissions a route was used
    transmission.routes.freq <- t(apply(transmission.routes.freq, 2, function(column){
      return(table(factor(column, levels = -1:length(x$d$contact))))
    }))

    ### Calculate per cycle the proportion of all transmissions that a specific route was used
    transmission.routes.props <- t(apply(transmission.routes.freq, 1, function(row){
      row / sum(row)
    }))
    colnames(transmission.routes.props) <- c("introduction", "unknown", names(x$d$contact))

    return(transmission.routes.props)
  }

  if (output == "infections"){
    ### Count number of infections for each observed route
    infections <- lapply(seq_len(nrow(transmission.routes[[1]])), function(i){
      do.call(rbind, lapply(transmission.routes, function(m){
        r <- m[i, , drop = FALSE]
        r[r > 0] <- 1
        return(r)
      }))
    })
    infections.matrix <- do.call(cbind, lapply(infections, colSums))

    ### Fraction of infections for which route observed
    res <- t(apply(infections.matrix, 1, function(row){
        row / obs
    }))
    colnames(res) <- c("introduction", "unknown", names(x$d$contact))
    return(res)
  }
}
