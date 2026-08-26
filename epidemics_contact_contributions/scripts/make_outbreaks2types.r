devtools::load_all("~/Documents/phybreak/R")

select_random_matrix_entries <- function(mat, prop = 0.1) {
    stopifnot(is.matrix(mat), length(dim(mat)) == 2)
    ones <- which(mat == 1)
    if (length(ones) == 0) return(mat)
    k <- max(1, round(length(ones) * prop))
    sel <- sample(ones, k)
    mat[sel] <- 0
    return(mat)
}

select_random_10pct_20x20 <- function(mat) {
    stopifnot(is.matrix(mat), all(dim(mat) == c(20, 20)))
    select_random_matrix_entries(mat, 0.1)
}

outbreak.list <- lapply(1:10, function(i){
    sim <- sim_phybreak(obsize = 20)
    print(sim$sample.times)
})
outbreak.list <- list.files("~/Documents/phd_files/contact/results/outbreaks/blank", full.names = TRUE) %>%
  gsub("\\D", "", .) %>%               # Remove non-numeric parts
  as.numeric() %>%                     # Convert to numeric
  order() %>%                          # Get sorted indices
  (\(x) list.files("~/Documents/phd_files/contact/results/outbreaks/blank", full.names = TRUE)[x])() %>%
  lapply(readRDS)

## Independent contact routes
#df <- data.frame(V1 = rep(c("0.1-0.1", "0.5-0.5", "1-1"), each = 3), V2 = rep(c("5-3", "15-4", "19-0"), 3))
df <- data.frame(V1 = rep(c("0.1-0.1", "0.5-0.5", "1-1"), each = 1), V2 = rep(c("4-0"), 3))

for (nr in seq_along(outbreak.list)){
    sim <- outbreak.list[[nr]]
    print(nr)
    for (i in seq_len(nrow(df))){
        infectors <- sapply(sim$sim.infectors, function(x){
            if (x == "index") return(0)
            else return(as.numeric(substr(x, 6, nchar(x))))
        })
        obs=20
        
        prop <- as.numeric(unlist(strsplit(df[i,1], "-")))
        probtrans <- as.numeric(unlist(strsplit(df[i,2], "-")))
        print(prop)
        p <- sim_contact_matrix(list(sim), probtrans, prop)
        toreturn <- with(p[[1]], phybreakdata(
            sequences = sequences,
            sample.times = sample.times,
            sample.names = names(sample.times),
            host.names = sample.hosts,
            sim.infection.times = sim.infection.times ,
            sim.infectors = sim.infectors ,
            sim.tree = sim.tree
        ))
        toreturn$contact.matrix <- p$contact.matrix
        saveRDS(toreturn, 
                sprintf("~/Documents/phd_files/contact/results/outbreaks_numbers/2types/contact_obsize20_intro1_2type_prop%s_probtrans%s_%s.RDS", 
                        paste(prop, collapse="-"), paste(probtrans, collapse="-"), nr))
    }
}

## Correlated infection routes
set.seed(123)
createCorrelationMatrix <- function(sim, nr, proportions, probtrans){    
    #sim <- outbreak.list[[nr]]
    infectors <- sapply(sim$sim.infectors, function(x){
        if (x == "index") return(0)
        else return(as.numeric(substr(x, 6, nchar(x))))
    })
    obs=20
    n = sum(infectors > 0)

    p1 <- proportions[1]
    p2 <- proportions[2]
    p12 <- proportions[3]
    
    prop <- c(1-(p1+p2-p12), p1-p12, p2-p12, p12)
    print(prop)
    prop.m <- c(prop[2]+prop[4], prop[3]+prop[4])

    # Sample route of transmission pairs with probability above
    random_vector <- rep(0, n)
    indices <- sample(1:n, sum(n))
    start <- 1
    category <- 1
    while (category <= length(probtrans)) {
      if (probtrans[category] == 0) {
        category <- category + 1
        next
      }
      end <- start + probtrans[category] - 1
      random_vector[indices[start:end]] <- category
      start <- end + 1
      category <- category + 1
      next
    }
    tp_contacts <- random_vector
    # Split random_vector into n vectors
    # tp_contacts <- lapply(seq_len(length(probtrans)), function(category) {
    #   x <- ifelse(random_vector == category, category, 0)
    #   ifelse(x > 0, 1, 0)
    # })

    # Create contact data matrix by fill in transmission pair contacts and sample nontransmission pair contact
    matrix.list <- lapply(seq_along(probtrans), function(i) {
        m <- matrix(NA, nrow = obs, ncol = obs)
        diag(m) <- 0
        return(m)
    })
    
    row=1
    while(row < obs){
        for (column in (row+1):obs){
            if ((infectors[row] == column || infectors[column] == row) & tp_contacts[row]>0){
                r = tp_contacts[row]
                rn = setdiff(seq_along(matrix.list), r)
                matrix.list[[r]][row,column] = 1
                matrix.list[[rn]][row,column] = sample(c(0,1), 1, prob = c(prop[4-r]/prop.m[r], prop[4]/prop.m[r]))
            } else {
                p <- sample(1:4, 1, prob = prop)
                if (p == 1){
                    matrix.list[[1]][row, column] = 0
                    matrix.list[[2]][row, column] = 0
                } else if (p == 2){
                    matrix.list[[1]][row, column] = 0
                    matrix.list[[2]][row, column] = 1
                } else if (p == 3){
                    matrix.list[[1]][row, column] = 1
                    matrix.list[[2]][row, column] = 0
                } else if (p == 4){
                    matrix.list[[1]][row, column] = 1
                    matrix.list[[2]][row, column] = 1
                }
            }
        }
        row <- row+1
    }

    for (r in seq_along(matrix.list)){
        m <- matrix.list[[r]]
        m[lower.tri(m)] <- t(m)[lower.tri(m)]
        matrix.list[[r]] <- m
    }

    toreturn <- with(sim, phybreakdata(
        sequences = sequences,
        sample.times = sample.times,
        sample.names = names(sample.times),
        host.names = sample.hosts,
        sim.infection.times = sim.infection.times ,
        sim.infectors = sim.infectors ,
        sim.tree = sim.tree
    ))
    toreturn$contact.matrix <- matrix.list
    saveRDS(toreturn, 
                sprintf("~/Documents/phd_files/contact/results/outbreaks_numbers/2types_corr/contact_obsize20_intro1_2type_prop%s_probtrans%s_%s.RDS", 
                        paste(prop, collapse="-"), paste(probtrans, collapse="-"), nr))
    # return(toreturn)
}
lapply(seq_along(outbreak.list), function(i) {
    print(i)
    createCorrelationMatrix(outbreak.list[[i]], i, proportions = c(0.1,0.1,0.1), probtrans = c(4,0))
})

## Block matrices
set.seed(123)

createBlockMatrix <- function(sim, nr){ 
    infectors <- sapply(sim$sim.infectors, function(x){
        if (x == "index") return(0)
        else return(as.numeric(substr(x, 6, nchar(x))))
    })

    df <- data.frame(from = infectors, to = 1:20)
    # Create an adjacency list
    tree <- split(df$to, df$from)

    find_paths <- function(node, path = c(), paths = list()) {
        path <- c(path, node)
        
        # If the node has no children, it's a leaf
        if (!as.character(node) %in% names(tree)) {
            paths[[length(paths) + 1]] <- path
        } else {
            for (child in tree[[as.character(node)]]) {
            paths <- find_paths(child, path, paths)
            }
        }
    
        return(paths)

    }
    paths <- find_paths(0)

    obs=20
        
    prop <- c(0.1, 0.1)
    prop1 <- c(2/3, 1/6, 1/6)
    prop2 <- c(3/4, 1/4)
    probtrans <- c(0.7, 0.2)

    tp_contact_prob <- sapply(seq_along(prop), function(i){
        prop[i]*probtrans[i]/sum(c(1,prop)*c(0.1,probtrans))
    })
    
    # Assign letters to nodes with probabilities 1/3, 1/3, 1/3
    node_labels <- list(list(), list())

    for (path in paths) {
        for (node in path) {
            if (!node %in% names(node_labels[[1]])) {
                if (node != 0) {
                    route = sample(c(0,1,2), 1, prob = c(1-sum(tp_contact_prob), tp_contact_prob[1], tp_contact_prob[2]))
                    if (route == 0 || path[which(path == node)-1] == 0) {
                        node_labels[[1]][[as.character(node)]] = sample(1:3, 1, prob = prop1)
                        node_labels[[2]][[as.character(node)]] = sample(1:2, 1, prob = prop2)
                    } else if (path[which(path == node)-1] != 0){
                        node_labels[[route]][as.character(node)] = node_labels[[route]][as.character(path[which(path==node)-1])]
                        if (route == 1) node_labels[[2]][[as.character(node)]] = sample(1:2, 1, prob = prop2)
                        if (route == 2) node_labels[[1]][[as.character(node)]] = sample(1:3, 1, prob = prop1)
                    }
                }
            }
        }
    }
    node_labels <- lapply(node_labels, function(lab){
        return(unlist(lab)[order(as.numeric(names(lab)))])
    })
    # Create matrices for node labels
    label_matrices <- lapply(node_labels, function(labels) {
        n <- length(labels)
        mat <- matrix(0, nrow = n, ncol = n)
        rownames(mat) <- names(labels)
        colnames(mat) <- names(labels)
        
        for (i in seq_along(labels)) {
            for (j in seq_along(labels)) {
                if (labels[i] == labels[j]) {
                    mat[i, j] <- 1
                }
            }
        }
        return(mat)
    })

    toreturn <- with(sim, phybreakdata(
        sequences = sequences,
        sample.times = sample.times,
        sample.names = names(sample.times),
        host.names = sample.hosts,
        sim.infection.times = sim.infection.times ,
        sim.infectors = sim.infectors ,
        sim.tree = sim.tree
    ))
    toreturn$contact.matrix <- label_matrices
    saveRDS(toreturn, 
                sprintf("~/Documents/phd_files/contact/results/outbreaks/2types_block/contact_obsize20_intro1_2type_prop%s_probtrans%s_largeblock_%s.RDS", 
                        paste(prop, collapse="-"), paste(probtrans, collapse="-"), nr))
}
rm(list = setdiff(ls(), 'outbreak.list'))

lapply(seq_along(outbreak.list), function(i) {
    print(i)
    createBlockMatrix(outbreak.list[[i]], i)
})

infectors <- sapply(sim$sim.infectors, function(x){
            if (x == "index") return(0)
            else return(as.numeric(substr(x, 6, nchar(x))))
        })
inf.remove <- which(infectors==0)
infectors <- infectors[-inf.remove]
infectees <- (1:20)[-inf.remove]  
cnt <- sapply(sim$contact.matrix, function(m){ 
    sum(sapply(seq_along(infectees), function(host){
    m[host,infectors[host]]
}))})
print(cnt/19)


files <- files %>% 
  gsub("\\D", "", .) %>%               # Remove non-numeric parts
  as.numeric() %>%                     # Convert to numeric
  order() %>%                          # Get sorted indices
  (\(x) files[x])()

lapply(seq_along(files), function(n){
    sim <- readRDS(sprintf("results/outbreaks/2types_swapped/%s",files[n]))
    sim$contact.matrix <- sim$contact.matrix[2:1]
    saveRDS(sim, sprintf("results/outbreaks/2types_swapped/contact_obsize20_intro1_2type_prop0.1-0.1_probtrans0.2-0.7_%s.RDS",n))
})

lapply(1:10, function(n){
    s1 <- readRDS(sprintf("results/outbreaks/2types_swapped/contact_obsize20_intro1_2type_prop0.1-0.1_probtrans0.2-0.7_%s.RDS", n))
    s2 <- readRDS(sprintf("results/outbreaks/2types_swapped/contact_obsize20_intro1_2type_prop0.1-0.1_probtrans0.7-0.2_%s.RDS", n))
    all(s1$contact.matrix[[1]] == s2$contact.matrix[[2]])
})

## Remove 20% of the contact matrix entries
files <- list.files("results/outbreaks_numbers/2types", pattern = "contact_obsize20_intro1_2type_prop0.1-0.1_probtrans15-4_\\d+\\.RDS$", full.names = FALSE)
lapply(seq_along(files), function(n){
    sim <- readRDS(sprintf("results/outbreaks_numbers/2types/%s",files[n]))
    sim$contact.matrix <- lapply(sim$contact.matrix, function(mat) {
        select_random_matrix_entries(mat, prop = 0.2)
    })
    saveRDS(sim, sprintf("results/outbreaks_numbers/2types_reportcov0.8/%s",files[n]))
})
