devtools::load_all("~/Documents/phybreak/R")

outbreak.list <- lapply(11:25, function(i){
    sim <- sim_phybreak(obsize = 20, introductions = 1)
    saveRDS(sim, sprintf("~/Documents/phd_files/contact/results/outbreaks/simulation_outbreak_nr%s.Rds", i))
    return(sim)
})
outbreak.list <- list.files("~/Documents/phd_files/contact/results/outbreaks/blank", full.names = TRUE) %>%
  gsub("\\D", "", .) %>%               # Remove non-numeric parts
  as.numeric() %>%                     # Convert to numeric
  order() %>%                          # Get sorted indices
  (\(x) list.files("~/Documents/phd_files/contact/results/outbreaks/blank", full.names = TRUE)[x])() %>%
  lapply(readRDS)

vec = c(0, 0.5, 1)
params <- data.frame()
for (i in vec[-1]){
    for (j in vec){
        params <- rbind(params, data.frame(prop = i, prob.trans = j))
    }
}

params <- data.frame(p = rep(c(0.1, 0.5, 1),3), x = rep(c(0, 10, 19),each=3))
params <- data.frame(p = rep(c(0.01, seq(0.05, 0.5, length.out = 10)), 5), x = rep(c(1, 5, 10, 15, 19), each = 11))


for (nr in seq_along(outbreak.list)){
    print(sprintf("Outbreak nr %s", nr))
    sim <- outbreak.list[[nr]]
    for (i in seq_len(nrow(params))){
        print(sprintf("Param set %s of %s", i, nrow(params)))
        #print(params[i,])
        infectors <- sapply(sim$sim.infectors, function(x){
            if (x == "index") return(0)
            else return(as.numeric(substr(x, 6, nchar(x))))
        })
        obs=20

        p <- params[i,1]
        x <- params[i,2]

        s <- sim_contact_matrix(list(sim), x, p)
        toreturn <- with(s[[1]], phybreakdata(
            sequences = sequences,
            sample.times = sample.times,
            sample.names = names(sample.times),
            host.names = sample.hosts,
            sim.infection.times = sim.infection.times ,
            sim.infectors = sim.infectors ,
            sim.tree = sim.tree
        ))
        toreturn$contact.matrix <- s$contact.matrix
        saveRDS(toreturn, sprintf("~/Documents/phd_files/contact/results/outbreaks_numbers/1type/contact_obsize20_intro1_1type_prop%s_probtrans%s_%s.RDS", 
                                  p, x, nr))
    }
}
