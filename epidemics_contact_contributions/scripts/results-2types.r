library(phybreak)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)

dir <- "~/Documents/phd_files/contact/results/results_numbers/2types"
files <- list.files(dir)
files <- files[grep("jeffreys", files)]

dir <- "~/Documents/phd_files/contact/results/results_numbers/2types"
files2 <- list.files(dir)
files2 <- files2[grepl("jeffreys", files2)]
files2 <- files2[!grepl("probtrans5-5", files2)]
files2 <- files2[c(1:10,21:30)]

nkeep = 5e3

# Function to extract numeric values from filenames
extract_numeric_values <- function(split, index, prefix_length) {
  as.numeric(unlist(strsplit(substring(split[index], prefix_length, nchar(split[index])), "-")))
}

# Function to calculate type probability estimates
calculate_type_prob_est <- function(contact.coeff, contact.prop, r) {
  do.call(rbind, lapply(1:2, function(k){
    sapply(seq_along(contact.coeff[k,]), function(j){
      return(contact.coeff[k,j] * contact.prop[k,j] / (r[j] + sum(contact.coeff[,j] * contact.prop[,j])))
    })
  }))
}

res <- do.call(rbind, lapply(seq_along(files), function(i){
    print(i)
    f <- files[i]
    split <- unlist(strsplit(f, "_"))

    prop <- extract_numeric_values(split, 5, 5)
    probtrans <- extract_numeric_values(split, 6, 10)
    n <- as.numeric(strsplit(strsplit(f, "_")[[1]][7], "\\.")[[1]][1])

    if(grepl("^informative", f)) prior <- "informative"
    else if (grepl("^jeffreys", f)) prior <- "jeffreys"
    else if (grepl("^uniform", f)) {
      prior <- "uniform"
      return(NULL)
    } else if (grepl("^strong", f)) {
      prior <- "strong_informative"
      return(NULL)
    } else prior <- NA
    print(prior)
    s <- readRDS(sprintf("%s/%s", dir, f))

    if (is.null(s$s$contact.fracs)) return(NULL)

    r <- tail(s$s$r, nkeep)
    #contact.coeff <- s$s$contact.coeff[, (ncol(s$s$contact.coeff)-nkeep+1):ncol(s$s$contact.coeff)]
    #contact.prop <- s$s$contact.prop[, (ncol(s$s$contact.prop)-nkeep+1):ncol(s$s$contact.prop)]
    #type.prob.est <- calculate_type_prob_est(contact.coeff, contact.prop, r)
    est <- s$s$contact.fracs[, (ncol(s$s$contact.fracs)-nkeep+1):ncol(s$s$contact.fracs)]

    sim <- readRDS(sprintf("~/Documents/phd_files/contact/results/outbreaks/blank/simulation_outbreaks_nr%s.RDS",n ))
    #sim <- readRDS(sprintf("~/Documents/phd_files/contact/results/outbreaks_numbers/2types/%s", f))
    inf.set <- infectorsets(s, output='matrix')
    accuracy <- sum(transtree(s)$infector == sim$sim.infectors)
    #print(accuracy)

    # corr = ifelse(prop[4] == 0.0, "anti-correlated", "correlated")
    # if (corr == "correlated" & prop[1] == 0.5) corr = "correlated_unequal"
    # if(grepl("block", f)){
    #   block <- ifelse(grepl("equal", f), "equal", "large")
    # } else {
    #   block <- NA
    # }

    return(data.frame(prior = prior,
                      prop1 = prop[1],
                      prop2 = prop[2], 
                      probtrans1 = probtrans[1],
                      probtrans2 = probtrans[2], 
                      nr = n, 
                      est1 = est[1,], 
                      est2 = est[2,],
                      intro = tail(s$s$introductions, nkeep),
                      #c1 = contact.coeff[1,], 
                      #c2 = contact.coeff[2,],
                      #r = r, 
                      #block = block,
                      accuracy = accuracy))
    # return(data.frame(prop1 = prop[2]+prop[4],
    #                   prop2 = prop[3]+prop[4], 
    #                   probtrans1 = probtrans[1],
    #                   probtrans2 = probtrans[2], 
    #                   nr = n, 
    #                   est1 = est[1,], 
    #                   est2 = est[2,],
    #                   accuracy = accuracy,
    #                   corr = corr))
}))

# Function to create a ggplot for each plot
create_plot <- function(plot_data, type) {
  if (type == "prob"){
    plot_data$subplot <- rep(LETTERS[1:length(unique(plot_data$prop1))], times = (plot_data %>% group_by(prop1, prop2) %>% summarize(n=n()))$n)

    ggplot(plot_data %>% pivot_longer(cols = c(probtrans1, probtrans2), names_to = "trans_name", values_to = "trans_value"), 
      aes(x = value, fill = name)) +
      geom_histogram(position=position_dodge()) +
      geom_vline(aes(xintercept = trans_value, color = trans_name), linetype = 'dashed') +
      scale_color_discrete(guide = 'none') +
      scale_fill_discrete(labels = c("Route 1", "Route 2")) +
      facet_wrap(prop2~ prop1, scales = 'fixed', labeller = label_both, nrow = 1) +
      ggtitle(sprintf("x1 = %s, x2 = %s", unique(plot_data$probtrans1), unique(plot_data$probtrans2))) +
      labs(x = "Estimated expected fraction of transmission via route", y = "Frequency", 
         fill = "Transmission routes", color = "") +
      theme_bw()
  } else if (type == "accuracy") {
    ggplot(plot_data, aes(x = accuracy)) +
      geom_histogram(position = position_dodge()) +
      facet_wrap(probtrans2 ~ probtrans1, scales = 'fixed') +
      #scale_x_continuous(sec.axis = sec_axis(~ . , name = "Proportion of transmission via route 1", labels = NULL, breaks = NULL)) +
      #scale_y_continuous(sec.axis = sec_axis(~ . , name = "Proportion of transmission via route 2", labels = NULL, breaks = NULL)) +
      ggtitle(sprintf("p1 = %s, p2 = %s", unique(plot_data$prop1), unique(plot_data$prop2))) +
      theme_bw()
  }
}

# Function to print final plot
print_plot <- function(data, type = 'prob'){
    if (type == 'prob'){
      data$plot <- rep(LETTERS[1:length(unique(data$probtrans1))], times = (data %>% group_by(probtrans1, probtrans2) %>% summarize(n=n()))$n)

#      data$plot <- c(rep("A",2*2.8e5), rep(LETTERS[2:4], each = 2*2e5))
      plots <- split(data, data$plot)
      plot_list <- lapply(plots, create_plot, type = 'prob')
    } else if (type == 'accuracy'){
      data$plot <- rep(LETTERS[1:4], each = 30)
      plots <- split(data, data$plot)
      plot_list <- lapply(plots, create_plot, type = "accuracy")
    }
    final_plot <- plot_grid(plotlist = plot_list, nrow = length(unique(data$plot)))
    print(final_plot)
}

# Split data by plot identifier and print plot
print_plot(res %>% 
           filter(probtrans1 %in% c(0.2, 0.7, 1) & probtrans2 %in% c(0, 0.2)) %>% 
           arrange(probtrans1, probtrans2) %>%
           pivot_longer(cols = c(est1,est2)), type = 'prob')
print_plot(res %>%
        group_by(prop1, prop2, probtrans1, probtrans2, nr) %>%
        summarize(accuracy = mean(accuracy)), type = 'accuracy')

### plot barplots from 0 to 1 giving fractions of transmission via routes

# order prior names
res$prior <- factor(res$prior, levels = c("simulation", "strong_informative", "informative", "jeffreys", "uniform"))

res.barplot <- res %>% filter(!(probtrans1 %in% c(2,3,5))) %>%
  # filter(prior != "informative") %>%
  group_by(prior, prop1, prop2, probtrans1, probtrans2)
  #summarize(est1 = mean(est1), est2 = mean(est2), mean_intro = median(intro, na.rm = T))

res.barplot <- rbind(res.barplot, data.frame(reportage = 'complete', prior = "simulation", 
                             prop1 = 0, prop2 = 0, 
                             probtrans1 = c(19,15,5,2,3,4), probtrans2 = c(0,4,3,2,0,0),
                             #nr = rep(1:25,9), 
                             est1 = c(1, 15/19, 5/19, 2/19, 3/19,4/19),
                             est2 = c(0, 4/19, 3/19, 2/19, 0,0), mean_intro = 1)) %>%
  mutate(est3 = 1-est1-est2) %>%
  pivot_longer(cols = c(est3, est2, est1), names_to = "route", values_to = "value") %>%
  mutate(value = value * 19)

res.barplot.sim <- res.barplot %>% filter(prior == "simulation") %>% 
mutate(mean = value)
res.barplot.est <- res.barplot %>% filter(prior != "simulation") %>%
  group_by(reportage, prior, prop1, prop2, probtrans1, probtrans2, route) %>%
  summarize(mean = mean(value), 
            lower = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upper = HPDinterval(as.mcmc(value), prob = 0.95)[,2])

res.barplot <- rbind(res.barplot.est, res.barplot.sim[-c(7:10, 12)])
  # summarize(value = mean(value), se = sqrt(mean(value)*(1-mean(value))/25), mean_intro = median(intro, na.rm=T)) %>%
  # mutate(ymin = cumsum(value) - value,  # Start of the bar
  #        ymax = ymin + value,            # End of the bar
  #        err_min = ymax - 1.96*se,         # Lower bound of error bar
  #        err_max = ymax + 1.96*se)

res.barplot$pb1pb2 <- factor(interaction(res.barplot$probtrans1, res.barplot$probtrans2), 
                             levels = c("19.0", "15.4", "5.3", "2.2", "3.0", "4.0"))
res.barplot$p1p2 <- factor(interaction(res.barplot$prop1, res.barplot$prop2), 
                            levels = c("0.1.0.1", "0.5.0.5", "1.1"))

# add simulated values
# res.sim <- data.frame(prior = "simulation", 
#                              prop1 = 0, prop2 = 0, 
#                              probtrans1 = rep(c(19,15,5), 3), probtrans2 = rep(c(0,4,3),3),
#                              #nr = rep(1:25,9), 
#                              est1 = rep(c(1, 15/19, 5/19),3),
#                              est2 = rep(c(0, 4/19, 3/19),3),
#                              accuracy = NA, intro = 1)
# res.sim <- res.sim %>% mutate(est0 = 1-est1-est2) %>%
#   pivot_longer(cols = c(est0, est1, est2), names_to = "route", values_to = "value")

# res.barplot <- rbind(res.barplot, res.sim)
vir.colors <- viridis::viridis(3)
names(vir.colors) <- c("est1", "est2", "est3")

res.barplot$route <- factor(res.barplot$route, levels = c("est3", "est2", "est1"))

trans.nr.barplot <- data.frame(res.barplot %>% filter(!(probtrans1 %in% c(5,2,3))) %>%  
  filter(!(prior == "simulation" & probtrans1 != 15)) %>% 
mutate(
    mean_scaled = mean,
    lower_scaled = lower,
    upper_scaled = upper,
    ymin_stack = cumsum(lag(mean_scaled, default = 0)),
    ymax_stack = ymin_stack + mean_scaled,
    lower_pos = ymin_stack - upper_scaled,
    upper_pos = ymin_stack + upper_scaled
  )) %>% mutate(upper_pos = pmin(upper_pos, 19), lower_pos = pmax(lower_pos, 0)) %>%  # Upper bound of error bar
  ggplot(aes(x = factor(prop1, levels = rev(c(0, 0.1, 0.5, 1))), y = mean_scaled, fill = route)) + 
  scale_y_continuous(limits = c(0, 20), breaks = c(0,5,10,15,20)) +
  coord_flip() +
  geom_bar(stat = 'identity') + 
  geom_errorbar(aes(ymin = ymax_stack, ymax = upper_pos, color = route),stat = 'identity', width = 0.2) +
  scale_fill_manual(labels = c("est1" = "Route 1", "est2" = "Route 2", "est3" = "Unknown route"),
                       values = vir.colors) +
  scale_color_manual(labels = c("est1" = "Route 1", "est2" = "Route 2", "est3" = "Unknown route"),
                       values = vir.colors) +
  labs(x = "Proportion of contact in population", y = "Number of transmissions via route", 
  fill = "Route", color = "Route") +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(reportage ~ 1, 
    labeller = labeller(reportage = c(complete = "Complete\ncontact data", incomplete = "Incomplete\ncontact data"))
  ) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Rotate x-axis text
    axis.text.y = element_text(size = 20),  # Y-axis text size
    axis.title = element_text(size = 20),  # Axis titles size
    legend.text = element_text(size = 20),  # Legend text size
    legend.title = element_text(size = 20),  # Legend title size
    plot.title = element_text(size = 20, face = "bold")  # Plot title size
  )
trans.nr.barplot

### correlation
### plot barplots from 0 to 1 giving fractions of transmission via routes
res.barplot <- res %>% filter(!(probtrans1==5 & probtrans2 == 5)) %>%
  # filter(prior != "informative") %>%
  group_by(prop1, prop2, probtrans1, probtrans2) %>%
  mutate(est3 = 1-est1-est2) %>%
  pivot_longer(cols = c(est3, est1, est2), names_to = "route", values_to = "value") %>%
  # mutate(value = value * 19) %>%
  group_by(prop1, prop2, probtrans1, probtrans2, route) %>%
  summarize(mean = mean(value), 
            lower = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upper = HPDinterval(as.mcmc(value), prob = 0.95)[,2])
# add simulated values
res.barplot <- rbind(res.barplot, data.frame(prop1 = 0, prop2 = 0, 
                             probtrans1 = rep(c(4,15), each = 3), probtrans2 = rep(c(0,4), each = 3),
                             route = c("est3", "est1", "est2"), 
                             mean = c(15/19, 4/19, 0/19, 
                                      0/19, 15/19, 4/19),
                             lower = 0, upper = 0))

# order correlation names
res$corr <- factor(res$corr, levels = c("simulation", "correlated", "correlated_unequal", "random", "anti-correlated"))

res.barplot$route <- factor(res.barplot$route, levels = c("est3", "est2", "est1"))
trans.nr.barplot <- res.barplot %>%       # Upper bound of error bar
  mutate(
    mean_scaled = mean * 19,
    lower_scaled = lower * 19,
    upper_scaled = upper * 19,
    ymin_stack = cumsum(lag(mean_scaled, default = 0)),
    ymax_stack = ymin_stack + mean_scaled,
    lower_pos = ymin_stack + lower_scaled,
    upper_pos = ymin_stack + upper_scaled
  ) %>%
  ggplot(aes(x = factor(prop1, levels = rev(c(0, 0.8, 0.1, 0.9, 0.5))), y = mean_scaled, fill = route)) + 
  geom_bar(stat = 'identity') +
  #geom_errorbar(aes(ymin = ymax_stack, ymax = upper_pos, col = route), stat = 'identity', width = 0.2) +
  coord_flip() +
  scale_fill_manual(
    labels = c("est1" = "Route 1", "est2" = "Route 2", "est3" = "Unknown route"),
    values = vir.colors
  ) +
  scale_color_manual(
    labels = c("est1" = "Route 1", "est2" = "Route 2", "est3" = "Unknown route"),
    values = vir.colors
  ) +
  labs(
    y = "Proportion of contact in population",
    x = "Number of transmissions via route",
    fill = "Route", col = "Route"
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(~ factor(probtrans1, levels = c(15, 4))) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Rotate x-axis text
    axis.text.y = element_text(size = 20),  # Y-axis text size
    axis.title = element_text(size = 20),  # Axis titles size
    legend.text = element_text(size = 20),  # Legend text size
    legend.title = element_text(size = 20),  # Legend title size
    plot.title = element_text(size = 20, face = "bold")  # Plot title size
  )
trans.nr.barplot


res %>% group_by(prop1, prop2, probtrans1, probtrans2,corr) %>% summarize(est1 = median(est1), est2=median(est2),accuracy = mean(accuracy))

res %>% pivot_longer(cols = c(est1, est2)) %>% 
  ggplot(aes(x=value, fill = factor(interaction(name, prop1),  levels = c("est1.0.1", "est1.0.5", "est2.0.1", "est2.0.5")))) + 
  #geom_histogram(position = position_dodge()) + 
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("est1.0.1" = 'red', "est1.0.5" = 'lightcoral', "est2.0.1" = 'blue', 'est2.0.5'= 'skyblue'),
                    labels = c("Route 1 - equal contact proportion", "Route 1 - different contact proportion", "Route 2 - equal contact proportion", "Route 2 - different contact proportion")) +
  facet_grid(~ corr) + 
  labs(x = "Estimated expected fraction of transmission via contact route",
       y = "Density",
       fill = "Transmission route") +
  theme_bw()

res %>% pivot_longer(cols = c(est1, est2)) %>% filter(corr == 'correlated') %>%
  ggplot(aes(x=value, fill = name)) + 
  geom_histogram(position = position_dodge()) + 
  #geom_density(alpha = 0.5) +
  # scale_fill_manual(values = c("est1.0.1" = 'red', "est1.0.5" = 'lightcoral', "est2.0.1" = 'blue', 'est2.0.5'= 'skyblue'),
  #                   labels = c("Route 1 - equal contact proportion", "Route 1 - different contact proportion", "Route 2 - equal contact proportion", "Route 2 - different contact proportion")) +
  scale_fill_discrete(labels = c("Route 1", "Route 2")) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Proportion of contact", labels = NULL, breaks = NULL)) +
  facet_grid(~ prop1) + 
  labs(x = "Estimated expected fraction of transmission via contact route",
       y = "Frequency",
       fill = "Transmission route") +
  theme_bw() + 
  theme(text = element_text(size = 20))

blank_files <- list.files("~/Documents/phd_files/contact/results/results_numbers/blank", full.names = T)
blank_files <- blank_files[grep("blank_contact",blank_files)]

blanks <- do.call(rbind, lapply(blank_files, function(f){
    #n <- as.numeric(substr(strsplit(strsplit(f, "_")[[1]][5], "\\.")[[1]][1],3,nchar(strsplit(strsplit(f, "_")[[1]][5], "\\.")[[1]][1])))
    n <- strsplit(f, "_")[[1]][6]
    n <- as.numeric(substr(n, 3, nchar(n)-4))
    s <- readRDS(f)
    sim <- readRDS(sprintf(sprintf("~/Documents/phd_files/contact/results/outbreaks/blank/simulation_outbreaks_nr%s.Rds",n )))
    # inf.set <- infectorsets(s, output='matrix')
    # accuracy <- sapply(seq_along(sim$sim.infectors), function(inf){
    #     infector <- sim$sim.infectors[inf]
    #     infectee <- names(sim$sim.infectors)[inf]
    #     support <- inf.set[rownames(inf.set) == infector, colnames(inf.set) == infectee]
    # })#)/length(sim$sim.infectors)
    accuracy = sum(transtree(s)$infector == sim$sim.infectors)
    return(data.frame(n = n, accuracy = accuracy))
}))
blanks <- blanks %>% 
    summarize(percentage = sum(accuracy)/500) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/500))


# Plot accuracy using ggplot for the combinations prop1, prop2, probtrans1, probtrans2 found in res
res <- res %>% filter(!(probtrans1 %in% c(5,2,3)))  %>% 
    filter(prior %in% c( 'informative','jeffreys')) %>% 
    group_by(prop1, prop2, probtrans1, probtrans2, prior) %>%
    summarize(percentage = sum(accuracy)/500) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/500))


pos <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.2)
colnames(res)[2] <- 'trans'
shapes.types <- c("19.0" = 0, "15.4" = 1, "4.0" = 2)
p2_acc_plot <- res  %>%
    ggplot(aes(x = factor(prop1), y = percentage, shape = factor(interaction(probtrans1, probtrans2), levels = c("19.0", "15.4", "4.0")))) + #group = interaction(factor(trans), factor(noGen)))) +
    geom_hline(yintercept = blanks$percentage, col = 'grey50', linetype = 'dashed') +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = blanks$percentage-1.96*blanks$se, ymax = blanks$percentage+1.96*blanks$se,
                   fill = "grey70", color = "grey70", alpha = 0.03) +
    geom_jitter(position = pos, size = 5) +
    geom_errorbar(aes(ymin = percentage-1.96*se, ymax = percentage+1.96*se, width = 0), position = pos) +
    #geom_line(linetype='dashed') +
    scale_y_continuous(limits = c(0.5,1)) +
    scale_shape_manual(values = shapes.types) +
    #facet_grid(~ factor(prior, levels = c("strong_informative", "informative", "jeffreys", "uniform"))) +
    #scale_x_continuous(sec.axis = sec_axis(~ . , name = "Expected fraction of transmission via known route", labels = NULL, breaks = NULL)) +
    labs(x = "Proportion of contact", y = "Fraction correctly\nidentified infectors",
         shape = "Simulated number of\ntransmission via route") +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    theme(legend.position = "top")
p2_acc_plot

p12_acc <- cowplot::plot_grid(p1_acc_plot, p2_acc_plot, nrow = 1, align = "h", labels = c("A", "B"), label_size = 20)
pdf("figures/p12_acc.pdf", width = 15, height = 8)
print(p12_acc)
dev.off()
