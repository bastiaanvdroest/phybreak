devtools::load_all("~/Documents/phybreak/R")
library(dplyr)
library(ggplot2)
library(viridis)

## 1 contact type, prop transmission plots:

dir <- "~/Documents/phd_files/contact/results/results_numbers/1type"
files <- list.files(dir)
files <- files[grepl("jeffreys.*.RDS", files)]
files <- files[grepl("jeffreyscontact.*.RDS", files)]

nkeep = 5e3

res <- do.call(rbind, lapply(seq_along(files), function(i){
    print(i)
    f <- files[i]

    prop <- strsplit(strsplit(f, "prop")[[1]][2], "_")[[1]][1]
    # if (prop == "01") prop <- 0.1
    # else if (prop == "05") prop <- 0.5
    # else if (prop == "09") prop <- 0.9
    prop <- as.numeric(prop)
    #if (prop == 0 | prop == 1) return(data.frame())

    probtrans <- strsplit(strsplit(f, "probtrans")[[1]][2], "_")[[1]][1]
    # if (probtrans == "01") probtrans <- 0.1
    # else if (probtrans == "05") probtrans <- 0.5
    # else if (probtrans == "09") probtrans <- 0.9
    probtrans <- as.numeric(probtrans)

    if(grepl("^informative", f)) prior <- "informative"
    else if (grepl("jeffreys", f)) prior <- "jeffreys"
    else if (grepl("uniform", f)) prior <- "uniform"
    else if (grepl("strong", f)) prior <- "strong_informative"
    else prior <- NA

    noGen <- ifelse (grepl("noGen", f), 1, 0) 

    split <- strsplit(f, "_")[[1]]
    n <- as.numeric(substr(tail(split,1), 1, nchar(tail(split,1))-4))

    s <- readRDS(sprintf("%s/%s", dir, f))
    r <- tail(s$s$r, nkeep)
    est <- tail(s$s$contact.fracs[1,], nkeep)
    if(is.null(est)) est <- NA
    intro <- tail(s$s$introductions, nkeep)
    # contact.prop <- tail(s$s$contact.prop[1,], nkeep )

    # type.prob.est <- (contact.coeff * contact.prop) / (r + contact.coeff * contact.prop)
    
    sim <- readRDS(sprintf("~/Documents/phd_files/contact/results/outbreaks/blank/simulation_outbreaks_nr%s.Rds",n ))
    # inf.set <- infectorsets(s, output='matrix')
    # accuracy <- sum(sapply(seq_along(sim$sim.infectors), function(inf){
    #     infector <- sim$sim.infectors[inf]
    #     infectee <- names(sim$sim.infectors)[inf]
    #     support <- inf.set[rownames(inf.set) == infector, colnames(inf.set) == infectee]
    # }))/length(sim$sim.infectors)
    accuracy = sum(transtree(s)$infector == sim$sim.infectors)
    # inf <- as.numeric(do.call(rbind, strsplit(sim$sim.infectors, "host."))[,2])
    # index <- which(is.na(inf))
    # coeff.mat <- sum(sapply(c(1:20)[-index], function(j) sim$contact.matrix[j, inf[j]])) / 20
    
    # prev <- sum(colSums(s$d$contact)) / 
    #     ncol(s$d$contact)^2

    # p1 = coeff.mat * (prev - prop) / (prev * (coeff.mat - prop))

    return(data.frame(prop = prop, probtrans = probtrans, prior = prior, noGen = noGen,
                      nr = n, 
                      # r = r,
                      # est1 = est,
                      # ntro = intro))#, 
                      accuracy = accuracy))
}))
# res %>% 
#     mutate(diff = abs(est - probtrans)) %>%
#     group_by(prop, probtrans, nr) %>%
#     summarize(diff.mean = mean(diff)) %>%
#     filter(prop == 0.5)

# res <- res %>% filter(prop > 0 & prop < 1) %>% filter(probtrans < 1)    

# hist_data <- do.call(rbind, lapply(unique(res$prop), function(prop){
#     do.call(rbind, lapply(unique(res$probtrans), function(probtrans){
#         est <- res$est[res$prop==prop & res$probtrans==probtrans]
#         accuracy <- res$accuracy[res$prop==prop & res$probtrans==probtrans]
#         breaks <- seq(min(est), max(est), length.out = 50)
#         bin <- cut(est, breaks = breaks)
#         bin_center <- (breaks[-1] + breaks[-length(breaks)]) / 2
#         bin_center <- bin_center[as.numeric(bin)]
#         return(data.frame(x = probtrans, y = prop, bin = bin, bin_center = bin_center))
#     }))
# }))
# hist_data <- hist_data %>% group_by(x, y, bin, bin_center) %>% summarize(n=n()) 

# dist.plot <- ggplot(res, aes(x = est)) +
#     geom_histogram(position=position_dodge()) +
#     #geom_vline(aes(xintercept = coeff.mat)) +
#     facet_grid(prop ~ probtrans, scales = 'fixed') +
#     scale_x_continuous(sec.axis = sec_axis(~ . , name = "Proportion of transmission via contact route", labels = NULL, breaks = NULL)) +
#     scale_y_continuous(sec.axis = sec_axis(~ . , name = "Proportion of contact route in population", labels = NULL, breaks = NULL)) +
#     theme_bw()
# dist.plot

# density_plot1 <- ggplot(res %>% filter(prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 10, 19)), aes(x = est * 19, fill = factor(prop), color = factor(prop))) +
#   geom_density(alpha = 0.4, linewidth = 1.2) +  # Transparent fill + density lines
# #   geom_vline(data = aggregate(est ~ probtrans + prop, res %>% filter(prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 0.5, 1)), median), 
# #              aes(xintercept = est, color = prop), 
# #              linetype = "dashed", size = 1) +  # Dashed line for mean
#   labs(
#     title = "",
#     x = "Estimated expected fraction of transmissions via contact route",
#     y = "Density",
#     fill = "Proportion\nof contacts",
#     col = "Proportion\nof contacts"
#   ) +
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "Expected fraction of transmission via contact route", labels = NULL, breaks = NULL)) +
#   facet_wrap(prior ~ probtrans, ) +
#   theme_bw()

density_plot2 <- res %>% filter(prior == "jeffreys" & prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 10, 19)) %>%
  filter(!(nr %in% c(2,4,8) & prop == 0.1)) %>%
  ggplot(aes(x = r, fill = factor(probtrans), color = factor(probtrans))) +
    #geom_density(alpha = 0.4, linewidth = 1.2) +  # Transparent fill + density lines
    geom_histogram(position = position_dodge()) +
  #  geom_vline(data = aggregate(est ~ probtrans + prop, res %>% filter(prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 0.5, 1)), median), 
  #              aes(xintercept = est, color = prop), 
  #              linetype = "dashed", size = 1) +  # Dashed line for mean
    labs(
      title = "",
      x = "Reproduction number",
      y = "Density",
      fill = "Simulated number\nof transmissions",
      col = "Simulated number\nof transmissions"
    ) +
    scale_x_continuous(sec.axis = sec_axis(~ . , name = "Proportion of contacts", labels = NULL, breaks = NULL)) +
    facet_grid( ~ prop, scales = 'free_y') +
    theme_bw()

res$prior <- factor(res$prior, levels = c("strong_informative", "informative", "jeffreys", "uniform"))
res <- res[!is.na(res$est),]
scatter_plot <- res %>% filter(prior %in% c("informative","jeffreys") & prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 10, 19) & noGen == 0) %>%
  #filter(!(nr %in% c(2,4,10))) %>%
  group_by(prop, probtrans, prior) %>%
  summarize(mean = mean(est), mean_intro = median(ntro)) %>%
  mutate(se = sqrt(mean*(1-mean)/1400)) %>%
  ggplot(aes(x = factor(prop), y = mean*(20-mean_intro) , color = factor(probtrans), group = factor(probtrans))) +
  geom_point() +
  geom_line(linetype = 'dashed') +
  geom_errorbar(aes(ymin = (20-mean_intro)*(mean - 1.96*se), ymax = (20-mean_intro)*(mean + 1.96*se)), width = 0.2) +
  geom_hline(yintercept = 19, linetype = 'dotted') +
  scale_y_continuous(limits = c(0, 19), breaks = c(0,5,10,15,20)) +
#  facet_grid(~ prior) +
  labs(x = "Proprotion of contact", y = "Estimated number of\ntransmissions via route", 
       color = "Simulated number of\ntransmissions via route")+
  theme_bw() + 
  theme(text = element_text(size = 20), legend.position = 'none')
scatter_plot

res.mean <- res %>% filter(!(prop %in% c(0.2, 0.3, 0.4))) %>%
  # filter(prior != "informative") %>%
  group_by(prop, probtrans) %>%
  mutate(est3 = 1-est1) %>%
  pivot_longer(cols = c(est1, est3), names_to = "route", values_to = "value") %>%
  group_by(prop, probtrans, route) %>%
  summarize(mean = mean(value), 
            lower = HPDinterval(as.mcmc(value), prob = 0.95)[,1],
            upper = HPDinterval(as.mcmc(value), prob = 0.95)[,2])

res.mean <- rbind(res.mean, data.frame(prop = 0, 
                                       probtrans = rep(c(0, 10, 19), each = 2),
                                       route = rep(c("est1", "est3"), 3),
                                       mean = c(0, 1, 10/19, 9/19, 1, 0)))



res.mean$route <- factor(res.mean$route, levels = c("est3","est1"))  

vir.colors <- viridis::viridis(3)
names(vir.colors) <- c("est1", "est2", "est3")

barplot <- res.mean %>%       # Upper bound of error bar
  mutate(
    mean_scaled = mean * 19,
    lower_scaled = lower * 19,
    upper_scaled = upper * 19,
    ymin_stack = cumsum(lag(mean_scaled, default = 0)),
    ymax_stack = ymin_stack + mean_scaled,
    lower_pos = ymin_stack + lower_scaled,
    upper_pos = ymin_stack + upper_scaled
  ) %>%
  ggplot(aes(x = factor(prop, levels = rev(c(0, 0.1, 0.5, 1))), y = mean_scaled, fill = route)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0, 20), breaks = c(0,5,10,15,20)) +
  coord_flip() +
  geom_errorbar(aes(ymin = ymax_stack, ymax = upper_pos, col = route),stat = 'identity', width = 0.2) +
  scale_fill_manual(values = vir.colors, labels = c("est3" = "Unknown route", "est1" = "Route 1")) +
  scale_color_manual(values = vir.colors, labels = c("est3" = "Unknown route", "est1" = "Route 1")) +
  labs(x = "Proportion of contact in population", y = "Number of transmissions via route", fill = "Route", col = "Route") +
  guides(fill = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE)) +
  facet_grid( ~ probtrans) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Rotate x-axis text
    axis.text.y = element_text(size = 20),  # Y-axis text size
    axis.title = element_text(size = 20),  # Axis titles size
    legend.text = element_text(size = 20),  # Legend text size
    legend.title = element_text(size = 20),  # Legend title size
    legend.position = "bottom",
    plot.title = element_text(size = 20, face = "bold")  # Plot title size
  )


res %>% filter(!is.na(prior) & prior != "strong_informative" & prop %in% c(0.1, 0.5, 1) & probtrans %in% c(0, 10, 19)) %>%
  #filter(!(nr %in% c(4,8,10))) %>%
  group_by(prop, probtrans, prior) %>%
  summarize(mean = median(est), mean_intro = median(intro)) %>%
  ggplot(aes(x = mean_intro)) %>%
  geom_point() +
  geom_line(linetype = 'dashed')

x <- seq(0,2,by=0.01)
a <- 1:19
p <- 1
#q <- sapply(1:18, function(n) n/19)
q <- seq(0,1,by=0.01)[-101]
#q <- 11/19

max_list <- vector("list", 18)
for (k in seq_along(a)){
  vec_list <- vector("list", length(q))
  for(i in seq_along(q)){
    for (j in seq_along(x)){
      vec_list[[i]] <- c(vec_list[[i]], a[k]*log(x[j] + x[j] * (q[i]/(p-p*q[i]))) + (19-a[k])*log(x[j]) - 20*(x[j]+x[j]*p*(q[i]/(p-p*q[i]))))
    }
  }
  max <- sapply(vec_list, function(y) y[which(y == max(y[!is.nan(y)]))])
  max_list[[k]] <- max
}
# max <- sapply(vec_list, function(y) y[which(y == max(y[!is.nan(y)]))])
# which(max(max) == max)
# y <- c()
# for (i in seq_along(a)){
#   y <- c(y, 19*log(max[i] + a[i]*max[i]) + 0*log(max[i]) - 19*(max[i]+a[i]*max[i]))
# }
# plot(1:18, y)

pdf(sprintf("figures/1type_likelihood_p%s.pdf",p), width = 15, height = 15)
par(mfrow = c(2,1))
plot(NULL, xlim = c(0, 19), ylim = range(unlist(max_list)[!is.infinite(unlist(max_list))]), 
     xlab = "Estimated number of transmission via route", ylab = "logLikelihood", main = sprintf("Proportion of contact = %s",p))
for (i in seq_along(max_list)) {
  lines(q*19, max_list[[i]], type = "o", col = i, pch = i)
}
legend("bottomleft", legend = sprintf("a = %s", a), col = 1:length(vec_list), pch = 1:length(vec_list), lty = 1)

metamax <- sapply(max_list, function(y){
  y <- unlist(y)
  return(q[which(y == max(y[!is.nan(y)]))])
})
plot(a,metamax*19,xlab = "Simulated number of transmissions via route", ylab="Most likely number of transmissions via route")
dev.off()

accuracy.plot <- ggplot(res %>% group_by(prop,  probtrans) %>% summarize(est = mean(est), accuracy = mean(accuracy)), 
      aes(x = factor(probtrans), y = factor(prop, levels = rev(unique(prop))), fill = accuracy)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', high='red',mid='white',midpoint=mean(blanks$accuracy)) +
    theme_minimal()

ggplot(res) +
    geom_histogram(aes(x = est, fill = accuracy), stat = "count", bins = 10) +
    facet_grid(prop ~ probtrans)

cowplot::plot_grid(dist.plot, accuracy.plot, ncol = 2)

acc <- res %>% group_by(probtrans,  prop) %>% summarize(accuracy = mean(accuracy))
colnames(acc) <- c("x", "y", "accuracy")
hist_data <- left_join(hist_data, acc, by = c("x", "y"))

tile_hist_plot <- ggplot(data = hist_data %>% group_by(x,y,bin,bin_center,accuracy) %>% summarize(n = mean(n)), 
                         aes(x = bin_center, y = n, fill= accuracy)) +
    # geom_tile(data = res %>% group_by(prop,  probtrans) %>% summarize(est = mean(est), accuracy = mean(accuracy)),
    #           aes(x = factor(probtrans), y = factor(prop, levels = rev(unique(prop))), fill = accuracy)) +
    # scale_fill_gradient2(low='blue', high='red',mid='white',midpoint=mean(blanks$accuracy)) +
    geom_bar(stat = 'identity',position = 'identity') +
    facet_grid(y ~ x, scales = 'free') +
    scale_fill_gradient2(low='blue', high='red',mid='grey80',midpoint=mean(blanks$accuracy)) +
        scale_x_continuous(sec.axis = sec_axis(~ . , name = "Probability of transmission", labels = NULL, breaks = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = "Proportion", labels = NULL, breaks = NULL)) +
    theme_bw()

### function to calculate number of transmissions with contact
contact_transmission <- sapply(seq_along(s$s$introductions), function(cycle){
  infectors <- s$s$infector[,cycle]
  fromto <- data.frame(from = 1:20, to = infectors)
  fromto <- fromto[fromto$to != 0,]
  return(sum(sapply(1:nrow(fromto), function(j) s$d$contact[fromto$from[j], fromto$to[j]])))
})


## 1 contact type, accuracy

dir <- "~/Documents/phd_files/contact/results/1type_accuracy/"
files <- list.files(dir)
nkeep = 1000

res <- do.call(rbind, lapply(seq_along(files), function(i){
    print(i)
    f <- files[i]
    prop <- strsplit(strsplit(f, "prop")[[1]][2], "_")[[1]][1]
    prop <- as.numeric(prop)

    trans <- strsplit(strsplit(f, "probtrans")[[1]][2], "_")[[1]][1]
    trans <- as.numeric(trans)

    n <- as.numeric(strsplit(strsplit(f, "_")[[1]][7], "\\.")[[1]][1])

    s <- readRDS(sprintf("%s/%s", dir, f))
    s <- phybreak:::thin.phybreak(s, nkeep = nkeep)

    sim.files <- list.files("~/Documents/phd_files/contact/results/outbreaks/1type_accuracy")
    if (f %in% sim.files){
        sim <- readRDS(sprintf("~/Documents/phd_files/contact/results/outbreaks/1type_accuracy/%s",f ))
        # inf.set <- infectorsets(s, output='matrix')
        # accuracy <- sapply(seq_along(sim$sim.infectors), function(inf){
        #     infector <- sim$sim.infectors[inf]
        #     infectee <- names(sim$sim.infectors)[inf]
        #     support <- inf.set[rownames(inf.set) == infector, colnames(inf.set) == infectee]
        # })#)/length(sim$sim.infectors)
        accuracy = sum(transtree(s)$infector == sim$sim.infectors)
        return(data.frame(prop = prop, trans = trans, nr = n, accuracy = accuracy))
    } else {
        return(data.frame(prop = NULL, trans = NULL, nr = NULL, accuracy = NULL))
    }
    
    
}))

res <- res %>% group_by(prop, trans) %>% summarize(low = quantile(accuracy, 0.025),
                                     med = quantile(accuracy, 0.5),
                                     mean = mean(accuracy),
                                     high = quantile(accuracy, 0.975))

res.inf <- res %>% filter(!is.na(prior)) %>% filter(prop %in% c(0.1, 0.5, 1)) %>% 
    filter(prior %in% c( 'informative')) %>% 
    group_by(prop, probtrans, prior, noGen) %>%
    summarize(percentage = sum(accuracy)/200) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/200))

res.jeff <- res %>% filter(!is.na(prior)) %>% filter(prop %in% c(0.1, 0.5, 1)) %>% 
    filter(prior %in% c( 'jeffreys')) %>% 
    group_by(prop, probtrans, prior, noGen) %>%
    summarize(percentage = sum(accuracy)/500) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/500))
res <- rbind(res.inf, res.jeff)

res.yesGen <- res %>% filter(!is.na(prior)) %>% filter(prop %in% c(0.1, 0.5, 1)) %>% 
    filter(prior %in% c( 'jeffreys'), noGen == 0) %>% 
    group_by(prop, probtrans, prior, noGen) %>%
    summarize(percentage = sum(accuracy)/500) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/500))

res.noGen <- res %>% filter(!is.na(prior)) %>% filter(prop %in% c(0.1, 0.5, 1)) %>% 
    filter(prior %in% c( 'jeffreys'), noGen == 1) %>% 
    group_by(prop, probtrans, prior, noGen) %>%
    summarize(percentage = sum(accuracy)/200) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/200))

res <- rbind(res.yesGen, res.noGen)

pos <- position_jitterdodge(jitter.width = 0.1, dodge.width = 0.2)
colnames(res)[2] <- 'trans'
p3_acc_plot <- res %>% filter(prop %in% c(0.1, 0.5, 1) ) %>%
    ggplot(aes(x = factor(prop), y = percentage, shape = factor(trans), color = factor(noGen), group = interaction(factor(trans), factor(noGen)))) +
    geom_hline(yintercept = blanks$percentage, col = 'grey50', linetype = 'dashed') +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = blanks$percentage-1.96*blanks$se, ymax = blanks$percentage+1.96*blanks$se,
                   fill = "grey70", color = "grey70", alpha = 0.03) +
    geom_jitter(position = pos, size = 5) +
    geom_errorbar(aes(ymin = percentage-1.96*se, ymax = percentage+1.96*se), position = pos, width = 0.01) +

    #geom_line(linetype='dashed') +
    scale_y_continuous(limits = c(0.2,1), breaks = (2:10)/10) +
    scale_color_manual(labels = c("0" = "Genetic + contact data", "1" = "Contact data only"), values = c("1" = "grey50", "0" = "black")) +
    # facet_grid(~ factor(prior, levels = c("informative", "jeffreys"))) +
    #scale_x_continuous(sec.axis = sec_axis(~ . , name = "Expected fraction of transmission via known route", labels = NULL, breaks = NULL)) +
    labs(x = "Proportion of contact", y = "Fraction correctly\nidentified infectors",
         shape = "Simulated number of\ntransmission via route",
         color = "Data types used") +
    theme_bw() +
    theme(text = element_text(size = 20), legend.position = 'top')
p1_acc_plot

pdf("figures/1type_p_vs_x_combined_jeffreys.pdf", width = 15, height = 10)
cowplot::plot_grid(scatter_plot, p1_acc_plot, ncol = 2, rel_widths = c(1, 1.4))
dev.off()

ggplot(res, aes(x = factor(prop), y = accuracy)) +
    geom_boxplot() +
    labs(x = "Proportion", y = "Accuracy") +
    facet_grid(~ trans) +
    theme_minimal()

blank_files <- list.files("~/Documents/phd_files/contact/results/results_numbers/blank", full.names = T)
blank_files <- blank_files[grep("blank_contact",blank_files)]
blanks <- do.call(rbind, lapply(blank_files, function(f){
    #n <- as.numeric(substr(strsplit(strsplit(f, "_")[[1]][5], "\\.")[[1]][1],3,nchar(strsplit(strsplit(f, "_")[[1]][5], "\\.")[[1]][1])))
    n <- do.call(rbind, strsplit(f, "_"))[,6]
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
blanks %>% summarize(low = quantile(accuracy, 0.025),
                                     med = quantile(accuracy, 0.5),
                                     mean = mean(accuracy),
                                     high = quantile(accuracy, 0.975))
blanks <- blanks %>% 
    summarize(percentage = sum(accuracy)/(25*20)) %>% 
    mutate(se = sqrt(percentage*(1-percentage)/500))

df <- data.frame(est = (contact.coeff * 0.1) / (r + contact.coeff * 0.1), sim = (2.25 * contact.prop) / (1.5 + 2.25 * contact.prop))
ggplot(df, aes(x = est)) +
    geom_histogram() +
    geom_vline(xintercept = (2.25 * contact.prop) / (1.5 + 2.25 * contact.prop))
