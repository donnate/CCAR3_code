library(tidyverse)
library(ggplot2)
library(pracma)
library(tidyverse)
theme_set(theme_bw(base_size = 14))

 file_list <- list.files(path = "experiments/simulations/new_results", 
                         pattern = "2026_newest_RRR_efficient_resultsnew_normalized*", full.names = TRUE)
results <- bind_rows(lapply(file_list, read.csv))
# 
# 
# 
summ = results %>% group_by(n, p1, p2, r, r_pca,
                            nnzeros,
                            overlapping_amount, noise,
                            #lambda_opt,
                            method,
                            theta_strength,
                            normalize_diagonal,
                            prop_missing) %>%
  summarise(distance_tot_mean = median(distance_tot),
            distance_U_mean = mean(distance_U),
            distance_V_mean = mean(distance_V),
            distance_tot_q50 = quantile(distance_tot, 0.5, na.rm=TRUE),
            distance_tot_q75 = quantile(distance_tot, 0.75, na.rm=TRUE),
            distance_tot_q25 = quantile(distance_tot, 0.25, na.rm=TRUE),
            distance_tot_q975 = quantile(distance_tot, 0.975, na.rm=TRUE),
            distance_tot_q025 = quantile(distance_tot, 0.025, na.rm=TRUE),
            distance_tot_q90 = quantile(distance_tot, 0.9, na.rm=TRUE),
            distance_tot_q10 = quantile(distance_tot, 0.1, na.rm=TRUE),
            prediction_tot_mean= mean(prediction_tot),
            prediction_tot_q50 = quantile(prediction_tot, 0.5, na.rm=TRUE),
            prediction_tot_q25 = quantile(prediction_tot, 0.75, na.rm=TRUE),
            prediction_tot_q75 = quantile(prediction_tot, 0.25, na.rm=TRUE),
            distance_U_q50 = quantile(distance_U, 0.5, na.rm=TRUE),
            distance_U_q25 = quantile(distance_U, 0.75, na.rm=TRUE),
            distance_U_q75 = quantile(distance_U, 0.25, na.rm=TRUE),
            prediction_U_mean= mean(distance_U),
            prediction_U_q50 = quantile(distance_U, 0.5, na.rm=TRUE),
            prediction_U_q25 = quantile(distance_U, 0.75, na.rm=TRUE),
            prediction_U_q75 = quantile(distance_U, 0.25, na.rm=TRUE),
            distance_V_q50 = quantile(distance_V, 0.5, na.rm=TRUE),
            distance_V_q25 = quantile(distance_V, 0.75, na.rm=TRUE),
            distance_V_q75 = quantile(distance_V, 0.25, na.rm=TRUE),
            prediction_V_mean= mean(distance_V),
            prediction_V_q50 = quantile(distance_V, 0.5, na.rm=TRUE),
            prediction_V_q25 = quantile(distance_V, 0.75, na.rm=TRUE),
            prediction_V_q75 = quantile(distance_V, 0.25, na.rm=TRUE),
            TPR_q50 = quantile(TPR, 0.5, na.rm=TRUE),
            TPR_q25 = quantile(TPR, 0.75, na.rm=TRUE),
            TPR_q75 = quantile(TPR, 0.25, na.rm=TRUE),
            FPR_mean = mean(FPR, na.rm=TRUE),
            FPR_q50 = quantile(FPR, 0.5, na.rm=TRUE),
            FPR_q25 = quantile(FPR, 0.75, na.rm=TRUE),
            FPR_q75 = quantile(FPR, 0.25, na.rm=TRUE),
            FNR_mean = mean(FNR, na.rm=TRUE),
            FNR_q50 = quantile(FNR, 0.5, na.rm=TRUE),
            FNR_q25 = quantile(FNR, 0.75, na.rm=TRUE),
            FNR_q75 = quantile(FNR, 0.25, na.rm=TRUE),
            time_med = quantile(time, 0.5, na.rm=TRUE),
            time_mean = mean(time),
            counts = n()

  ) %>%
  ungroup()
write_csv(summ, "~/Downloads/results_sparse_rrr_experiments_fall.csv")


summ <- read_csv( "~/Downloads/results_sparse_rrr_experiments_fall.csv")

unique(summ$method)



legend_order <- c("Oracle", 
                  "FIT_SAR_CV",
                  #"FIT_SAR_BIC", 
                  "Witten_Perm", 
                  "Witten.CV",
                  "SCCA_Parkhomenko",
                  "Waaijenborg-CV", 
                  "ccar3",
                 # "cca_rrr_-ADMM-one-iteration",
                  "Fantope",
                  "SGCA",
                  "Chao" )
my_colors  <- c(           "#999999",
                "#009E73",
                #"#6EE212",
                 "#0072B2",
                 "#56B4E9",
                "#F0E442",
                "#6EE212",
                "#000000",
                "#D55E00",
                "#CC79A7",
                "#E69F00"
                )

unique(summ$method)
labels_n <-    c("Oracle", 
                 "SAR CV (Wilms et al)",
                 #"SAR BIC (Wilms et al)",
                 "Sparse CCA, permuted\n(Witten et al)",
                 "Sparse CCA with CV\n(Witten et al)",
                 "SCCA (Parkhomenko et al)",
                 "Sparse CCA with CV\n(Waaijenborg et al)",
                 #"Sparse CCA (Waaijenborg et al)",
                 "CCAR3 (this paper)",
                 "Fantope-CCA  (Gao et al)",
                 "SGCA (Gao and Ma)",
                 "Sparse CCA (Gao et al)")
theme_set(theme_bw(base_size = 14))


unique(summ$method)
colnames(results)
colnames(summ)
unique(summ$nnzeros)
unique(summ$r)
unique(summ$r_pca)
unique(summ$p1)
unique(summ$p2)
unique(summ$overlapping_amount)

summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))

unique(summ$p2)
ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        ((p2==5) & p1 < 2000)| p1<3000,
                        p2!=10,
                        n==500,
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=p1, 
    y = distance_tot_mean, 
    colour =method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  #geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
  #                  colour =method), width=0.05, alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(p2 ~ theta_strength, scales = "free",labeller = as_labeller(c(`5` = "q = 5",
                                                                          `10` = "q = 10",
                                                                          `20` = "q = 20",
                                                                         `30` = "q = 30",
                                                                         `50` = "q = 50",
                                                                         `80` = "q = 80",
                                                                         `100` = "n = 100",
                                                                         `200` = "n = 200",
                                                                         `300` = "n = 300",
                                                                         `500` = "n = 500",
                                                                         `high` = "High",
                                                                         `1000` = "n = 1,000",
                                                                         `2000` = "n = 2,000",
                                                                         `10000` = "n = 10,000",
                                                                         `medium` = "Medium",
                                                                         `low` = "Low"
                                                                         
  ))) +
  xlab("p") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  #scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top")


ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        #((p2==5) & p1 < 2000)| p1<3000,
                        p2 %in% c(50, 80),
                        p1==1000,
                        n< 25000,
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=n, 
    y = distance_tot_mean, 
    colour =method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  #geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
  #                  colour =method), width=0.05, alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`5` = "q = 5",
                                                                          `10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `80` = "q = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("n") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  #scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top")
uniq
unique(summ$r)
summ$time_mean
ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        ((p2==5) & p1 < 2000)| p1<3000,
                        p2!=10,
                        n==500,
                        method %in% c( 
                                      "FIT_SAR_CV",
                                      #"FIT_SAR_BIC", 
                                      "Witten_Perm", 
                                      "Witten.CV",
                                      "ccar3",
                                      # "cca_rrr_-ADMM-one-iteration",
                                      "Fantope",
                                      "SGCA",
                                      "Chao" ),
                        method!="Oracle",
                        overlapping_amount==1
),
aes(x=p1, 
    y = time_mean, 
    colour =method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  #geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
  #                  colour =method), width=0.05, alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`5` = "q = 5",
                                                                          `10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `80` = "q = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Time Elapsed (s)")) +
  labs(colour="Method") + 
  scale_y_log10()+
  #scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top")
unique(summ$n)

ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        p2 > 10,
                        p1 == 1000,
                        theta_strength == "medium",
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=n, 
    y = distance_tot_mean, 
    colour =method)) +
  geom_point(size=2, position = position_dodge(width = 0.2),)+
  geom_line(linewidth=0.8, position = position_dodge(width = 0.2),)+
  geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
                    colour =method), position = position_dodge(width = 0.2), size = 0.6)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`5` = "q = 5",
                                                                          `10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `80` = "q = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("n") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top")



ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        n==500,
                        p2==80,
                        theta_strength == "medium",
                        p1<3000,
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=p1, 
    y = time_mean, 
    colour =method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
                    colour =method), width=0.05, alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`5` = "q = 5",
                                                                          `10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `80` = "q = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Time (s)")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+
  theme(legend.position = "top")

test = summ %>% filter( r_pca == 5, r==3,
                 nnzeros==10, 
                 p1 == 1000,
                 n==500,
)

ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==10, 
                        n==500,
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=p1, 
    y = distance_U_q50, 
    colour =method)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=distance_U_q25, ymax=distance_U_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `80` = "q = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Distance U")) +
  labs(colour="Method") + 
  #scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(summ %>% filter( r_pca == 5, r==3,
                        nnzeros==10, 
                        method %in% legend_order,
                        overlapping_amount==1
),
aes(x=n, 
    y = distance_tot_q50, 
    colour =method)) +
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(p2~ theta_strength, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `70` = "q = 70",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("n") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==5, 
                        method %in% legend_order
),
aes(x=p1, 
    y = time_mean, 
    colour =method)) +
  geom_point()+
  geom_line()+
  #geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
  #                  colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `70` = "q = 70",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Time (in seconds)")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==5, 
                        method %in% legend_order
),
aes(x=p1, 
    y = FNR_mean, 
    colour =method)) +
  geom_point(size=1.6)+
  geom_line()+
  #geom_jitter(width = 0.2, height = 0) + # Jitter only in the x-direction
  geom_errorbar(aes(ymin=FNR_q25, ymax=FNR_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                          `30` = "q = 30",
                                                                          `50` = "q = 50",
                                                                          `70` = "q = 70",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("False Negative Rate")) +
  labs(colour="Method") + 
  #scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


r_pca = 0
ggplot(summ %>% filter( r_pca == r_pca, 
                        nnzeros==5,
                        str_detect(method, "group-RRR"), !str_detect(method, "opt")),
       aes(x=lambda, 
       y = distance_tot_q50, 
       colour ="group")) +
  geom_point()+
  geom_line() +
  geom_point(data = summ %>% filter( r_pca == r_pca, 
                                      nnzeros==5,
                                      str_detect(method, "RRR"), 
                                      !str_detect(method, "opt"),
                                      !str_detect(method, "CVX")),
                     aes(x=lambda, 
                         y = distance_tot_q50, 
                         colour ="RRR"))+
  facet_grid(theta_strength~ p1 + n + nnzeros, scales = "free") +
  xlab("lambda") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

  geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) 


