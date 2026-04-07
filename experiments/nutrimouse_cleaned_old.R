library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("~/Documents/CCAR3_code")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/evaluation.R')
source("experiments/simulations/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/alternative_methods/scca_chao.R')

source("~/Documents/ccar3/R/helpers.r")
source("~/Documents/ccar3/R/globals.R")
source("~/Documents/ccar3/R/utils.R")
source("~/Documents/ccar3/R/graph_utils.R")
source("~/Documents/ccar3/R/reduced_rank_regression.R")
source("~/Documents/ccar3/R/graph_reduced_rank_regression.R")



#Load X and Y from nutrimouse data
data(nutrimouse)
X = as.matrix(nutrimouse$gene %>% mutate_all(~rank(.)))
Y = as.matrix(nutrimouse$lipid %>% mutate_all(~rank(.)))
X = scale(X, center = T, scale = T)
Y = scale(Y, center = T, scale = T)
n = nrow(X)
p = ncol(X)
q = ncol(Y)

props = seq(0, 0.3, 0.05)
seeds = 1:30
lambda1s = c(0.001, 0.01, 0.1, 1) 
lambda2s = c(0.001, 0.01, 0.1, 1)
set.seed(0)
nfold = 10

r = 5
#### Apply CCA
X = scale(X)
Y =scale(Y)
lambda_Kx = NULL
param_lambda=c(10^seq(-3, 3, length.out = 100))
kfolds=5
penalty="l21"
solver="admm"
LW_Sy = TRUE

### Let's do a cross validation setting

p = ncol(X)
Y = scale(Y)
nb_experiments = 10
set.seed(12345)
correlation <-c()
for (exp in 1:nb_experiments){
  folds = createFolds(1:nrow(X), 8)
  order = 1:length(folds)
  for (i in  1:length(folds)){
    index = order[ifelse(i < length(folds), i + 1, (i+1)%%length(folds))]
    index2 =order[ifelse(i < length(folds) -1, i + 2, (i+2)%%length(folds))]
    print(c(i, index, index2))
    for (lambda in 10^seq(from=-3, 0, length.out=20)){
      final = cca_rrr(as.matrix(X)[-c(folds[[index]], folds[[index2]]),], 
                      as.matrix(Y)[-c(folds[[index]],
                                      folds[[index2]]),],  
                      lambda = lambda, 
                      r=r,
                      rho=1, niter=2 * 1e4,
                      preprocess=FALSE, 
                      thresh=1e-6,
                      solver= "ADMM",
                      LW_Sy = LW_Sy)
      
      correlation <- rbind(
        correlation,
        c("CCA_rrr",
          lambda,
          exp,
          i,
          diag(cov(as.matrix(X)[-c(folds[[index]],
                                   folds[[index2]]),] %*% final$U,
                   as.matrix(Y)[-c(folds[[index]],
                                   folds[[index2]]),] %*%  final$V)),
          apply(((as.matrix(X)[-c(folds[[index]],
                                  folds[[index2]]),] %*% final$U) -
                   (as.matrix(Y)[-c(folds[[index]],
                                    folds[[index2]]),] %*%  final$V))^2, 2, mean),
          diag(t(as.matrix(X)[folds[[index]],] %*% final$U) %*%
                 (as.matrix(Y)[folds[[index]],] %*%  final$V)),
          diag(cor(as.matrix(X)[folds[[index]],] %*% final$U, (as.matrix(Y)[folds[[index]],] %*%  final$V))),
          apply(((as.matrix(X)[folds[[index]],] %*% final$U) -
                   (as.matrix(Y)[folds[[index]],] %*%  final$V))^2,2,mean),
          subdistance(as.matrix(X)[folds[[index]],] %*% final$U,
                      as.matrix(Y)[folds[[index]],] %*%  final$V),
          diag(t(as.matrix(X)[folds[[index2]],] %*% final$U) %*%
                 (as.matrix(Y)[folds[[index2]],] %*%  final$V)),
          diag(cor(as.matrix(X)[folds[[index2]],] %*% final$U, (as.matrix(Y)[folds[[index2]],] %*%  final$V))),
          apply(((as.matrix(X)[folds[[index2]],] %*% final$U) -
                   (as.matrix(Y)[folds[[index2]],] %*%  final$V))^2,2,mean),
          subdistance(as.matrix(X)[folds[[index2]],] %*% final$U, 
                      as.matrix(Y)[folds[[index2]],] %*%  final$V)
        ))
      
    }
    correlation_df = data.frame(correlation)
    colnames(correlation_df) = c("method", "lambda", "exp", "fold",
                                 "train_cov1",  "train_cov2",    "train_cov3",  "train_cov4",   "train_cov5", 
                                 "train_mse1",  "train_mse2",   "train_mse3",  "train_mse4", "train_mse5", 
                                 "test_prod1",  "test_prod2",   "test_prod3",  "test_prod4", "test_prod5", 
                                 "test_cor1",  "test_cor2",          "test_cor3",  "test_cor4", "test_cor5", 
                                 "test_mse1", "test_mse2",  "test_mse3", "test_mse4",  "test_mse5", "test_dist",
                                 "val_prod1",  "val_prod2",     "val_prod3",  "val_prod4",    "val_prod5",  
                                 "val_cor1",  "val_cor2",    "val_cor3",  "val_cor4",    "val_cor5",  
                                 "val_mse1",  "val_mse2", "val_mse3",  "val_mse4","val_mse5","val_dist")
    write_csv(correlation_df, "~/Downloads/nutrimouse_new_trial.csv")
    for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                     "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                     "SCCA_Parkhomenko", "Fantope", "SGCA")){
      
      print(paste0("Starting ", method))
      tryCatch({
        test1<-additional_checks(as.matrix(X)[-c(folds[[index2]]),], as.matrix(Y)[-c(folds[[index2]]),],
                                 S=NULL, 
                                 rank=r, kfolds=8, 
                                 method.type = method,
                                 lambdax= 10^seq(-3,1, length.out = 30),
                                 lambday =  10^seq(-3,1, length.out = 30))
        correlation <- rbind(
          correlation,
          c(method,
            lambda,
            exp,
            i,
            diag(cov(as.matrix(X)[-c(folds[[index]],
                                     folds[[index2]]),] %*% test1$u,
                     as.matrix(Y)[-c(folds[[index]],
                                     folds[[index2]]),] %*%  test1$v)),
            apply(((as.matrix(X)[-c(folds[[index]],
                                    folds[[index2]]),] %*% test1$u) -
                     (as.matrix(Y)[-c(folds[[index]],
                                      folds[[index2]]),] %*%  test1$v))^2, 2, mean),
            diag(t(as.matrix(X)[folds[[index]],] %*% test1$u) %*%
                   (as.matrix(Y)[folds[[index]],] %*%  test1$v)),
            diag(cor(as.matrix(X)[folds[[index]],] %*% test1$u, (as.matrix(Y)[folds[[index]],] %*%  test1$v))),
            apply(((as.matrix(X)[folds[[index]],] %*% test1$u) -
                     (as.matrix(Y)[folds[[index]],] %*%  test1$v))^2,2,mean),
            subdistance(as.matrix(X)[folds[[index]],] %*% test1$u, 
                        as.matrix(Y)[folds[[index]],] %*%  test1$v),
            diag(t(as.matrix(X)[folds[[index2]],] %*% test1$u) %*%
                   (as.matrix(Y)[folds[[index2]],] %*%  test1$v)),
            diag(cor(as.matrix(X)[folds[[index2]],] %*% test1$u, (as.matrix(Y)[folds[[index2]],] %*%  test1$v))),
            apply(((as.matrix(X)[folds[[index2]],] %*% test1$u) -
                     (as.matrix(Y)[folds[[index2]],] %*%  test1$v))^2,2,mean),
            subdistance(as.matrix(X)[folds[[index2]],] %*% test1$u, 
                        as.matrix(Y)[folds[[index2]],] %*%  test1$v)
          ))
      }, error = function(e) {
        # Print the error message
        cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
        # Skip to the next iteration
      })
      correlation_df = data.frame(correlation)
      colnames(correlation_df) = c("method", "lambda", "exp", "fold",
                                   "train_cov1",  "train_cov2",    "train_cov3",  "train_cov4",   "train_cov5", 
                                   "train_mse1",  "train_mse2",   "train_mse3",  "train_mse4", "train_mse5", 
                                   "test_prod1",  "test_prod2",   "test_prod3",  "test_prod4", "test_prod5", 
                                   "test_cor1",  "test_cor2",          "test_cor3",  "test_cor4", "test_cor5", 
                                   "test_mse1", "test_mse2",  "test_mse3", "test_mse4",  "test_mse5", "test_dist",
                                   "val_prod1",  "val_prod2",     "val_prod3",  "val_prod4",    "val_prod5",  
                                   "val_cor1",  "val_cor2",    "val_cor3",  "val_cor4",    "val_cor5",  
                                   "val_mse1",  "val_mse2", "val_mse3",  "val_mse4","val_mse5","val_dist")
      write_csv(correlation_df, "~/Downloads/nutrimouse_new_trial3.csv")
    }
    
  }
  
}
STOP
correlation_df =  read_csv( "~/Downloads/nutrimouse_new_trial2.csv")# data.frame(correlation)
colnames(correlation_df) = c("method", "lambda", "exp", "fold",
                             "train_cov1",  "train_cov2",    "train_cov3",  "train_cov4",   "train_cov5", 
                             "train_mse1",  "train_mse2",   "train_mse3",  "train_mse4", "train_mse5", 
                             "test_prod1",  "test_prod2",   "test_prod3",  "test_prod4", "test_prod5", 
                             "test_cor1",  "test_cor2",          "test_cor3",  "test_cor4", "test_cor5", 
                             "test_mse1", "test_mse2",  "test_mse3", "test_mse4",  "test_mse5", "test_dist",
                             "val_prod1",  "val_prod2",     "val_prod3",  "val_prod4",    "val_prod5",  
                             "val_cor1",  "val_cor2",    "val_cor3",  "val_cor4",    "val_cor5",  
                             "val_mse1",  "val_mse2", "val_mse3",  "val_mse4","val_mse5","val_dist")


summ1  = correlation_df %>%
  filter(((lambda<1)& (method == "CCA_rrr")) | (method!="CCA_rrr")) %>%
  mutate(test_mse = (test_mse1 + test_mse2 + test_mse3 + test_mse4 +test_mse5 )/5,
         train_mse = train_mse1 + train_mse2 +train_mse3 + train_mse4+train_mse5,
         val_mse =  (val_mse1 + val_mse2 + val_mse3 + val_mse4+val_mse5 )/5,
         test_prod = (test_prod1 + test_prod2 +  test_prod3+ test_prod4+ test_prod5)/5,
         test_cor = (test_cor1 + test_cor2  + test_cor3  + test_cor4  + test_cor5)/5,
         train_cov = (train_cov1 + train_cov2 + train_cov3+train_cov4+train_cov5)/5,
         val_prod = (val_prod1 + val_prod2 + val_prod3+ val_prod4 + val_prod5)/5,
         val_cor = (val_cor1 + val_cor2 + val_cor3 + val_cor4 + val_cor5)/5
  ) %>%
  group_by(method, lambda, exp) %>%
  summarise_all(mean)  %>%
  drop_na(val_cor) %>%
  slice_min(test_mse, n=1) %>%
  ungroup() 


#drop_na()

ggplot(summ1,
       aes(x=method, y=val_mse )) +
  geom_boxplot()

summary_correlation = summ1 %>%
  group_by(method) %>%
  summarise(counts =c(),
            test_cor= median(test_cor),
            test_mse = median(test_mse),
            test_dist = median(test_dist),
            val_cor= mean(val_cor),
            val_mse = mean(val_mse),
            val_dist = mean(val_dist))%>%
  arrange((val_mse)) %>% ungroup()

correlation_df$exp

summary_correlation2 = summ1 %>%
  group_by(method) %>%
  summarise(counts =c(),
            test_cor= median(test_cor),
            test_mse = median(test_mse),
            test_dist = median(test_dist),
            val_cor= median(val_cor),
            val_mse = median(val_mse),
            val_dist = median(val_dist))%>%
  arrange((val_mse)) %>% ungroup()

summary_correlation$lambda = as.numeric(summary_correlation$lambda )
ggplot(summary_correlation %>% 
         filter(method == "CCA_graph_rrr")%>% ungroup())+
  geom_line(aes(x = as.numeric(lambda), y=test_mse))+
  geom_point(aes(x = as.numeric(lambda), y=test_mse))+
  scale_x_log10()


lambda_opt = summary_correlation$lambda[which.min(summary_correlation$train_mse[which(summary_correlation$method ==  "CCA_graph_rrr")])]
lambda_opt = 0.009236709

relevant_correlations = summary_correlation %>% 
  group_by(method) %>%
  dplyr::select(method, test_mse, test_cor, val_mse, val_cor)



relevant_correlations = summary_correlation %>% 
  filter( (method == "CCA_graph_rrr" & lambda  < 0.12 & lambda >0.09 ) | 
            method!= "CCA_graph_rrr" ) %>%
  dplyr::select(method, test_mse, test_cor, val_mse, val_cor)


library(knitr)
library(kableExtra)
relevant_correlations
latex_table <- kable(relevant_correlations, format = "latex", booktabs = TRUE)


r = 5


test1<-additional_checks(as.matrix(X), as.matrix(Y),  
                         S=NULL, 
                         rank=r, kfolds=5, 
                         method.type =  "Fantope",
                         lambdax= 10^seq(-3,1, length.out = 30),
                         lambday = c(0))

df = data.frame(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u) #)
#df = data.frame(1/sqrt(nrow(X)) * as.matrix(X)%*% final$U) #)
#df = data.frame(1/sqrt(nrow(X)) * as.matrix(Y)%*% test1$v) #)
gen =  nutrimouse$genotype
diet = nutrimouse$diet
df["diet"] = diet
df["gen"] = gen



set.seed(107)
inTrain <- createDataPartition(
  y = df$diet,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
training <- df[ inTrain,]
testing  <- df[-inTrain,]
svm_model <- train(diet ~ .-gen, data = training, method = "lda",
                   trControl = trainControl(method = "cv", number = 5),
                   preProcess = c("center", "scale"),
                   tuneLength = 5)
# Evaluate the model
predictions <- predict(svm_model, newdata = testing)
confusionMatrix(predictions, testing$diet)

test_cluster = kmeans(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u[,1:5], 2) #)
library(mclust)
model <- Mclust(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u, G=5:5)
confusionMatrix(factor(model$classification), factor(as.numeric(df$diet)))



# inTrain2 <- createDataPartition(
#   y = df$gen,
#   ## the outcome data are needed
#   p = .75,
#   ## The percentage of data in the
#   ## training set
#   list = FALSE
# )
training2 <- df[ inTrain2,]
testing2  <- df[-inTrain2,]
svm_model <- train(gen ~ .-diet, data = training2, 
                   method = "lda",
                   trControl = trainControl(method = "cv", number = 5),
                   preProcess = c("center", "scale"),
                   tuneLength = 5)
predictions <- predict(svm_model, newdata = testing2)
confusionMatrix(predictions, testing2$gen)


test_cluster2 = kmeans(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u[,1:5], 2) #)
model <- Mclust(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u, G=2:2)
confusionMatrix(factor(model$classification), factor(as.numeric(df$gen)))



library(ellipse)
legend_order <- c("lin", "sun", "fish",
                  "ref", "coc")
my_colors <- c(  "red","orange",  "dodgerblue", 
                 "black", "brown"
                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
                 #"navyblue", #"cyan", 
                 #"dodgerblue"
)
test_cluster = kmeans(1/sqrt(nrow(X)) * as.matrix(X)%*% test1$u, 5) #)
confusionMatrix(factor(test_cluster$cluster), factor(as.numeric(df$diet)))
labels_n <-   c("LIN", "SUN", "FISH",
                "REF", "COC")

ellipse.level =0.95
theme_set(theme_bw(base_size = 18))
ggplot(df, aes(x=X1, y=X2 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -0.4, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="lin")], 
                                                df$X2[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X1[which(diet=="lin")], 
                                                            df$X2[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="ref")], 
                                                df$X2[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="ref")], 
                          df$X2[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="sun")], 
                                                df$X2[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="sun")], 
                          df$X2[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="fish")], 
                                                df$X2[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="fish")], 
                          df$X2[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="coc")], 
                                                df$X2[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="coc")], 
                          df$X2[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-1") + 
  ylab("CD-2")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 2)))


ggplot(df, aes(x=X3, y=X4 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -0.4, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="lin")], 
                                                df$X4[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X3[which(diet=="lin")], 
                                                            df$X4[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="ref")], 
                                                df$X4[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="ref")], 
                          df$X4[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="sun")], 
                                                df$X4[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="sun")], 
                          df$X4[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="fish")], 
                                                df$X4[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="fish")], 
                          df$X4[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="coc")], 
                                                df$X4[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="coc")], 
                          df$X4[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-3") + 
  ylab("CD-4")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 2)))


ggplot(df, aes(x=X1, y=X5 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -0.4, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="lin")], 
                                                df$X5[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X1[which(diet=="lin")], 
                                                            df$X5[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="ref")], 
                                                df$X5[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="ref")], 
                          df$X5[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="sun")], 
                                                df$X5[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="sun")], 
                          df$X5[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="fish")], 
                                                df$X5[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="fish")], 
                          df$X5[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="coc")], 
                                                df$X5[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="coc")], 
                          df$X5[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-1") + 
  ylab("CD-5")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 2)))


Vhat_comp = test1$v
rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)
df_V = data.frame(Vhat_comp)
colnames(df_V) = c("CD-1", "CD-2", "CD-3","CD-4", "CD-5")
df_V["Hepatic Fatty Acids"] = colnames(Y)
df_V = pivot_longer(df_V, cols=-c("Hepatic Fatty Acids"))



# Example new labels

# Assuming your questions are in a column named 'question'
theme_set(theme_bw(base_size = 18))
ggplot(df_V, aes(x = value, y = reorder(`Hepatic Fatty Acids`, value), fill=name)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ name, scale="free_y", nrow=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Hepatic Fatty Acids", x = "Response Intensity") + 
  labs(fill = "Canonical\nDirection")


ggplot(df_V %>% filter(abs(value) > 1e-1) %>% mutate(row = ifelse(name %in% c("CD-4", "CD-5"), 2,1)), aes(x = value, y = reorder(`Hepatic Fatty Acids`, value), fill=name)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(.~name, scale="free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Hepatic Fatty Acids", x = "Loading Value") + 
  labs(fill = "Canonical\nDirection")

Uhat_comp = test1$u
rownames(Uhat_comp) = colnames(X)
print(Uhat_comp)
heatmap(Uhat_comp)
df_U = data.frame(Uhat_comp)
colnames(df_U) = c("CD-1", "CD-2", "CD-3","CD-4", "CD-5")
df_U["Gene"] = colnames(X)
df_U = pivot_longer(df_U, cols=-c("Gene"))



# Example new labels
# Reorder the levels of the x-axis variable based on the ordering column
#df_U$Gene <- factor(df_U$Gene, levels = levels(df$Gene)[order(df_U$name)])
# Assuming your questions are in a column named 'question'
theme_set(theme_bw(base_size = 18))
ggplot(df_U %>% filter(abs(value) > 1e-1), aes(x = value, y = reorder(`Gene`, value), 
                                               fill=name)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~name, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Gene", x = "Response Intensity") + 
  labs(fill = "Canonical\nDirection")


r = 5
LW_Sy = TRUE
final = CCA_rrr(as.matrix(X),
                as.matrix(Y),  
                lambda = 0.1, 
                Kx=NULL, r=r,
                rho=1, niter=2 * 1e4,
                do.scale = FALSE, lambda_Kx=0,
                thresh=1e-6,
                solver= "ADMM",
                LW_Sy = LW_Sy)

df = data.frame(1/sqrt(nrow(X)) * as.matrix(X)%*% final$U) #)
#df = data.frame(1/sqrt(nrow(X)) * as.matrix(X)%*% final$U) #)
#df = data.frame(1/sqrt(nrow(X)) * as.matrix(Y)%*% test1$v) #)
gen =  nutrimouse$genotype
diet = nutrimouse$diet
df["diet"] = diet
df["gen"] = gen


set.seed(107)
inTrain <- createDataPartition(
  y = df$diet,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
#training <- df[ inTrain,]
#testing  <- df[-inTrain,]
svm_model <- train(diet ~ .-gen, data = training, method = "lda",
                   trControl = trainControl(method = "cv", number = 5),
                   preProcess = c("center", "scale"),
                   tuneLength = 5)
# Evaluate the model
predictions <- predict(svm_model, newdata = testing)
confusionMatrix(predictions, testing$diet)

library(mclust)
model <- Mclust(1/sqrt(nrow(X)) * as.matrix(X)%*% final$U[, 1:2], G=5:5)
confusionMatrix(factor(model$classification), factor(as.numeric(df$diet)))



inTrain2 <- createDataPartition(
  y = df$gen,
  ## the outcome data are needed
  p = .80,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
training2 <- df[ inTrain2,]
testing2  <- df[-inTrain2,]
svm_model <- train(gen ~ .-diet, data = training2, 
                   method = "lda",
                   trControl = trainControl(method = "cv", number = 5),
                   #preProcess = c("center", "scale"),
                   tuneLength = 5)
predictions <- predict(svm_model, newdata = testing2)
confusionMatrix(predictions, testing2$gen)
model <- Mclust(1/sqrt(nrow(X)) * as.matrix(X)%*% final$U, G=2:2)
confusionMatrix(factor(model$classification), factor(as.numeric(df$gen)))

library(ellipse)
legend_order <- c("lin", "sun", "fish",
                  "ref", "coc")
my_colors <- c(  "red","orange",  "dodgerblue", 
                 "black", "brown"
                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
                 #"navyblue", #"cyan", 
                 #"dodgerblue"
)

labels_n <-   c("LIN", "SUN", "FISH",
                "REF", "COC")

ellipse.level =0.95
theme_set(theme_bw(base_size = 18))
ggplot(df, aes(x=X1, y=X2 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -1, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="lin")], 
                                                df$X2[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X1[which(diet=="lin")], 
                                                            df$X2[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="ref")], 
                                                df$X2[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="ref")], 
                          df$X2[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="sun")], 
                                                df$X2[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="sun")], 
                          df$X2[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="fish")], 
                                                df$X2[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="fish")], 
                          df$X2[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="coc")], 
                                                df$X2[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="coc")], 
                          df$X2[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-1") + 
  ylab("CD-2")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 3)))


ggplot(df, aes(x=X3, y=X4 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -1, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="lin")], 
                                                df$X4[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X3[which(diet=="lin")], 
                                                            df$X4[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="ref")], 
                                                df$X4[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="ref")], 
                          df$X4[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="sun")], 
                                                df$X4[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="sun")], 
                          df$X4[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="fish")], 
                                                df$X4[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="fish")], 
                          df$X4[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X3[which(diet=="coc")], 
                                                df$X4[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X3[which(diet=="coc")], 
                          df$X4[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-3") + 
  ylab("CD-4")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 2)))


ggplot(df, aes(x=X1, y=X5 ,colour=diet))+
  geom_point(aes( shape=gen), size = 4)+
  geom_text(aes(label = gen), vjust = -1, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="lin")], 
                                                df$X5[which(diet=="lin")]))), 
                                    centre=colMeans(t(rbind(df$X1[which(diet=="lin")], 
                                                            df$X5[which(diet=="lin")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="lin")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="ref")], 
                                                df$X5[which(diet=="ref")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="ref")], 
                          df$X5[which(diet=="ref")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="ref")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="sun")], 
                                                df$X5[which(diet=="sun")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="sun")], 
                          df$X5[which(diet=="sun")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="sun")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="fish")], 
                                                df$X5[which(diet=="fish")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="fish")], 
                          df$X5[which(diet=="fish")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="fish")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(diet=="coc")], 
                                                df$X5[which(diet=="coc")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(diet=="coc")], 
                          df$X5[which(diet=="coc")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="coc")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid"
  )))) +
  xlab("CD-1") + 
  ylab("CD-5")+
  labs(colour = "Diet", shape = "Genotype")+
  guides(colour = guide_legend(override.aes = list(size = 3)))


Vhat_comp = final$V
rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)
df_V = data.frame(Vhat_comp)
colnames(df_V) = c("CD-1", "CD-2", "CD-3","CD-4", "CD-5")
df_V["Hepatic Fatty Acids"] = colnames(Y)
df_V = pivot_longer(df_V, cols=-c("Hepatic Fatty Acids"))



# Example new labels

# Assuming your questions are in a column named 'question'
theme_set(theme_bw(base_size = 18))
ggplot(df_V %>% filter(abs(value) > 1e-1) %>% mutate(row = ifelse(name %in% c("CD-4", "CD-5"), 2,1)), aes(x = value, y = reorder(`Hepatic Fatty Acids`, value), fill=name)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(.~name, scale="free_y", nrow=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Hepatic Fatty Acids", x = "Loading Value") + 
  labs(fill = "Canonical\nDirection")


Uhat_comp = final$U
rownames(Uhat_comp) = colnames(X)
print(Uhat_comp)
heatmap(Uhat_comp)
df_U = data.frame(Uhat_comp)
colnames(df_U) = c("CD-1", "CD-2", "CD-3","CD-4", "CD-5")
df_U["Gene"] = colnames(X)
df_U = pivot_longer(df_U, cols=-c("Gene"))



# Example new labels
# Reorder the levels of the x-axis variable based on the ordering column
#df_U$Gene <- factor(df_U$Gene, levels = levels(df$Gene)[order(df_U$name)])
# Assuming your questions are in a column named 'question'
theme_set(theme_bw(base_size = 18))
ggplot(df_U %>% filter(abs(value) > 1e-1), aes(x = value, y = reorder(`Gene`, value), 
                                               fill=name)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~name, scale="free_y", nrow=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Gene", x = "Loading Value") + 
  labs(fill = "Canonical\nDirection")
