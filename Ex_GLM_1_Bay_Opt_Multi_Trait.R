# Import SKM library
library(SKM)

# Load the dataset
load("Chis_toy.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
geno_order <- sort(rownames(Geno))
#########Ordering markers
Markers=Markers
dim(Markers)
Markers_sort<-Markers[geno_order,]
head(Markers_sort[,1:5])
Markers_S=scale(Markers_sort[,-1])
# Design matrix of genotypes
Z_g<- model.matrix(~ 0 + Line, data = Pheno)
#K_g=Z_g%*% Geno%*%t(Z_g) ###Linear kernel
X_g=Z_g%*%Markers_S
dim(X_g)
# Env design matrix without the first column
X_E<- model.matrix(~ 0 + Env, data = Pheno)[, -1]
#K_e=X_E%*%t(X_E)/ncol(X_E)
#K_ge <- K_e* K_g

X_ge<- model.matrix(~0+X_g:X_E)
head(X_ge)
# Bind all matrices in a single one
X<-cbind(X_E,X_g,X_ge)

# Retrieve the two continuos response variables
head(Pheno)
# Retrieve the two continuos response variables
y <- Pheno[,c("GY", "Height")]

folds <- cv_random_line(lines =as.character(Pheno$Line),folds_number=5, testing_proportion =0.2)
#

# Folds generation for random cross validation
#folds <- SKM::cv_random(
#  records_number = nrow(Pheno),
#  folds_number = 5,
#  testing_proportion = 0.25)

# List for storing the individual predictions at each fold
Predictions <- list(GY = data.frame(), Height = data.frame())


for (i in seq_along(folds)) {
 fold <- folds[[i]]
  
  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training,]
  y_testing <- y[fold$testing,]
  
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <-SKM::generalized_linear_model(
    # Predictor variables
    X_training,
    # Response variable
    y_training,
    
    # Tunable hyperparameters
    alpha = list(min = 0, max = 1),
    
    # Tune configuration parameters
    tune_type = "bayesian_optimization",
    tune_bayes_samples_number = 5,
    tune_bayes_iterations_number = 10,
    
    # Other algorithm's parameters
    lambdas_number = 100,
    standardize = TRUE,
    
    # Seed for reproducible results
    seed =323,
    verbose = TRUE
  )
  
  # Predict over testing set
  predictions <- predict(model,X_testing )
  
  # Save the predictions along with environment and line information
  Predictions$GY <- rbind(
    Predictions$GY,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing$GY,
      Predicted = predictions$GY$predicted
    )
  )
  
  # Save the predictions along with environment and line information
  Predictions$Height <- rbind(
    Predictions$Height,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing$Height,
      Predicted = predictions$Height$predicted
    )
  )
}
## Errors metrics for prediction performance
GY_Summaries<-SKM::gs_summaries(Predictions$GY)

###Print summaries by line, environment and fold
#GY_Summaries$line
GY_Summaries$env
#GY_Summaries$fold

write.csv(GY_Summaries$env,file="Summary_Chis_GLM_Bay_Opt_Multitrait_GY.csv")

## Errors metrics for prediction performance
Height_Summaries<-SKM::gs_summaries(Predictions$Height)

###Print summaries by line, environment and fold
#Height_Summaries$line
Height_Summaries$env
#Height_Summaries$fold

write.csv(Height_Summaries$env,file="Summary_Chis_GLM_Bay_Opt_Multitrait_Height.csv")
