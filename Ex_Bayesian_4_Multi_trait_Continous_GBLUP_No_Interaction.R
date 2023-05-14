# Import SKM library
library(SKM)

# Load the dataset
load("Chis_toy.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
geno_order <- sort(rownames(Geno))
# Order Geno by Line
Geno <- Geno[geno_order, geno_order]
Geno=as.matrix(Geno)

# Design matrix of lines (genotypes)
Z_g<- model.matrix(~ 0 + Line, data = Pheno)
K_g=Z_g%*%Geno%*%t(Z_g) ###Linear kernel

# Env design matrix without the first column
X_E<- model.matrix(~ 0 + Env, data = Pheno)[, -1]
K_e=X_E%*%t(X_E)/ncol(X_E)

###########Genotype by Environment interaction
K_ge <- K_e* K_g

# Put the matrices in the expected format with the model
# to use with each one
ETA <- list(
  list(x =K_e, model = "BGBLUP"),
  list(x =K_g, model = "BGBLUP")
)

# Retrieve the two continuos response variables
y <- Pheno[,c("GY", "Height")]

# Folds generation for random cross validation
#folds <- SKM::cv_random(
#  records_number = nrow(Pheno),
#  folds_number = 5,
#  testing_proportion = 0.25)
folds <- cv_random_line(lines =as.character(Pheno$Line),folds_number=5, testing_proportion =0.2)
#

# List for storing the individual predictions at each fold
Predictions <- list(GY = data.frame(), Height = data.frame())

for (i in seq_along(folds)) {
  fold <- folds[[i]]
  
  # Set testing indices to NA in the response variables
  y_na <- y
  y_na[fold$testing,] <- NA
  y_testing <- y[fold$testing,]
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model<-SKM::bayesian_model(
    ETA,
    y_na,
    iterations_number = 1000,
    burn_in = 500,
    seed =323
  )
  
  # Predict over testing set
  predictions <- predict(model)
  
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

write.csv(GY_Summaries$env,file="Summary_Chis_GBLUP_Multitrait_GY_NO_GE.csv")

## Errors metrics for prediction performance
Height_Summaries<-SKM::gs_summaries(Predictions$Height)

###Print summaries by line, environment and fold
#Height_Summaries$line
Height_Summaries$env
#Height_Summaries$fold

write.csv(Height_Summaries$env,file="Summary_Chis_GBLUP_Multitrait_Height_NO_GE.csv")
