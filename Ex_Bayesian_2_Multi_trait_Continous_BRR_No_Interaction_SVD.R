# Import SKM library
library(SKM)

# Load the dataset
load("Chis_toy.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
geno_order <- sort(rownames(Geno))
# Ordering Markers by Line
#Geno <- Geno[geno_order, geno_order]
#Geno=as.matrix(Geno)
Markers=Markers
dim(Markers)
Markers_sort<-Markers[geno_order,]
head(Markers_sort[,1:5])
Markers_SA=scale(Markers_sort[,-1])
svd_Markers=svd(Markers_SA)
Markers_S=svd_Markers$u%*%diag(svd_Markers$d)

# Design matrix of lines
Z_g<- model.matrix(~ 0 + Line, data = Pheno)
#K_g=Z_g%*% Geno%*%t(Z_g) ###Linear kernel
X_g=Z_g%*%Markers_S
dim(X_g)
# Env design matrix without the first column
X_E<- model.matrix(~ 0 + Env, data = Pheno)[, -1]
#K_e=X_E%*%t(X_E)/ncol(X_E)
#K_ge <- K_e* K_g

X_ge<- model.matrix(~0+X_g:X_E)
#head(X_ge)
# Put the matrices in the expected format with the model
# to use with each one
ETA <- list(
  list(x =X_E, model = "FIXED"), ###Effects of environments
  list(x =X_g, model = "BRR") #####Effects of lines (Genotypes)

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
    iterations_number = 100,
    burn_in = 50,
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
#DTHDSummaries$line
Summary_Env_GY=GY_Summaries$env
#DTHDSummaries$fold
Summary_Env_GY
write.csv(Summary_Env_GY,file="Summary_Chis_BRR_Multitrait_GY_NO_GE.csv")

## Errors metrics for prediction performance
Height_Summaries<-SKM::gs_summaries(Predictions$Height)

###Print summaries by line, environment and fold
#DTHDSummaries$line
Summary_Env_Height=Height_Summaries$env
#DTHDSummaries$fold
Summary_Env_Height
write.csv(Summary_Env_Height,file="Summary_Chis_BRR_Multitrait_Height_NO_GE.csv")
