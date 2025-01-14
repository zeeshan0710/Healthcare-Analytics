library(tidyverse)
library(ggplot2)
library(DT)
library(corrplot)
library(GGally)
library(ggforce)
library(Boruta)
library(caret)
library(randomForest)
library(reshape2)
library(pROC)
library(MASS)
library(Rtsne)
library(lubridate)
library(rgl)
library(tidyr)
library(cluster)
library(dendextend)
library(recipes)
library(smotefamily)


diabetic_data <- read.csv("diabetic_data.csv")
#ids_mapping <- read.csv("IDS_mapping.csv")
str(diabetic_data)

#data cleaning
diabetic_data[diabetic_data == "?"] <- NA
# Set missing values to NA according to mapping data
diabetic_data[(diabetic_data$admission_type_id %in% c(5,6,8)),"admission_type_id"] <- NA

diabetic_data[(diabetic_data$discharge_disposition_id %in% c(18,25,26)), "discharge_disposition_id"] <- NA

diabetic_data[(diabetic_data$admission_source_id %in% c(9,15,17,20,21)),"admission_source_id"] <- NA

#data type conversion
diabetic_data$race <- as.factor(diabetic_data$race)
diabetic_data$gender <- as.factor(diabetic_data$gender)
diabetic_data$age <- as.factor(diabetic_data$age)
diabetic_data$payer_code <- as.factor(diabetic_data$payer_code)
diabetic_data$max_glu_serum <- as.factor(diabetic_data$max_glu_serum)
diabetic_data$A1Cresult <- as.factor(diabetic_data$A1Cresult)
diabetic_data$admission_type_id <- as.factor(diabetic_data$admission_type_id)
diabetic_data$admission_source_id <- as.factor(diabetic_data$admission_source_id)
diabetic_data$discharge_disposition_id <- as.factor(diabetic_data$discharge_disposition_id)

for (i in 25:50){
  diabetic_data[,i] <- as.factor(diabetic_data[,i])
}

convert_diag <- function(feature,diabetic_data)
{
  val_name <- paste(feature,"_group",sep='')
  diabetic_data[,val_name] <- NA
  for (i in 1:length(diabetic_data[,feature]))
  {
    if (!is.na(diabetic_data[i,feature]))
    {
      if (grepl("[E-Ve-v]", diabetic_data[i,feature], perl = T)) {diabetic_data[i,val_name] <- "Other"}
      else{
        cur_val <- as.numeric(diabetic_data[i,feature])
        if (cur_val %in% c(390:459,785)) {diabetic_data[i,val_name] <- "Circulatory"}
        else if (cur_val %in% c(460:519,786)) {diabetic_data[i,val_name] <- "Respiratory"}
        else if (cur_val %in% c(520:579,787)) {diabetic_data[i,val_name] <- "Digestive"}
        else if (cur_val %in% c(800:999)) {diabetic_data[i,val_name] <- "Injury"}
        else if (cur_val %in% c(710:739)) {diabetic_data[i,val_name] <- "Musculoskeletal"}
        else if (cur_val %in% c(580:629,788)) {diabetic_data[i,val_name] <- "Genitourinary"}
        else if (cur_val %in% c(140:239)) {diabetic_data[i,val_name] <- "Neoplasms"}
        else if (cur_val>=250 & cur_val<251) {diabetic_data[i,val_name] <- "Diabetes"}
        else {diabetic_data[i,val_name] <- "Other"}
      }
    }
  }
  return(diabetic_data)
}

for (i in c("diag_1","diag_2","diag_3")){
  diabetic_data <- convert_diag(i,diabetic_data)
}
diabetic_data$diag_1_group <- as.factor(diabetic_data$diag_1_group)
diabetic_data$diag_2_group <- as.factor(diabetic_data$diag_2_group)
diabetic_data$diag_3_group <- as.factor(diabetic_data$diag_3_group)

diabetic_data <- subset(diabetic_data, select = -c(diag_1,diag_2,diag_3))

#data filtering
hospiceORdeathORnull_id <- c(11,13,14,19,20,21,18,25,26)
diabetic_data <- group_by(diabetic_data, patient_nbr) %>% slice(1) %>% 
  filter(!discharge_disposition_id %in% hospiceORdeathORnull_id)
diabetic_data <- as.data.frame(diabetic_data)

#feature pruning
miss_data <- data.frame(features = colnames(diabetic_data),
                        missing_rates = round(colSums(is.na(diabetic_data))/nrow(diabetic_data),3))
miss_data <- miss_data[order(-miss_data$missing_rates),]
# datatable(miss_data)
miss_data$features = factor(miss_data$features, levels = unique(miss_data$features))
ggplot(miss_data, aes(features,missing_rates,fill=missing_rates))+
  geom_bar(stat="identity")+
  ggtitle("Missing value rates of 50 features") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
diabetic_data <- subset(diabetic_data, select = -c(medical_specialty,encounter_id,weight,patient_nbr,payer_code))

#EDA
#dist of categorical variables
data_factor <- diabetic_data[, unlist(lapply(diabetic_data, is.factor))]
num_subplots <- 8

for (col_start in seq(1, ncol(data_factor),num_subplots)){
  long <- pivot_longer(data_factor[,col_start:min(col_start+num_subplots-1,ncol(data_factor))], everything(),
                       names_to = "features", 
                       values_to = "categories")
  
  g <- ggplot(long, aes(x=categories, fill = features))
  
  print(g + geom_histogram(stat = "count") +
          facet_wrap(~features,scales = "free", ncol = 4) +
          ggtitle("Histogram Plots") + coord_flip() + 
          scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
          theme(legend.position = "top", strip.text.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

#Distributions of Numerical Features
data_numeric <- diabetic_data[, unlist(lapply(diabetic_data, is.numeric))]
num_subplots <- 8
for (col_start in seq(1, ncol(data_numeric),num_subplots)){
  long <- pivot_longer(data_numeric[,col_start:min(col_start+num_subplots-1,ncol(data_numeric))], everything(),
                       names_to = "features", 
                       values_to = "values")
  
  g <- ggplot(long, aes(x=values, fill = features))
  
  print(g + geom_histogram(stat = "count") +
          facet_wrap(~features,scales = "free", ncol = 4) +
          ggtitle("Histogram Plots") + coord_flip() + 
          scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
          theme(legend.position = "top", strip.text.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

#Correlation between Numerical Variables
data_numeric_cor <- cor(data_numeric,use="pairwise.complete.obs")
corrplot(data_numeric_cor, method="color", type="upper", tl.cex = 0.8, mar=c(0,0,1,0), 
         tl.pos='d',tl.offset = 1)

#Data preprocessing
# Load necessary libraries
library(dplyr)
library(smotefamily)

#Preprocessing Data
#Replace "?" with NA
diabetic_data[] <- lapply(diabetic_data, function(x) ifelse(x == "?", NA, x))

#Handle Specific Values
diabetic_data$admission_type_id[diabetic_data$admission_type_id %in% c(5, 6, 8)] <- NA
diabetic_data$discharge_disposition_id[diabetic_data$discharge_disposition_id %in% c(18, 25, 26)] <- NA
diabetic_data$admission_source_id[diabetic_data$admission_source_id %in% c(9, 15, 17, 20, 21)] <- NA

#Impute Missing Values
#For numeric columns, replace NA with the median
numeric_cols <- sapply(diabetic_data, is.numeric)
diabetic_data[numeric_cols] <- lapply(diabetic_data[numeric_cols], function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

#For character columns, replace NA with "Unknown"
char_cols <- sapply(diabetic_data, is.character)
diabetic_data[char_cols] <- lapply(diabetic_data[char_cols], function(x) ifelse(is.na(x), "Unknown", x))

#Convert Non-Numeric Features to Numeric
#Convert character columns to factors and then to numeric
diabetic_data[char_cols] <- lapply(diabetic_data[char_cols], function(x) as.numeric(as.factor(x)))

#Verify No NA Values Remain
if (any(is.na(diabetic_data))) stop("Dataset still contains NA values!")

#Apply SMOTE
library(smotefamily)

X <- diabetic_data[, -which(names(diabetic_data) %in% c("readmitted"))]  # Features
Y <- as.factor(diabetic_data$readmitted)  # Target variable (ensure it's a factor)

# Over-sampling percentage
percOver <- max(0, (sum(Y == 0) - sum(Y == 1)) / sum(Y == 1) * 100)
table(Y)
smote_result <- SMOTE(X = X, target = Y, K = 5, dup_size = percOver / 100)

#Combine Balanced Data
balanced_data <- data.frame(smote_result$data)
colnames(balanced_data)[ncol(balanced_data)] <- "readmitted"
balanced_data$readmitted <- as.factor(balanced_data$readmitted)

# Check the results
print(table(balanced_data$readmitted))

#RF
library(randomForest)
library(caTools)

# Calculate class distribution
counts <- table(data_predict$readmitted)

# Plot the class distribution
barplot(counts, 
        main = "Class Distribution in Balanced Data",
        names.arg = c("No", "Yes"),
        xlab = "Classes",
        ylab = "Number of Cases",
        col = c("orange", "darkgreen"),
        legend = TRUE)

# Assuming data_predict is already balanced
# Calculate the total number of samples
num_samples <- nrow(data_predict)

# Randomly split the data into training (90%) and test (10%) sets
set.seed(123)  # Set seed for reproducibility
test_id <- sample(num_samples, size = 0.1 * num_samples)

# Split the data
data_train <- data_predict[-test_id, ]  # 90% Training data
data_test  <- data_predict[ test_id, ]  # 10% Testing data

rf <- randomForest(readmitted~., 
                   data = data_train,
                   na.action = na.omit,
                   importance = TRUE,
                   do.trace=10)
print(rf)
test_pred <- predict(rf, data_test)
rf_cfm <- caret::confusionMatrix(test_pred, data_test$readmitted)
rf_cfm
fourfoldplot(rf_cfm$table, color = c("#CC6666", "#99CC99"),
             conf.level = 0, margin = 1, main = "RF - Confusion Matrix (ACC: 90.86%)")
test_prob <- predict(rf, data_test, type = "prob")
par(pty="s")
roc_rf <- roc(data_test$readmitted ~ test_prob[,2], plot = TRUE, print.auc = TRUE, main="ROC curve of random forest")

#feature selection RFE
set.seed(123)
# use 10% of the training data for feature selection
num_train <- nrow(data_train)
sele_id <- sample(num_train, 0.1*num_train)
data_train.sele <- data_train[sele_id,]

control <- rfeControl(functions = rfFuncs, method = "cv", number=10,verbose=TRUE)
rf.select <- rfe(data_train.sele[complete.cases(data_train.sele), ][, -length(colnames(data_train.sele))], 
                 data_train.sele[complete.cases(data_train.sele), ][, length(colnames(data_train.sele))], 
                 sizes=seq(1,45,by=10), rfeControl=control, verbose=TRUE)
rf.select

plot(rf.select, type=c("g", "o"), cex=1,main="RFE")
RFE_feature <- predictors(rf.select)
RFE_feature

#Feature Selection Boruta
set.seed(123)
data_boruta <- Boruta(readmitted~., data=data_train.sele[complete.cases(data_train.sele), ], doTrace=0)
print(data_boruta)
par(mar=c(8,8,4,4),oma=c(0.1,0.1,0.1,0))
plot(data_boruta, xlab="", xaxt="n",cex.lab=3,cex.axis=2)
title("Boruta Feature Selection")
lz <- lapply(1:ncol(data_boruta$ImpHistory), function(i)
  data_boruta$ImpHistory[is.finite(data_boruta$ImpHistory[, i]), i])
names(lz) <- colnames(data_boruta$ImpHistory)
lb <- sort(sapply(lz, median))
axis(side=1, las=2, labels=names(lb), at=1:ncol(data_boruta$ImpHistory), cex.axis=0.5, font = 4)
boruta_feature <- getSelectedAttributes(data_boruta, withTentative = F)

#LR
# select important features and remove rows with missing value
data_train.prune <- data_train[,c(boruta_feature,"readmitted")]
data_test.prune <- data_test[,c(boruta_feature,"readmitted")]
logit <- glm(readmitted ~ ., data = data_train.prune, family = "binomial", 
             na.action = na.omit)
summary(logit)
selectedFeatures <- rownames(summary(logit)$coefficients)[summary(logit)$coefficients[,4] < 0.001 & summary(logit)$coefficients[,1] > 0]
selectedFeatures
selectedFeatures <- rownames(summary(logit)$coefficients)[summary(logit)$coefficients[,4] < 0.001 & summary(logit)$coefficients[,1] < 0]
selectedFeatures

# remove rows with missing value
data_train.prune <- data_train.prune[complete.cases(data_train.prune), ]
data_test.prune <- data_test.prune[complete.cases(data_test.prune), ]
logit.prob <- predict(logit, data_test.prune, type="response")
logit.pred <- factor(ifelse(logit.prob>0.5,"YES","NO"))
levels(logit.pred) <- levels(data_test.prune$readmitted)
logit.pred <- factor(logit.pred, levels = levels(data_test.prune$readmitted))

logit_test_cfm <- caret::confusionMatrix(logit.pred, data_test.prune$readmitted)
logit_test_cfm
lr_cfm <- caret::confusionMatrix(test_pred, data_test$readmitted)
lr_cfm
fourfoldplot(lr_cfm$table, color = c("#CC6666", "#99CC99"),
             conf.level = 0, margin = 1, main = "LR - Confusion Matrix (ACC: 90.72%)")
par(pty="s")
roc_lg <- roc(data_test.prune$readmitted ~ logit.prob, plot = TRUE, print.auc = TRUE, main="ROC curve of logistic regression")

#comparison
compare_model <- data.frame(Model=c("Random Forest", "Logistic Regression"),
                            Accuracy = c(0.9086, 0.9072),
                            #Kappa = c(0.861, 0.5889),
                            Sensitivity = c(1.00, 0.9982),
                            #Specificity = c(0.9044, 0.8298),
                            AUC = c(0.659, 0.674))
knitr::kable(compare_model)

par(pty="s")
plot(roc_rf, col = "green", main = "ROC for two models")
lines(roc_lg, col = "blue")
legend("bottomright", c("Random Forest (0.659)", "Logistic Regression (0.674)"), fill=c("green", "blue"),bty = "n")
