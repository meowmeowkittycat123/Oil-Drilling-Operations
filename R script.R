library(rpart)
library(rpart.plot)
library(data.table)
library(ggplot2)
install.packages("car")
library(car)
library(caTools)

setwd("/Users/jesslynlee/Documents/UNIVERSITY/Year 2/BC2406 Analytics 1/AY23 BC2406 CBA")
dt<-fread("marun_sample2.csv")
summary(dt)

#METHOD 1 
rows_with_missing_values <- rowSums(is.na(dt)) > 0
missingvalues<-which(rows_with_missing_values)
missingvalues

#METHOD 2 
# Calculate the median of FAN600, FAN300, MIN10GEL, and MUDLOSSU
median_FAN600 <- median(dt$FAN600, na.rm = TRUE)
median_FAN300 <- median(dt$FAN300, na.rm = TRUE)
median_MIN10GEL <- median(dt$MIN10GEL, na.rm = TRUE)
median_MUDLOSSU <- median(dt$MUDLOSSU, na.rm = TRUE)

# Replace missing values with the calculated medians
dt$FAN600[is.na(dt$FAN600)] <- median_FAN600
dt$FAN300[is.na(dt$FAN300)] <- median_FAN300
dt$MIN10GEL[is.na(dt$MIN10GEL)] <- median_MIN10GEL
dt$MUDLOSSU[is.na(dt$MUDLOSSU)] <- median_MUDLOSSU

median_FAN600
median_FAN300
median_MIN10GEL
median_MUDLOSSU


# Verify that missing values have been imputed

summary(dt)
nrow(dt)
cor(dt)

# Count the number of data points for each formation type
formation_counts <- table(dt$Formation)
print(formation_counts)

#Bar Chart
ggplot(dt, aes(x = Formation, y = MUDLOSSU, fill= Formation)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2)

#Scatter Plot
ggplot(dt, aes(x = Formation, y = MUDLOSSU, colour=Formation)) +
  geom_point()

# 70% trainset 30% testset
set.seed(2)
train<-sample.split(dt$MUDLOSSU,SplitRatio=0.7)
trainset<-subset(dt, train==TRUE)
testset<-subset(dt, train==FALSE)
dim(trainset) #Means that there are 1885 data points and 20 columns

#Linear Regression
# Step 1: Start with the full model
m0 <- lm(MUDLOSSU ~ ., data = dt)
vif(m0)

# Calculate Cook's distance of initial model
cooksd <- cooks.distance(m0)
plot(cooksd, pch = 19, cex = 1, main = "Cook's Distance Plot")
abline(h = 4 / length(cooksd), col = "red")

# Step 2: Backward elimination
set.seed(2)
dt$Formation <- factor(dt$Formation)
dt$Formation
m <- lm(MUDLOSSU ~ ., data = trainset)  # . means all other variables
m0 <- step(m)  # Akaike Information Criterion

# Step 3: Linear regression with selected variables after elimination
subset_data <- trainset[, c("MUDLOSSU", "Formation", "Pore pressure","Fracture pressure", "Mud pressure (psi)", "Hole size (in)", "METERAGE", "WOB", "Pump flow rate", "MFVIS", "RETSOLID", "FAN300","FAN600")]

# Fit the linear regression model with the selected predictors
m1 <- lm(MUDLOSSU ~ ., data = subset_data)
summary(m1)
vif(m1)


#Removing the variables based on high vif  - Removed FAN600 and Fracture Pressure and Mud pressure (psi)
subset_data2 <- trainset[, c("MUDLOSSU", "Easting", "Formation", "Pore pressure",
                            "Hole size (in)", "METERAGE", "WOB", "Pump flow rate", 
                            "MFVIS", "RETSOLID", "FAN300")]
m2 <- lm(MUDLOSSU ~ ., data = subset_data2)
vif(m2)
summary(m2)

#Removing more variables based on significance - Removed Easting, Pore Pressure and WOB
subset_data3 <- trainset[, c("MUDLOSSU", "Formation",
                             "Hole size (in)", "METERAGE", "Pump flow rate", 
                             "MFVIS", "RETSOLID", "FAN300")]
m3 <- lm(MUDLOSSU ~ ., data = subset_data3)
vif(m3)
summary(m3)

# Calculate Cook's distance of Optimal model 
cooksd <- cooks.distance(m3)

# Create a Cook's distance plot
plot(cooksd, pch = 19, cex = 1, main = "Cook's Distance Plot")
abline(h = 4 / length(cooksd), col = "red")

# RMSE for the training set
RMSE.m3.train <- sqrt(mean(residuals(m3)^2))  # RMSE on trainset based on m3 model.
summary(abs(residuals(m3)))  # Check Min Abs Error and Max Abs Error.

# RMSE for the test set
predict.m3.test <- predict(m3, newdata = testset)
testset.error <- testset$MUDLOSSU - predict.m3.test

RMSE.m3.test <- sqrt(mean(testset.error^2))
summary(abs(testset.error))

RMSE.m3.train 
RMSE.m3.test 


#CART MODEL 
cart<-rpart(MUDLOSSU~. ,method = 'anova', trainset, control = rpart.control(minsplit = 2,cp=0))
dt <- data.table(cart$cptable)
dt[, index := 1:nrow(dt)]
min_cp_index <- min(dt[(xerror+xstd == min(xerror+xstd)), index])
errorcap<- dt[min_cp_index, xerror+xstd]
optimal_cp_index<- min(dt[(xerror< errorcap), index])
cp.opt=sqrt(dt[index== optimal_cp_index, CP]*dt[index==optimal_cp_index-1,CP])
m.opt<-prune(cart, cp=cp.opt)
printcp(m.opt)
plot_obj<-rpart.plot(m.opt, cex=0.9)

# RMSE for the training set
trainset.predict <- predict(m.opt, newdata = trainset, type = "vector")
rmse_train <- sqrt(mean((trainset$MUDLOSSU - trainset.predict)^2))

# RMSE for the test set
testset.predict <- predict(m.opt, newdata = testset, type = "vector")
rmse_test <- sqrt(mean((testset$MUDLOSSU - testset.predict)^2))

rmse_train
rmse_test   

#Variable Importance
m.opt$variable.importance
summary(m.opt)








