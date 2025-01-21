
######cohort1######
#DEgene预测特征筛选：DEgene：|log2FC|>1.5，p<0.05
####0.数据特征描述####
library(readxl)
library(tableone)
library(summarytools)

V1_clin <- read_excel("seq/V1-V2-delta-clin - 副本.xlsx")
str(V1_clin)
V1_clin$Remission <- ifelse(V1_clin$Remission == "NR", 0, ifelse(V1_clin$Remission == "R", 1, V1_clin$Remission))
V1_clin$Remission <- as.factor(V1_clin$Remission)
V1_clin$Gender <- as.factor(V1_clin$Gender)

library(dplyr)
V1_clin_filled <- V1_clin %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
data <- V1_clin_filled
sum(is.na(data))

# 存储不符合正态分布的列名
biomarkers <- c()
# 进行正态性检验，排除Remission和Gender列
for (column in names(data)) {
  if (!(column %in% c("Remission", "Gender"))) {
    # 进行Shapiro-Wilk正态性检验
    shapiro_test <- shapiro.test(data[[column]])
    # 检查p值，小于0.05则认为不符合正态分布
    if (shapiro_test$p.value < 0.05) {
      biomarkers <- c(biomarkers, column)
    }
  }
}


library(tableone)
library(survival)

catVars <- c("Gender")
myVars <- setdiff(names(data), catVars)

tab1 <- CreateTableOne(vars = myVars, data = data, factorVars = catVars)
print(tab1, nonnormal = biomarkers,showAllLevels = TRUE)
tab1Mat <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal = biomarkers)

tab2 <- CreateTableOne(vars = myVars, strata = "Remission" , data = data, factorVars = catVars)
print(tab2, nonnormal = biomarkers,showAllLevels = TRUE)
tab2Mat <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal = biomarkers)

tabMat <- cbind(tab1Mat, tab2Mat)
library(openxlsx)
write.xlsx(tabMat, file = "output.xlsx", rowNames = TRUE)

####1.lasso特征筛选####
library(readxl)
library(glmnet)
library(foreign)
library(survival)

seq_V1_V2_TPM <- t(read_excel("seq/seq-V1-V2-TPM.xlsx"))
colnames(seq_V1_V2_TPM) <- unlist(seq_V1_V2_TPM[1, ])
data <- seq_V1_V2_TPM[-1, ]

numeric_matrix_data <- matrix(as.numeric(data), nrow = nrow(data))
colnames(numeric_matrix_data) <- colnames(data)
rownames(numeric_matrix_data) <- rownames(data)
colnames(numeric_matrix_data) <- gsub("-", ".", colnames(numeric_matrix_data))

#numeric_matrix_data <- as.data.frame(numeric_matrix_data[-9,])


#列均值插补
library(mice)
library(Hmisc)
for (i in 1:ncol(numeric_matrix_data)) {
  numeric_matrix_data[, i] <- impute(numeric_matrix_data[, i], median)}
sum(is.na(numeric_matrix_data))


#列中位数插补
library(dplyr)
numeric_matrix_data_filled <- numeric_matrix_data %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
numeric_matrix_data <- numeric_matrix_data_filled


x <-numeric_matrix_data[,2:89]
X = scale(x,center = T,scale = T)
#X变量应该用均值零和单位方差进行标准化，因为变量的尺度差异往往会使惩罚分配不均。
y <- numeric_matrix_data[,1]

cvfit <- glmnet(x=X, y=y, family="binomial",nlambda = 100,alpha = 1,nfolds = 8)
print(cvfit)
plot(cvfit,xvar="lambda", label=F,lwd = 1.75)

cv_fit <- cv.glmnet(x=X, y=y, family="binomial",nlambda = 100,alpha = 1,nfolds =8)#默认为10折交叉验证
plot(cv_fit)
cv_fit$lambda.min
cv_fit$lambda.1se

model_lasso_min <- glmnet(x=X, y=y, family="binomial",alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=X, y=y, family="binomial",alpha = 1, lambda=cv_fit$lambda.1se)

choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
model_lasso_min$beta
model_lasso_1se$beta
length(choose_gene_min)
length(choose_gene_1se)
choose_gene_min
choose_gene_1se

#将choose_gene_min、choose_gene_1se合并为一个2列的数据框输出
choose_gene_min_df <- data.frame(choose_gene_min)
choose_gene_1se_df <- data.frame(choose_gene_1se)
#将choose_gene_min、choose_gene_1se合并输出
choose_gene_df <- cbind(choose_gene_min_df, choose_gene_1se_df)
write.csv(choose_gene_df, "LSAAO-choose_gene_df.csv")


lasso.prob <- predict(cv_fit,  #predict是预测，是一个泛型函数
                      newx=X ,  #这个是预测对象，cv_fit是预测数据
                      s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
#将每个样本的生死和预测结果放在一起
re=cbind(y ,lasso.prob) 
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)  #作为因子来变成分类变量
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2


lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se

# Fit models with lambda.min and lambda.1se
model_min <- glmnet(x = X, y = y, family = "binomial", alpha = 1, lambda = lambda_min)
model_1se <- glmnet(x = X, y = y, family = "binomial", alpha = 1, lambda = lambda_1se)
lasso_prob <- predict(cv_fit, newx = X, s = c(lambda_min, lambda_1se))
re <- data.frame(event = as.factor(y), prob_min = lasso_prob[, 1], prob_1se = lasso_prob[, 2])
# Calculate ROC curves and AUC
pred_min <- prediction(re$prob_min, re$event)
perf_min <- performance(pred_min, "tpr", "fpr")
auc_min <- performance(pred_min, "auc")@y.values[[1]]

pred_1se <- prediction(re$prob_1se, re$event)
perf_1se <- performance(pred_1se, "tpr", "fpr")
auc_1se <- performance(pred_1se, "auc")@y.values[[1]]

# Plot ROC curves using ggplot2
ggplot() +
  geom_line(data = data.frame(tpr = perf_min@y.values[[1]], fpr = perf_min@x.values[[1]]),
            aes(x = fpr, y = tpr), color = "#0072B5FF", size = 1.5) +
  geom_line(data = data.frame(tpr = perf_1se@y.values[[1]], fpr = perf_1se@x.values[[1]]),
            aes(x = fpr, y = tpr), color = "#BC3C29FF", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.6, y = 0.3, label = paste0("AUC.min = ", round(auc_min, 3)), color = "#0072B5FF", size = 4.5, hjust = 0) +
  annotate("text", x = 0.6, y = 0.2, label = paste0("AUC.1se = ", round(auc_1se, 3)), color = "#BC3C29FF", size = 4.5, hjust = 0) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "ROC Curves for Lasso Models") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)  # 居中标题，设置字体大小为20
  )


####2.RF特征筛选######
library(randomForest)
library(varSelRF)
library(pROC)
library(caret)
library(readxl)

raw_data <- t(read_excel("seq/seq-V1-V2-TPM.xlsx"))
colnames(raw_data) <- unlist(raw_data[1, ])
raw_data <- raw_data[-1, ]
df <- as.data.frame(raw_data)
df[] <- lapply(df, as.numeric)

#列均值插补
library(mice)
library(Hmisc)
for (i in 1:ncol(df)) {
  df[, i] <- impute(df[, i], median)}
sum(is.na(df))
df[] <- lapply(df, as.numeric)

data <- df
colnames(data) <- gsub("-", ".", colnames(data))
data$Remission <- factor(data$Remission)

sum(is.na(data))
str(data)

set.seed(1)
train <- sample(nrow(data),0.7*nrow(data))
train <- createDataPartition(data$Remission,p=0.7,list = F)
test <- data[-train,]
train <- data[train,]

set.seed(2)
rf.train <- randomForest(Remission ~ ., data=train, importance=TRUE,mtry=17, proximity=TRUE) 
#数据集特征少，尝试较小的mtry值，如sqrt(p)或者log2(p),p=特征数量,以确保每棵树的分裂是基于数据集中的少量特征
#数据集特征多，尝试较大的mtry值，如p/3或者p/2，这样可以增加每棵树的随机性，防止过拟合
#如果不确定应该选择多大的mtry值，可以通过交叉验证尝试不同的mtry值，然后选择在交叉验证中表现最好的那个
#importance参数是看样本重要性的，proximity是样本之间距离

print(rf.train)
plot(rf.train,main = "RF origin")


# 从随机森林模型中提取变量重要性
importance_values <- rf.train$importance
variable_names <- rownames(importance_values)
importance_df <- data.frame(importance_values)


library(openxlsx)
excel_file <- "variable_importance_v1+2.xlsx"
write.xlsx(importance_df, file = excel_file, rowNames = TRUE)



importance_scores <- importance(rf.train)
feature_importance_gini <- data.frame(Gene = rownames(importance_scores), 
                                      Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features_gini <- feature_importance_gini[order(-feature_importance_gini$Importance), ]

feature_importance_accuracy <- data.frame(Gene = rownames(importance_scores), 
                                          Importance = importance_scores[, "MeanDecreaseAccuracy"])
ordered_features_accuracy <- feature_importance_accuracy[order(-feature_importance_accuracy$Importance), ]

importance_gini <- importance(rf.train, type = 2)
importance_accuracy <- importance(rf.train, type = 1)

library(ggplot2)
library(ggthemes)


top_genes <- head(ordered_features_accuracy, 15)
color_range <- c("#ADD8E6", "#00008B")

ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) + # 使用渐变填充
  theme_few() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    plot.caption = element_text(size = 11), 
    axis.text = element_text(size = 11), 
    axis.title = element_text(size = 14)
  ) +
  labs(title = "Top 15 important miRNA in RF Prediction", x = "miRNA", y = "MeanDecreaseAccuracy") +
  coord_flip() # 横向柱状图更易读

top_genes <- head(ordered_features_gini, 15)
color_range <- c("#D8BFD8", "#8B008B")

ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) + # 使用渐变填充
  theme_few() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    plot.caption = element_text(size = 11), 
    axis.text = element_text(size = 11), 
    axis.title = element_text(size = 14)
  ) +
  labs(title = "Top 15 important miRNA in RF Prediction", x = "miRNA", y = "MeanDecreaseGini") +
  coord_flip() # 横向柱状图更易读


# 对训练数据进行预测
predictions <- predict(rf.train, train, type = "response")
prediction_data <- data.frame(Actual = train$Remission, Predicted = predictions)# 创建一个数据框，包含实际值和预测值
correct_predictions <- prediction_data$Actual == prediction_data$Predicted# 计算预测准确率
accuracy <- sum(correct_predictions) / nrow(prediction_data)
hist(accuracy, main = "Prediction Accuracy Histogram", xlab = "Accuracy")# 绘制预测准确率的直方图

set.seed(4)
pred_class<-predict(rf.train,newdata=train,type = "response")  #输出类别
pred_caret <- caret::confusionMatrix(as.factor(pred_class),as.factor(train$Remission))
#一次性输出AUC、混淆矩阵、敏感度特异度PPV NPV等模型评价指标
pred_caret

# 对测试数据进行预测
set.seed(5)
pred_class<-predict(rf.train,newdata=test,type = "response")  #输出类别
pred_caret <- caret::confusionMatrix(as.factor(pred_class),as.factor(test$Remission))
#一次性输出AUC、混淆矩阵、敏感度特异度PPV NPV等模型评价指标
pred_caret

rf.train <- randomForest(Remission ~ ., data=train, importance=TRUE,mtry=17, proximity=TRUE) 
# 对所有数据进行预测
set.seed(6)
pred_class<-predict(rf.train,newdata=data,type = "response")  #输出类别
pred_caret <- caret::confusionMatrix(as.factor(pred_class),as.factor(data$Remission))
#一次性输出AUC、混淆矩阵、敏感度特异度PPV NPV等模型评价指标
pred_caret

pred_prob<-predict(rf.train,newdata=data,type="prob")  #输出概率
library(pROC)
roc_curve <- roc(data$Remission, pred_prob[,2])
roc_data <- data.frame(
  fpr = roc_curve$specificities,
  tpr = roc_curve$sensitivities
)

auc_value <- auc(roc_curve)

ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "#0072B5FF", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.6, y = 0.3, label = paste0("AUC = ", round(auc_value, 3)), color = "#0072B5FF", size = 4.5, hjust = 0) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curve for Random Forest Model") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  )

####3.smv-REF特征选择####
#方法三：SVM-RFE
# 标准化或归一化数据
scaled_data <- scale(data[, -which(names(data) == "Remission")])
scaled_data <- cbind(scaled_data, data["Remission"])

# 定义 SVM-RFE 函数
svmFuncs <- caretFuncs
svm_rfe <- rfeControl(functions = svmFuncs, method = "cv", number = 10)
# 使用 SVM-RFE 进行特征选择
feature_select <- rfe(x = scaled_data[, -which(names(scaled_data) == "Remission")],
                      y = scaled_data$Remission,
                      sizes = c(1:ncol(scaled_data) - 1),
                      rfeControl = svm_rfe)
print(feature_select)
summary(feature_select)

library(openxlsx)
write.xlsx(results, file = "SMV-REF-results.xlsx")


results <- data.frame(Features = rownames(varImp(feature_select)),
                      Importance = varImp(feature_select)$Overall)
results <- results[order(results$Importance, decreasing = TRUE), ][1:9, ]
ggplot(results, aes(x = reorder(Features, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#ADD8E6", width = 0.85) +
  labs(title = "Top 9 Feature Importance using SVM-RFE",
       y = "Importance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,face = "bold",size = 9),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.5))


# 绘制变量个数与误差之间的关系图
plot(feature_select, type=c("g","o"))

results <- feature_select$results
# 绘制变量个数与误差之间的关系图
ggplot(results, aes(x = Variables, y = Accuracy)) +
  geom_line(color = "#BC3C29FF",linewidth = 1) +
  geom_point(color = "#BC3C29FF",size=1.5) +
  labs(title = "Number of Features vs Accuracy",
       x = "Number of Features",
       y = "RMSE (Root Mean Squared Error)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)
  )


# 获取特征选择后的数据集
selected_indices <- feature_select$optVariables
selected_data <- scaled_data[, c(selected_indices)]
remission_column <- data[["Remission"]]
merged_data <- cbind(selected_data, Remission = remission_column)
selected_data <- merged_data

set.seed(12)
train_index <- createDataPartition(selected_data$Remission, p = 0.7, list = FALSE)
train_data <- selected_data[train_index, ]
test_data <- selected_data[-train_index, ]

svm_model <- svm(Remission ~ ., data = train_data, kernel = "linear", probability = TRUE)
#kernel参数用于指定核函数的类型。核函数决定了数据在高维空间中的映射方式，进而影响了模型在原始特征空间中的表现。常见的核函数包括线性核（linear）、多项式核（polynomial）、高斯径向基函数核（radial basis function, RBF）等
summary(svm_model)
predicted_labels <- predict(svm_model, newdata = test_data[, -which(names(test_data) == "Remission")])

conf_matrix <- table(test_data$Remission, predicted_labels)
conf_matrix
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
accuracy

# 计算精确率
precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
print(paste("Precision:", precision))
# 计算召回率
recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
print(paste("Recall:", recall))
# 计算 F1 分数
f1_score <- 2 * precision * recall / (precision + recall)
print(paste("F1 Score:", f1_score))


library(ROCR)
# 获取SVM模型的预测概率
svm_prob <- attr(predict(svm_model, newdata = test_data[, -which(names(test_data) == "Remission")], probability = TRUE), "probabilities")
svm_pred <- prediction(svm_prob[,2], test_data$Remission)
svm_perf <- performance(svm_pred, "tpr", "fpr")

# 提取 ROC 曲线数据
roc_curve <- data.frame(fpr = unlist(svm_perf@x.values[[1]]), tpr = unlist(svm_perf@y.values[[1]]))
svm_auc <- performance(svm_pred, "auc")@y.values[[1]]
# 绘制 ROC 曲线图
ggplot(roc_curve, aes(x = fpr, y = tpr)) +
  geom_line(color = "#0072B5FF", size = 1.5) +  # 第一条 ROC 曲线，设置颜色和线条粗细
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +  # 对角线，表示随机分类器
  annotate("text", x = 0.75, y = 0.3, label = paste0("AUC = ", round(svm_auc, 3)), color = "#0072B5FF", size = 4.5, hjust = 0) +  # 注释 AUC 值，设置文本颜色、大小和位置
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curve for SVM Model") +  # 设置坐标轴标签和标题
  theme_minimal() +  # 使用最小化的主题风格
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)  # 居中标题，设置字体大小
  )


######cohort2临床模型构建######
library(readxl)
V1_clin <- read_excel("PCR/V1-V2-Delta - 副本.xlsx")
str(V1_clin)
V1_clin$Remission <- ifelse(V1_clin$Remission == "NR", 0, ifelse(V1_clin$Remission == "R", 1, V1_clin$Remission))
V1_clin$Remission <- as.factor(V1_clin$Remission)
V1_clin$Gender <- as.factor(V1_clin$Gender)

library(dplyr)
V1_clin_filled <- V1_clin %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
data <- V1_clin_filled[, -1]
sum(is.na(data))

# 存储不符合正态分布的列名
biomarkers <- c()
# 进行正态性检验，排除Remission和Gender列
for (column in names(data)) {
  if (!(column %in% c("Remission", "Gender"))) {
    # 进行Shapiro-Wilk正态性检验
    shapiro_test <- shapiro.test(data[[column]])
    # 检查p值，小于0.05则认为不符合正态分布
    if (shapiro_test$p.value < 0.05) {
      biomarkers <- c(biomarkers, column)
    }
  }
}


library(tableone)
library(survival)

catVars <- c("Gender")
myVars <- setdiff(names(V1_clin), catVars)

tab1 <- CreateTableOne(vars = myVars, data = data, factorVars = catVars)
print(tab1, showAllLevels = TRUE)
print(tab1, nonnormal = biomarkers)
tab1Mat <- print(tab1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal = biomarkers)

tab2 <- CreateTableOne(vars = myVars, strata = "Remission" , data = data, factorVars = catVars)
print(tab2, nonnormal = biomarkers)
tab2Mat <- print(tab2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal = biomarkers)

tabMat <- cbind(tab1Mat, tab2Mat)
library(openxlsx)
write.xlsx(tabMat, file = "output.xlsx", rowNames = TRUE)



#####1.临床模型构建#####
library(readxl)
V1_clin <- read_excel("PCR/V1-V2-Delta.xlsx")
#将remission NR改为0 R改为1
V1_clin$Remission <- ifelse(V1_clin$Remission == "NR", 0, ifelse(V1_clin$Remission == "R", 1, V1_clin$Remission))
#转化为因子变量
V1_clin$Remission <- as.factor(V1_clin$Remission)
str(V1_clin$Remission)
V1_clin$Gender <- as.factor(V1_clin$Gender)
V1_clin$Group <- as.factor(V1_clin$Group)
V1_clin$`family history of T2DM` <- as.factor(V1_clin$`family history of T2DM`)
colnames(V1_clin) <- gsub("/", ".", colnames(V1_clin))
names(V1_clin)

#使用mutate_if函数对所有数值变量进行中值插补
library(dplyr)
V1_clin_filled <- V1_clin %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
#将V1_clin_filled中的姓名列列名改为ID
colnames(V1_clin_filled)[1] <- "ID"
data <- V1_clin_filled
sum(is.na(data))


#t检验*：V2-glu0 ，V2-AUCglu ，V2-glu120 ，V1-Weight ，V1-UA 
#SLR*：V2-glu0，V2-AUCglu，V2-ISI，delta-ISI

formular1 <- Remission ~ `V1-UA` +`V1-Weight` +`delta-ISI`+`V2-ISI`+`V2-glu0`+`V2-AUCglu`
fit1 <- glm(formular1, 
           data = data, 
           family = binomial(link = "logit"))
summary(fit1)

library(car)
vif_values <- vif(fit1)
print(vif_values)

library(MASS)
model.both <- stepAIC(fit1,direction = "both")
summary(model.both)$coefficients

#删除与最关键变量V2-ISI共线的delta-ISI后，再次拟合模型
formular2 <- Remission ~ `V1-UA` +`V1-Weight` +`V2-ISI`
formular2 <- Remission ~ `V1-UA` +`V1-Weight` +`V2-ISI`+`V1-BMI`+`Age`+`Gender`+`V1-HBA1C`
fit2 <- glm(formular2, 
           data = data, 
           family = binomial(link = "logit"))
summary(fit2)

#是否加入V2-FBG相关变量对于模型拟合无影响，遂无需考虑引入冗余变量
formular3 <- Remission ~ `V1-UA` +`V1-Weight` +`V2-ISI`+`V2-glu0`+`V2-AUCglu`
fit3 <- glm(formular3, 
            data = data, 
            family = binomial(link = "logit"))
summary(fit3)

formular4 <- Remission ~ `V1-UA` +`V1-Weight` +`V2-ISI`+`V2-glu0`
fit4 <- glm(formular4, 
            data = data, 
            family = binomial(link = "logit"))
summary(fit4)

#最终临床特征模型为fit2：AUC0.7496
library(pROC)
predicted_probabilities <- predict(fit2, type = "response")
roc_data <- roc(data$Remission, predicted_probabilities)
auc(roc_data)


#####2.miR模型构建#####
library(readxl)
PCR_result_5 <- read_excel("D:/F盘迁移文件/研究生/7.课题/microRNA/文章/PCR-result-5.xlsx")

data_pcr <- PCR_result_5[,c(1:2,8:ncol(PCR_result_5))]
data_pcr$Remission <- factor(data_pcr$Remission)
str(data_pcr)

##对data_pcr数据标准化，除外Remission列
features <- data_pcr[, !(names(data_pcr) %in% c("Remission","ID"))]
normalized_features <- scale(features)
data_pcr[, !(names(data_pcr) %in% c("Remission","ID"))] <- normalized_features

#将data_pcr中所有列名中包含“-”替换为“.”
colnames(data_pcr) <- gsub("-", ".", colnames(data_pcr))


#对data_pcr中所有的变量做单因素logistic回归,y=Remission
independent_vars <- setdiff(names(data_pcr), c("Remission","ID"))
results <- data.frame()
for (var in independent_vars) {
  # 构建模型公式，使用逻辑回归 glm 函数
  formula <- as.formula(paste("Remission ~", var))
  model <- glm(formula, data = data_pcr, family = binomial)
  
  # 提取回归结果的汇总信息
  summary_info <- summary(model)$coefficients
  # 提取估计值和p值
  estimate <- summary_info[2, 1]
  p_value <- summary_info[2, 4]
  
  # 将结果添加到结果数据框中
  results <- rbind(results, data.frame(
    Variable = var,
    Estimate = estimate,
    P_value = p_value
  ))
}
print(results)

formula <- as.formula(paste("Remission ~", paste(independent_vars, collapse = " + ")))
#回归法寻找最佳模型
m1 <- glm(formula,family = binomial, data = data_pcr)
summary(m1)

library(MASS)
model.both <- stepAIC(m1,direction = "both")
summary(model.both)
model_miR <- model.both
#AUC=0.702

m2 <- glm(Remission ~  V1.miR.502.3p, data = data_pcr, family = binomial())
summary(m2)
m2$coefficients
model_miR <- m2
#m2为仅包含p显著的3个miRNA的模型，AUC=0.6709
#若仅包含502-3p，则AUC=0.6207

predicted_probabilities <- predict(model_miR, type = "response")
library(pROC)
roc_data <- roc(data_pcr$Remission, predicted_probabilities)
auc(roc_data)

#####3.多个模型比较#####
#97人的临床数据
library(readxl)
V1_clin <- read_excel("PCR/V1-V2-Delta.xlsx")
#将remission NR改为0 R改为1
V1_clin$Remission <- ifelse(V1_clin$Remission == "NR", 0, ifelse(V1_clin$Remission == "R", 1, V1_clin$Remission))
#转化为因子变量
V1_clin$Remission <- as.factor(V1_clin$Remission)
str(V1_clin$Remission)
V1_clin$Gender <- as.factor(V1_clin$Gender)
V1_clin$Group <- as.factor(V1_clin$Group)
V1_clin$`family history of T2DM` <- as.factor(V1_clin$`family history of T2DM`)
colnames(V1_clin) <- gsub("/", ".", colnames(V1_clin))
names(V1_clin)

#使用mutate_if函数对所有数值变量进行中值插补
library(dplyr)
V1_clin_filled <- V1_clin %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .))
#将V1_clin_filled中的姓名列列名改为ID
colnames(V1_clin_filled)[1] <- "ID"
data <- V1_clin_filled
sum(is.na(data))

#查看data的列名
names(data)
#将V1-matsyda index和V2-matsyda index都除以18再重新赋值
data$`V1-matsyda index` <- data$`V1-matsyda index`/18
data$`V2-matsyda index` <- data$`V2-matsyda index`/18

#将 "V1-ISI"重新赋值为"V1-matsyda index"*"V1-AUCins.glu"
data$`V1-ISI` <- data$`V1-matsyda index`*data$`V1-AUCins.glu`
#将 "V2-ISI"重新赋值为"V2-matsyda index"*"V2-AUCins.glu"
data$`V2-ISI` <- data$`V2-matsyda index`*data$`V2-AUCins.glu`


PCR_result_5 <- read_excel("D:/F盘迁移文件/研究生/7.课题/microRNA/文章/PCR-result-5.xlsx")
data_pcr <- PCR_result_5[,c(1:2,8:ncol(PCR_result_5))]
data_pcr$Remission <- factor(data_pcr$Remission)
str(data_pcr)

##对data_pcr数据标准化，除外Remission列
features <- data_pcr[, !(names(data_pcr) %in% c("Remission","ID"))]
normalized_features <- scale(features)
data_pcr[, !(names(data_pcr) %in% c("Remission","ID"))] <- normalized_features

#将data_pcr中所有列名中包含“-”替换为“.”
colnames(data_pcr) <- gsub("-", ".", colnames(data_pcr))

#将V1_clin_filled与data_pcr根据ID列合并，且只保留一列Remission
merged_data <- merge(V1_clin_filled, data_pcr, by = "ID")
merged_data <- merged_data[, -which(names(merged_data) == "Remission.y")]
colnames(merged_data)[which(names(merged_data) == "Remission.x")] <- "Remission"

#将merged_data数据划分
library(caret)
set.seed(101)
train_index <- createDataPartition(merged_data$Remission, p = 0.7, list = FALSE)
train_merged_data <- merged_data[train_index, ]
test_merged_data <- merged_data[-train_index, ]

colnames(train_merged_data) <- gsub("-", ".", colnames(train_merged_data))
colnames(test_merged_data) <- gsub("-", ".", colnames(test_merged_data))



model_clinic <- glm(Remission ~ `V1.UA` +`V1.Weight` +`V2.ISI`, 
                    data = train_merged_data, family = binomial())
summary(model_clinic)

model_miR <- glm(Remission ~ V1.miR.502.3p ,data = train_merged_data, family = binomial())
summary(model_miR)


formular2 <- Remission ~ `V1.UA` +`V1.Weight` +`V2.ISI`+`V1.BMI`+`Age`+`Gender`+`V1.HBA1C`+ V1.miR.502.3p
model <- glm(formular2,
             data =train_merged_data, family = binomial())
summary(model)

model <- glm(Remission ~`V1.UA`+`V1.Weight` +`V2.ISI`+ V1.miR.502.3p  ,
             data =train_merged_data, family = binomial())
summary(model)


odds_ratios <- exp(coef(model_clinic))
odds_ratios
# 计算置信区间
confint_model <- confint(model_clinic)
or_confint <- exp(confint_model)
or_confint


#似然比检验
anova(model_clinic, model, test = "Chisq")
anova(model_miR, model, test = "Chisq")

library(report)
#生成模型的详细报告
report(model_clinic)
report(model)

library(easystats)
model_dashboard(model)


library(performance)
model_performance(model)
compare_performance(model_clinic, model_miR, model, rank = TRUE, verbose = FALSE)

model_clinical_factor <- model_clinic
model_miR <- model_miR
model_combined <- model

plot(compare_performance(model_clinical_factor, model_miR, model_combined, rank = TRUE, verbose = FALSE))


#####3.1ROC曲线#####
library(pROC)
library(ggplot2)
# 计算ROC曲线
roc_clinic <- roc(response = train_merged_data$Remission, predictor = predict(model_clinic), smooth = TRUE)
roc_miR <- roc(response = train_merged_data$Remission, predictor = predict(model_miR), smooth = TRUE)
roc_model <- roc(response = train_merged_data$Remission, predictor = predict(model), smooth = TRUE)

# 提取数据
roc_clinic_df <- data.frame(fpr = roc_clinic$specificities, tpr = roc_clinic$sensitivities, model = "Clinic")
roc_miR_df <- data.frame(fpr = roc_miR$specificities, tpr = roc_miR$sensitivities, model = "miR")
roc_model_df <- data.frame(fpr = roc_model$specificities, tpr = roc_model$sensitivities, model = "Model")

roc_df <- rbind(roc_clinic_df, roc_miR_df, roc_model_df)

# 计算AUC
auc_clinic <- auc(roc_clinic)
auc_miR <- auc(roc_miR)
auc_model <- auc(roc_model)

auc_train_clinic <- auc_clinic
auc_train_miR <- auc_miR
auc_train_model <- auc_model

# 绘图
ggplot(roc_df, aes(x = 1 - fpr, y = tpr, color = model)) +
  geom_line(size = 1.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.43, y = 0.35, label = paste0("AUC - Clinical Factor Model = ", round(auc_clinic, 3)), color = "#377DB8", size = 4.4, hjust = 0) +
  annotate("text", x = 0.43, y = 0.25, label = paste0("AUC - miR Model = ", round(auc_miR, 3)), color = "#4DAE49", size = 4.4, hjust = 0) +
  annotate("text", x = 0.43, y = 0.15, label = paste0("AUC - Combined Model = ", round(auc_model, 3)), color =  "#E31A1C", size = 4.4, hjust = 0) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curves Comparison") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 13),          
        axis.title.y = element_text(size = 13),          
        legend.title = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1)) +
  scale_color_manual(values = c("Clinic" = "#377DB8", "miR" = "#4DAE49", "Model" =  "#E31A1C"))+
  guides(color = "none") 

#####3.2校准曲线绘制####
library(rms)
dd=datadist(train_merged_data)
options(datadist="dd")

library(riskRegression)
par(mar = c(5, 4, 4, 2) + 0.1)
xb <- Score(list("Clinical Factor Model"=model_clinic,
                 "miR Model"=model_miR,
                 "Combined Model"=model),
            formula = Remission~1,
            null.model = F,
            conf.int = T,
            plots = c("calibration","ROC"),
            metrics = c("auc","brier"),
            B=1000,M=10,
            data=train_merged_data)

plotCalibration(xb,cex = 0.8)

#loess平滑曲线
library(ggplot2)
library(rms)
library(runway)

y <- train_merged_data[, "Remission"]

pr_clinic<- predict(model_clinic,type = c("response"))#得出预测概率
pr_miR<- predict(model_miR,type = c("response"))#得出预测概率
pr <- predict(model,type = c("response"))#得出预测概率

model_name <- rep(c("model_clinical"), length.out = 69)
df_c <- data.frame(
  outcomes = y,
  predictions = pr_clinic,
  model_name = model_name
)

model_name <- rep(c("model_miR"), length.out = 69)
df_m <- data.frame(
  outcomes = y,
  predictions = pr_miR,
  model_name = model_name
)

model_name <- rep(c("model"), length.out = 69)
df <- data.frame(
  outcomes = y,
  predictions = pr,
  model_name = model_name
)

multi_model_dataset <- rbind(df_c,df_m,df)


threshperf_plot_multi(multi_model_dataset,
                      outcome = 'outcomes',
                      prediction = 'predictions',
                      model = 'model_name',
                      positive = 1)

cal_plot_multi(multi_model_dataset,
               outcome = 'outcomes',
               prediction = 'predictions',
               model = 'model_name',
               n_bins = 0,
               show_loess = TRUE,
               positive = 1)



#####3.3.DCA曲线、CIC曲线绘制#######
library(rmda)

formula_clinic <- as.formula(Remission ~ `V1.UA` +`V1.Weight` +`V2.ISI`)
formula_miR <- as.formula(Remission ~ V1.miR.502.3p)
formula <- as.formula(Remission ~`V1.UA`+`V1.Weight` +`V2.ISI`+ V1.miR.502.3p)

#将merged_data$Remission.x列转化为数值
train_merged_data$Remission <- as.character(train_merged_data$Remission)
str(train_merged_data$Remission)
train_merged_data$Remission <- as.numeric(train_merged_data$Remission)


model_clinic_DCA <- decision_curve(formula = formula_clinic,data = train_merged_data,
                                   family = binomial(link = 'logit'),
                                   confidence.intervals = 0.95,study.design = 'cohort',)

model_miR_DCA <- decision_curve(formula = formula_miR,data = train_merged_data,
                                family = binomial(link = 'logit'),
                                confidence.intervals = 0.95,study.design = 'cohort',)

model_DCA <- decision_curve(formula = formula,data = train_merged_data,
                            family = binomial(link = 'logit'),
                            confidence.intervals = 0.95,study.design = 'cohort',)


model_all <- list(model_clinic_DCA,model_miR_DCA,model_DCA)

plot_decision_curve(model_all,
                    curve.names = c('model_Clinical Factor','model_miR','model_Combined'),
                    xlim = c(0,0.8),
                    cost.benefit.axis = F,
                    col = c("#377DB8",  "#4DAE49",  "#E31A1C"),
                    confidence.intervals = F,
                    lwd = 3.5,
                    standardize = T)



##临床影响曲线
model_DCA <- decision_curve(formula = formula,data = train_merged_data,
                            family = binomial(link = 'logit'),
                            confidence.intervals = 0.95,study.design = 'cohort',)

plot_clinical_impact(model_DCA,population.size = 97,cost.benefit.axis = T,
                     n.cost.benefits = 8,col = c("#BC3C29FF","#3A5FCD"),
                     confidence.intervals = T)



####5.combined模型可视化#####
 ####nomogram图绘制：变量名注意不含特殊字符
#动态列线图
library(DynNom)
dat=train_merged_data


colnames(dat) <- gsub("V1.", "pre.", colnames(dat))
colnames(dat) <- gsub("V2.", "post.", colnames(dat))

#将dat中Gender这一列的1替换为字符male，2替换为字符female
#dat$Gender <- ifelse(dat$Gender == 1, "male", ifelse(dat$Gender == 2, "female", dat$Gender))
                       
colnames(dat) <- gsub("pre.miR.502.3p", "pre.hsa.miR.502.3p", colnames(dat))
colnames(dat) <- gsub("pre.HBA1C", "pre.HbA1c", colnames(dat))
colnames(dat) <- gsub("post.ISI", "post.ISSI.2", colnames(dat))


formular2 <- Remission ~ `pre.UA` +`pre.Weight` +`post.ISSI.2`+`pre.BMI`+`Age`+`Gender`+`pre.HbA1c`+ `pre.hsa.miR.502.3p`
model <- glm(formular2,
             data =dat, family = binomial())

DynNom(model,DNtitle="Nomogram",DNxlab="Probability")

DNbuilder(model) 

rsconnect::setAccountInfo(name='siit1yearrediction-tangswandlhhj',
                          token='6EC6ECCF3754151B14C2E1F8CE8C6349',
                          secret='IAV8M/8mTgomBUAHphY8xXznlhxWlU1a7T0hs5mF')
library(rsconnect)
rsconnect::deployApp('path/to/your/app')



#静态列线图
library(rms)
dd <- datadist(train_merged_data)
options(datadist = "dd")

model_rms <- lrm(Remission ~ `V1.UA` + `V1.Weight` + `V2.ISI` + V1.miR.502.3p, 
                 data = train_merged_data)

model_rms <- model
library(regplot)
nom.glm <- regplot(model_rms, 
                   observation = train_merged_data[1,], 
                   center = TRUE, 
                   title = "Nomogram",
                   points = TRUE, 
                   odds = FALSE, 
                   showP = TRUE, 
                   rank = "sd", 
                   clickable = FALSE)


######7.模型效能验证#######
####7.1抽样验证区分度bootstrap####
library(pROC)
train_merged_data$Remission <- as.factor(train_merged_data$Remission)
levels(train_merged_data$Remission) <- make.names(levels(train_merged_data$Remission))
levels(test_merged_data$Remission) <- make.names(levels(test_merged_data$Remission))
str(train_merged_data$Remission)
str(test_merged_data$Remission)

#traindata
form.bestglm<-as.formula(Remission ~`V1.UA`+`V1.Weight` +`V2.ISI`+ V1.miR.502.3p )
train.control <- trainControl(method = "boot", number = 1000)

LogMod <- train(form.bestglm, 
                data=train_merged_data, 
                trControl=train.control, 
                method='glm')

predicted_probs <- predict(LogMod, newdata = train_merged_data, type = "prob")
pred <- predicted_probs$X1
labels <- train_merged_data$Remission

n_bootstraps <- 1000
roc_boot <- NULL

#重抽样
for (i in 1:n_bootstraps) {
  boot_indices <- sample(length(pred), replace = TRUE)
  boot_pred <- pred[boot_indices]
  boot_labels <- labels[boot_indices]
  roc_boot[[i]] <- roc(boot_labels, boot_pred)}

#平均ROC
roc_mean <- roc(labels,pred,ci=TRUE)

#每一次ROC
plot(roc_boot[[1]], type = "n", main = "Bootstrap ROC Curve", 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

for (i in 1:n_bootstraps) {
  lines(roc_boot[[i]], col = "grey", alpha = 0.2)}

lines(roc_mean, col = "chartreuse4", lwd = 3)

auc_mean <- auc(roc_mean)
legend_text <- c("Bootstrap ROC", paste("Mean ROC (AUC =", round(auc_mean, 3), ")"))
legend("bottomright", legend = legend_text, 
       col = c("grey", "chartreuse4"), lwd = c(3, 3), bty = "n",cex = 1)


#testdata
LogMod <- train(form.bestglm, 
                data=test_merged_data, 
                trControl=train.control, 
                method='glm')
predicted_probs <- predict(LogMod, newdata = test_merged_data, type = "prob")
pred <- predicted_probs$X1
labels <- test_merged_data$Remission

n_bootstraps <- 1000
roc_boot <- NULL

#重抽样
for (i in 1:n_bootstraps) {
  boot_indices <- sample(length(pred), replace = TRUE)
  boot_pred <- pred[boot_indices]
  boot_labels <- labels[boot_indices]
  roc_boot[[i]] <- roc(boot_labels, boot_pred)}

#平均ROC
roc_mean <- roc(labels,pred,ci=TRUE)

#每一次ROC
plot(roc_boot[[1]], type = "n", main = "Bootstrap ROC Curve", 
     xlab = "False Positive Rate", ylab = "True Positive Rate")

for (i in 1:n_bootstraps) {
  lines(roc_boot[[i]], col = "grey", alpha = 0.2)}

lines(roc_mean, col = "chartreuse4", lwd = 3)

auc_mean <- auc(roc_mean)
legend_text <- c("Bootstrap ROC", paste("Mean ROC (AUC =", round(auc_mean, 3), ")"))
legend("bottomright", legend = legend_text, 
       col = c("grey", "chartreuse4"), lwd = c(3, 3), bty = "n")


####7.2.抽样验证校准度####
library(rms)
dd=datadist(train_merged_data)
options(datadist="dd")

library(riskRegression)

xb_train <- Score(list("model_clinic"=model_clinic,
                 "model_miR"=model_miR,
                 "model"=model),
            formula = Remission~1,
            null.model = F,
            conf.int = T,
            plots = c("calibration","ROC"),
            metrics = c("auc","brier"),
            B=1000,M=5,
            data=train_merged_data)

plotCalibration(xb_train,cex = 0.8)

xb_test <- Score(list("model_clinic"=model_clinic,
                 "model_miR"=model_miR,
                 "model"=model),
            formula = Remission~1,
            null.model = F,
            conf.int = T,
            plots = c("calibration","ROC"),
            metrics = c("auc","brier"),
            B=1000,M=5,
            data=test_merged_data)

plotCalibration(xb_test,cex = 0.8)


#loess平滑曲线
library(ggplot2)
library(rms)
library(runway)

y <- test_merged_data[, "Remission"]

pr_clinic <- predict(model_clinic, newdata = test_merged_data, type = c("response"))
pr_miR <- predict(model_miR, newdata = test_merged_data, type = c("response"))
pr <- predict(model, newdata = test_merged_data, type = c("response"))


model_name <- rep(c("model_clinical"), length.out = 28)
df_c <- data.frame(
  outcomes = y,
  predictions = pr_clinic,
  model_name = model_name
)

model_name <- rep(c("model_miR"), length.out = 28)
df_m <- data.frame(
  outcomes = y,
  predictions = pr_miR,
  model_name = model_name
)

model_name <- rep(c("model"), length.out = 28)
df <- data.frame(
  outcomes = y,
  predictions = pr,
  model_name = model_name
)

multi_model_dataset <- rbind(df_c,df_m,df)

multi_model_dataset$outcomes <- ifelse(multi_model_dataset$outcomes == "X0", 0, 
                                       ifelse(multi_model_dataset$outcomes == "X1", 1, 
                                              multi_model_dataset$outcomes))

threshperf_plot_multi(multi_model_dataset,
                      outcome = 'outcomes',
                      prediction = 'predictions',
                      model = 'model_name',
                      positive = 1)


cal_plot_multi(multi_model_dataset,
               outcome = 'outcomes',
               prediction = 'predictions',
               model = 'model_name',
               n_bins = 0,
               show_loess = TRUE,
               positive = 1)



####8.靶器官寻找 相关性计算####
library(correlation)
library(dplyr)
library(see)
library(readxl)
library(writexl)

X97_baseline_临床信息 <- read_excel("PCR/97-baseline-临床信息.xlsx")
X97_follow_up_临床信息 <- read_excel("PCR/97-follow-up-临床信息.xlsx")
X97_delta_临床信息 <- read_excel("PCR/97-delta-临床信息.xlsx")

data <- X97_baseline_临床信息[,-c(1,2)]
data <- X97_follow_up_临床信息[,-c(1,2)]
data <- X97_delta_临床信息[,-c(1,2)]

data <- na.omit(data)
str(data)

data <- data %>%
  mutate(across(everything(), as.numeric))

check_normality_and_correlate <- function(var1, var2) {
  normal_test_var1 <- shapiro.test(var1)
  normal_test_var2 <- shapiro.test(var2)
  
  # Initialize the variable to store the test type
  test_type <- NA
  
  # If both variables are normally distributed, use Pearson's correlation
  if (normal_test_var1$p.value > 0.05 && normal_test_var2$p.value > 0.05) {
    correlation_result <- cor.test(var1, var2, method = "pearson")
    test_type <- "Pearson"  # Indicate that Pearson's test was used
  } else {
    # Otherwise, use Spearman's correlation
    correlation_result <- cor.test(var1, var2, method = "spearman")
    test_type <- "Spearman"  # Indicate that Spearman's test was used
  }
  
  # Return a data frame with correlation, p-value, and test type
  return(data.frame(
    correlation = correlation_result$estimate,
    p_value = correlation_result$p.value,
    test_used = test_type  # Add the test type to the data frame
  ))
}
results <- data.frame()

for (i in 1:(ncol(data) - 1)) {
  for (j in (i + 1):ncol(data)) {
    var1 <- data[[i]]
    var2 <- data[[j]]
    correlation_result <- check_normality_and_correlate(var1, var2)
    
    # Store the results in the data frame, including the test used
    results <- rbind(results, data.frame(
      variable1 = colnames(data)[i],
      variable2 = colnames(data)[j],
      correlation = correlation_result$correlation,
      p_value = correlation_result$p_value,
      test_used = correlation_result$test_used  # Add the test type to the results
    ))
  }
}
rownames(results) <- NULL

#selected_vars <- c("V1-miR-502-3p")#基线的只需要做基线gene与基线临床指标

selected_vars <- c("V1-miR-502-3p", "V2-miR-502-3p")
selected_vars <- c("V1-miR-502-3p", "V2-miR-502-3p","delta")

filtered_results <- results %>%
  filter(variable1 %in% selected_vars | variable2 %in% selected_vars)

write_xlsx(filtered_results, "filtered_results.xlsx")

library(ggplot2)
library(ggpubr)
library(ggpmisc)

colnames(data) <- gsub("-", ".", colnames(data))
df <- data

#AUCglu和基线miR-502-3p的相关性
ggscatter(df, 
          x = "V1.miR.502.3p", 
          y = "AUCglu",
          add = "reg.line", 
          conf.int = TRUE,    
          add.params = list(fill = "lightgray", color = "#006400"),  # 设置回归线的填充和颜色
          color = "#006400",  # 设置点的颜色
          size = 2.5,  # 设置点的大小，调整数值可以改变点的大小
          alpha = 0.8  # 设置点的透明度，范围 0 到 1
) +
  stat_cor(method = "pearson", 
           label.x = 2, 
           label.y = 30) +
  labs(x = "miR-502-3p baseline",
       y = "AUC.glu  after  SIIT") +
  theme_pubr() +  # 使用 ggpubr 提供的主题
  theme(legend.position = "top",  # 设置图例位置
        panel.grid.major = element_line(color = "lightgray", size = 0.5),  # 大网格线
        panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # 小网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        plot.background = element_rect(fill = "white"),  # 整体背景颜色 
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # 设置边框
)

#ISSI2和基线miR-502-3p的相关性
ggscatter(df, 
          x = "V1.miR.502.3p", 
          y = "ISSI.2",
          add = "reg.line", 
          conf.int = TRUE,    
          add.params = list(fill = "lightgray", color = "#006400"),  # 设置回归线的填充和颜色
          color = "#006400",  # 设置点的颜色
          size = 2.5,  # 设置点的大小，调整数值可以改变点的大小
          alpha = 0.8  # 设置点的透明度，范围 0 到 1
) +
  stat_cor(method = "spearman", 
           label.x = 2, 
           label.y = 37) +
  labs(x = "miR-502-3p baseline",
       y = "ISSI-2 after SIIT") +
  theme_pubr() +  # 使用 ggpubr 提供的主题
  theme(legend.position = "top",  # 设置图例位置
        panel.grid.major = element_line(color = "lightgray", size = 0.5),  # 大网格线
        panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # 小网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        plot.background = element_rect(fill = "white"),  # 整体背景颜色 
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # 设置边框
)

#DI-60和基线miR-502-3p的相关性
ggscatter(df, 
          x = "V1.miR.502.3p", 
          y = "DI.60",
          add = "reg.line", 
          conf.int = TRUE,    
          add.params = list(fill = "lightgray", color = "#006400"),  # 设置回归线的填充和颜色
          color = "#006400",  # 设置点的颜色
          size = 2.5,  # 设置点的大小，调整数值可以改变点的大小
          alpha = 0.8  # 设置点的透明度，范围 0 到 1
) +
  stat_cor(method = "spearman", 
           label.x = 2, 
           label.y = 60) +
  labs(x = "miR-502-3p baseline",
       y = "DI after SIIT") +
  theme_pubr() +  # 使用 ggpubr 提供的主题
  theme(legend.position = "top",  # 设置图例位置
        panel.grid.major = element_line(color = "lightgray", size = 0.5),  # 大网格线
        panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # 小网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        plot.background = element_rect(fill = "white"),  # 整体背景颜色 
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # 设置边框
)

#体重和基线miR-502-3p的相关性
ggscatter(df, 
          x = "V1.miR.502.3p", 
          y = "体重",
          add = "reg.line", 
          conf.int = TRUE,    
          add.params = list(fill = "lightgray", color = "#006400"),  # 设置回归线的填充和颜色
          color = "#006400",  # 设置点的颜色
          size = 2.5,  # 设置点的大小，调整数值可以改变点的大小
          alpha = 0.8  # 设置点的透明度，范围 0 到 1
) +
  stat_cor(method = "spearman", 
           label.x = 2, 
           label.y = 10) +
  labs(x = "miR-502-3p baseline",
       y = "Weight change") +
  theme_pubr() +  # 使用 ggpubr 提供的主题
  theme(legend.position = "top",  # 设置图例位置
        panel.grid.major = element_line(color = "lightgray", size = 0.5),  # 大网格线
        panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # 小网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        plot.background = element_rect(fill = "white"),  # 整体背景颜色 
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # 设置边框
)

#V2-ins60和V2-miR-502-3p的相关性
ggscatter(df, 
          x = "V2.miR.502.3p", 
          y = "ins60",
          add = "reg.line", 
          conf.int = TRUE,    
          add.params = list(fill = "lightgray", color = "#006400"),  # 设置回归线的填充和颜色
          color = "#006400",  # 设置点的颜色
          size = 2.5,  # 设置点的大小，调整数值可以改变点的大小
          alpha = 0.8  # 设置点的透明度，范围 0 到 1
) +
  stat_cor(method = "spearman", 
           label.x = -6.5, 
           label.y = 0.2) +
  labs(x = "miR-502-3p follow-up",
       y = "ins60") +
  theme_pubr() +  # 使用 ggpubr 提供的主题
  theme(legend.position = "top",  # 设置图例位置
        panel.grid.major = element_line(color = "lightgray", size = 0.5),  # 大网格线
        panel.grid.minor = element_line(color = "lightgray", size = 0.25),  # 小网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        plot.background = element_rect(fill = "white"),  # 整体背景颜色 
        panel.border = element_rect(color = "black", fill = NA, size = 1)  # 设置边框
)



#####9.配对箱线图######
library(readxl)
PCR_result_5 <- read_excel("PCR/PCR-result-5.xlsx")
df <- PCR_result_5
df$Group <- factor(df$Group, levels = c("pre", "post"))

library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggsci)
library(ggbeeswarm)
library(gghalves)

df %>%
  ggplot(aes(Group, Expression)) +  
  geom_half_violin(aes(split = Group), side = 2, alpha = 0.8, size = 1,color="grey") +  
  stat_boxplot(geom = "errorbar", width = 0.2, size = 1,color="grey") +
  geom_boxplot(width = 0.4, size = 0.95,color="grey") +
  geom_line(aes(group = paired), color = "grey80") +
  geom_point(aes(fill = Group, group = paired, size = Expression, alpha = Expression), 
             pch = 21, position = position_dodge(0.2), size = 3.5, color = "black") +  # 黑色边框
  
  scale_fill_manual(values = c( "#a48ad3", "#7ca0d4")) +  # 自定义颜色
  scale_size_continuous(range = c(1, 3)) +
  
  geom_signif(comparisons = list(c("pre", "post")),
              map_signif_level = FALSE, vjust = 0.1, color = "black",
              textsize = 4, test = t.test, step_increase = 0.1) +
  
  facet_wrap(. ~ Type, nrow = 1) +
  scale_y_continuous(limits = c(-9, 9), minor_breaks = seq(-9, 9, 4)) +
  labs(x = NULL, y = NULL) +
  theme_test() +
  theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        axis.line = element_line(color = "black", size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2, color = "#e5e5e5"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(margin = margin(t = 2), color = "black", size = 12),
        strip.text = element_text(size = 11.5),  # 放大 Type 标签
        legend.position = "none",
        panel.spacing = unit(0.5, "lines")) +
  coord_cartesian()
