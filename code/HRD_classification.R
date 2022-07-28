##upload models
output_AM = readRDS("output_AM.rds")
output_TAI= readRDS("output_TAI.rds")
output_LST= readRDS("output_LST.rds")
output_HRD= readRDS("output_HRD.rds")


##Paired run

#Import ASCAT output data for classification

test_comb <- read.delim("all_paired_tum_jun2019.tsv", header = T, sep = "\t") ##paired_input result
rownames(test_comb) <- test_comb$sample_id
test.data.comb <- test_comb[,c(1,13:14,11:12)]
test.data.comb$Arithmetic.mean <- rowMeans(subset(test.data.comb, select = c("Telomeric.AI", "LST", "HRD")), na.rm = TRUE)

#check model on ILMN data: predictions
P.comb.test.AM <- round(predict(output_AM,newdata = data.frame(Arithmetic.mean = test.data.comb$Arithmetic.mean), type = "response"))
P.comb.test.TAI <- round(predict(output_TAI,newdata = data.frame(Telomeric.AI.1 = test.data.comb$Telomeric.AI), type = "response"))
P.comb.test.LST <- round(predict(output_LST,newdata = data.frame(LST.1 = test.data.comb$LST), type = "response"))
P.comb.test.HRD <- round(predict(output_HRD,newdata = data.frame(HRD.LOH.1 = test.data.comb$HRD), type = "response"))

#save as table
df_pred_comb <- cbind(test.data.comb, P.comb.test.AM, P.comb.test.TAI, P.comb.test.LST, P.comb.test.HRD)
#kable(df_pred_comb, format = "markdown")
df_pred_comb$ID = rownames(df_pred_comb)
write.table(df_pred_comb, file="pred_res_paired.tsv", row.names = F, sep = "\t", quote = F)


##Single sample run

#Import ASCAT output data for classification

test_comb <- read.delim("all_unpaired_tum_Jun2019.tsv", header = T, sep = "\t") ##paired_input result
rownames(test_comb) <- test_comb$sample_id
test.data.comb <- test_comb[,c(1,13:14,11:12)]
test.data.comb$Arithmetic.mean <- rowMeans(subset(test.data.comb, select = c("Telomeric.AI", "LST", "HRD")), na.rm = TRUE)

#check model on ILMN data: predictions
P.comb.test.AM <- round(predict(output_AM,newdata = data.frame(Arithmetic.mean = test.data.comb$Arithmetic.mean), type = "response"))
P.comb.test.TAI <- round(predict(output_TAI,newdata = data.frame(Telomeric.AI.1 = test.data.comb$Telomeric.AI), type = "response"))
P.comb.test.LST <- round(predict(output_LST,newdata = data.frame(LST.1 = test.data.comb$LST), type = "response"))
P.comb.test.HRD <- round(predict(output_HRD,newdata = data.frame(HRD.LOH.1 = test.data.comb$HRD), type = "response"))

#save as table
df_pred_comb <- cbind(test.data.comb, P.comb.test.AM, P.comb.test.TAI, P.comb.test.LST, P.comb.test.HRD)
df_pred_comb$ID = rownames(df_pred_comb)
#kable(df_pred_comb, format = "markdown")
write.table(df_pred_comb, file="pred_res_unpaired.tsv", row.names = F, sep = "\t", quote = F)


#Check concordancy of single sample runs with paired samples.

