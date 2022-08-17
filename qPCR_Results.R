## Script for converting Quant Studio 3 qPCR txt output (results tab) into results for managers

## Install required R packages
##install.packages("dplyr")
##install.packages("tidyr")
##install.packages("ggplot2")
library("dplyr")
library("tidyr")
library("ggplot2")

#Import data 
Results_all <- read.delim("zq_test_results.csv", sep = ",")
#select required columns
df <- Results_all[,c(3:6, 9, 11, 12, 18, 24)]
# Remove undetermined Ct 
df$CT <- gsub("Undetermined",0, df$CT)

#Use Filter to create data frame with just target results and 1 with just IPC results
IPC_Res <- filter(df, Target.Name == "IPC" & Task == "UNKNOWN" & Omit == "FALSE")
IPC_sample_sum <- group_by(IPC_Res, Sample.Name) %>%
  summarize(Ct.mean = mean(as.numeric(CT)), Ct.SD = mean(Ct.SD))
#Change Column name so IPC Ct mean is distinguishable from Target Ct mean
colnames(IPC_sample_sum) <- c('Sample.Name', 'IPC.Ct', 'IPC.Ct.SD')

## Calculate Ct Delays (IPC Mean Ct - IPC_POS Mean Ct)
#Get IPC POS Ct value and repeat for length of IPC_sample_sum
IPC_POS_delay <- filter(IPC_sample_sum, Sample.Name == "IPC POS")
IPC <- IPC_POS_delay[,2]
IPC_length = nrow(IPC_sample_sum)
IPC_delays <- rep(IPC, each = IPC_length)
Ct_delays = IPC_sample_sum$IPC.Ct - IPC_delays$IPC.Ct
my_Ct_delay = data.frame(Ct_delays)
Ct_delay = my_Ct_delay[-1,]

# Use Filter to create data frame with just target results and 1 with just IPC results 
Target_Res <- filter(df, Target.Name != "IPC" & Task == "UNKNOWN" & Omit == "FALSE" & Sample.Name !="IPC POS")
Target_Res$Quantity <-as.numeric(as.character(Target_Res$Quantity))

#Group by Sample.Name to get summary stats for qPCR replicates by sample
Target_sample_sum <- group_by(Target_Res, Sample.Name) %>%
  summarize(Sample.N = n(), Amp.N = sum(Amp.Status=="Amp"), 
            Copy.mean = mean(Quantity))
Target_df <- Target_sample_sum

df_2 <- cbind(Target_df, Ct_delay)

Target_length = nrow(Target_df)
## get IPC POS Ct SD
IPC_Ct_SD <- data.frame(IPC_POS_delay[,3])
IPC_SD <- rep(IPC_Ct_SD, each = Target_length)

## get efficiency coefficient for adjusting quantity
Efficiency_coefficient <- data.frame(rep((Target_Res[1,8]/100)*2, each = Target_length))

## calculate adjusted quantity
df_3 <- cbind(df_2, IPC_Ct_SD, Efficiency_coefficient)
df_3[is.na(df_3)] <- 0
colnames(df_3) <- c('Sample.Name', 'Sample.N', 'Amp.N','Copy.mean', 
                    'Ct.Delay', 'IPC.Ct.SD', 'Efficiency Coefficient')
adj_Ct = (df_3$Ct.Delay - df_3$IPC.Ct.SD)

df_4 <- cbind(df_3[,-6], adj_Ct)

adj_quant <- (df_4$adj_Ct*df_4$Copy.mean*df_4$`Efficiency Coefficient`)+ df_4$Copy.mean
## combine with d4 (without adj_Ct and Efficiency coefficient columns)
my_results <- cbind(df_4[,-c(6,7)], adj_quant)  

# Add Yes if Copy.mean is above LOD (10) and Amp.N is greater than 1
target_detection = ifelse(test= my_results$adj_quant > 10 & my_results$Amp.N > 1,
                          yes = "Yes",
                          no = "No")
## create final result table and export CSV
result_table <- cbind(my_results, target_detection)

write.csv(result_table, "C:/Users/qPCR_results_test.csv")





