#############################################

#Change the path to your output files, which will have names specific to the model you labeled
#These are for the example output files provided:
sim_data <- read.delim("/Goodness_of_Fit/Empirical/Simulation_Results.txt", header = TRUE, sep = "\t")
emp_data <- read.delim("/Goodness_of_Fit/Empirical/Empirical.sym_mig.optimized.txt", header = TRUE, sep = "\t")

#check headers
ls(sim_data)
#should be: Simulation	Best_Replicate	log-likelihood	theta	sfs_sum	chi-squared	optimized_params
ls(emp_data)
#should be: Model	Replicate	log-likelihood	theta	sfs_sum	chi-squared

#quick summary of files:
summary(sim_data)
summary(emp_data)

#summary of variables for simulated sfs fits:
summary(sim_data$log.likelihood)
summary(sim_data$chi.squared)

#summary of empirical values:
emp_data$log.likelihood
emp_data$chi.squared


#####################################
#FREQUENCY DISTRIBUTIONS
#Edit below numbers to change the frequency distribution set up for your data,
#by checking the summary functions above.

#syntax for setting bin number is: seq(low number, high number, increment)
#-----------------------------------------------------------------------------------------

#likelihood plot
#change the bin settings below
ll_seq<- seq(-100,-20,2)
hist(sim_data$log.likelihood, breaks=ll_seq, main = "Simulation Results - Log-likelihood distribution", xlab="log-likelihood", col="grey")
abline(v=emp_data$log.likelihood, lwd = 3, col='blue')

#chi-squared plot
#change the bin settings below
chi_seq<- seq(0,100,5)
hist(sim_data$chi.squared, breaks=chi_seq, main = "Simulation Results - Chi-squared distribution ", xlab="Chi-squared test statistic", col="grey")
abline(v=emp_data$chi.squared, lwd = 3, col='blue')

#####################################
#For log-transformed chi.squared

#transform chi-squared test stat
log_chi <- log(sim_data$chi.squared)
emp_log_chi <- log(emp_data$chi.squared)
summary(log_chi)


#log transformed chi-squared plot
#change the bin settings below
lchi_seq<- seq(0,8,0.2)
hist(log_chi, breaks=lchi_seq, main = "Simulation Results - Log Chi-squared distribution ", xlab="log Chi-squared test statistic", col="grey")
abline(v=emp_log_chi, lwd = 3, col='blue')



