# Executed 4/__/2025
# Set working directory
# No packages required
# setwd("~/Desktop")

# Read in raw counts from csv file for a single recombination experiment 
# Data table with structure "Class, Interval_1, Interval_2, ... , Count"
# Recombination pattern for each class encoded as 1,0 in interval column
#observed_dataset <- read.csv("../Raw Data/Pilots/Recombination Rate/counting_model_dataset.csv")

# Extract counts and calculate total
observed_count <- observed_dataset$Count
observed_total <- sum(observed_count)

# Define number of marker based intervals and max p&k value in an interval
markers <- 5
intervals <- 4
max_p_value <- 10
max_k_value <- 3

# Enter estimate of genetic length (x) in each interval (index) in Morgans
# Use to calculate initial value of model parameter y optimized by optim()
# Can use published genetic maps or recombination fractions from raw data
x1 <- 0.12
x2 <- 0.21
x3 <- 0.25
x4 <- 0.06
genetic_lengths <- c(x1,x2,x3,x4)

# Generate combinations of the model parameter p in marker-based intervals 
# Parameter p gives strength of renewal process (p=m+1) in the CxCoM model
# Evaluate values of p from "no interference" (p=1) up to (p=max_p_value)
# p_uniform is matrix of all possible p's +/- 1 from previous estimate p
for (l in 1:intervals) {
  assign(paste0("p", l), c(1:max_p_value))
}
p_uniform <- cbind(p1,p2,p3,p4)
colnames(p_uniform) <- c("p1","p2","p3","p4")

# Generate starting point for parameter y in each interval to be optimized
# Calculate as 2*p*x using the p_uniform matrix and genetic_length vector
y_uniform <- as.matrix(cbind(2*p_uniform[,1]*x1,
                             2*p_uniform[,2]*x2,
                             2*p_uniform[,3]*x3,
                             2*p_uniform[,4]*x4))
colnames(y_uniform) <- c("y1","y2","y3","y4")

# Generate combinations of the viability model parameter v for each marker 
# Parameter v gives relative viability of marker homozygote to heterozygote
# v_uniform is matrix of all v in order from centromere to telomere
v_uniform <- matrix(1,max_p_value,markers)
colnames(v_uniform) <- c("v1","v2","v3","v4","v5")

# Custom function to calculate the likelihood of a specific CxCoM model
# Corresponds to the analysis of chi-square model from Zhao et al. 1995
# Required inputs are "parameters" a vector of y and v values to be optimized
# Probabilities for classes will need to be tailored to each experiment
Uniform.CxCoM.Likelihood <- function(parameters) {
  
  # Vectors for internal use in the custom function
  # Dimensions auto-adjusted for combinations of p
  y <- parameters[1:4]
  v <- parameters[5:9]
  stationary_vector <- matrix(data=1, nrow=1, ncol=p[1])
  summation_vector <- matrix(data=1, nrow=p[intervals], ncol=1)
  
  # Matrices for internal use in custom function, stored as lists
  # Dimensions auto-adjusted to allow valid matrix multiplication
  kclass <- list(matrix(data=1, nrow=ncol(stationary_vector), ncol=p[1]))
  for (l in 2:intervals) {
    kclass[[l]] <- matrix(data=1, nrow=p[l-1], ncol=p[l])
  }
  
  k0class <- list(matrix(data=1, nrow=ncol(stationary_vector), ncol=p[1]))
  for (l in 2:intervals) {
    k0class[[l]] <- matrix(data=1, nrow=p[l-1], ncol=p[l])
  }
  
  # CASE OF k > 0
  #-------------------------------------------------------------------------
  # For loops giving probability of "k" number of crossovers in interval "l"
  #   1) assuming "s" is the number of "C" precursor events in interval "l" 
  #      and events are produced by poisson process with rate parameter "y"
  #   2) assuming "C" precursor events mature to "k" number of "Cx" events
  #      with "s-k" number of "Co" events in the counting process "Cx(Co)^m"
  #   3) assuming the "Cx(Co)^m" counting process is stationary in that the
  #      first "C" equally likely to mature into "Cx" or any one of the "Co"
  #   4) defining "i" as the identity of the first "C" in interval "l" that
  #      numerically is position in series CxCoCoCo plus 1 (i for Cx = 1)
  #   5) defining "j" as the identity of the last "C" in interval "l" that
  #      numerically is p minus position in series CoCoCoCx (j for Cx = 0)
  #
  # Check internal inconsistency of Zhao et al. 1995 main text and appendix
  # R script uses pk-p+i+j definition from Zhao et al. 1995 appendix proof
  #
  # Dk_yl is p(l-1) by p(l) matrix of crossover probabilities in interval l
  # Summing over all k>0 (max_k_value=3) and assuming s>0 in the interval l
  # Probability of k=0 treated separately
  #
  # Matrix multiplication using interval and parameter notation requires 
  # First interval treated separately for valid matrix and multiplication
  # First interval is square matrix by virtue of stationary assumption but
  # p[l-1] not defined (would have been number of columns in 0th interval)
  
  # Nested for loops for the first interval
  interval_matrix <- as.matrix(kclass[[1]])
    for (k in 1:max_k_value) {
      for (i in 1:p[1]) {
        for (j in (p[1]-1):0) {
          interval_matrix[i,(p[1]-j)] <- ((exp(-y[1]))*(y[1]^((p[1]*k)-p[1]+i+j)))/(factorial((p[1]*k)-p[1]+i+j))
        } 
      }
      assign(paste0("Dk", k), interval_matrix)
    }
  assign(paste0("Dk_y", 1), Dk1+Dk2+Dk3)
  
  # Nested for loops for all other intervals
  for (l in 2:intervals) {
    interval_matrix <- kclass[[l]]
      for (k in 1:max_k_value) {
        for (i in 1:p[l-1]) {
          for (j in (p[l]-1):0) {
            interval_matrix[i,(p[l]-j)] <- ((exp(-y[l]))*(y[l]^((p[l]*k)-p[l]+i+j)))/(factorial((p[l]*k)-p[l]+i+j))
          } 
        }
        assign(paste0("Dk", k), interval_matrix)
      }
    assign(paste0("Dk_y", l), Dk1+Dk2+Dk3)
  }
  
  # R_intervals is a p(l-1) by p(l) matrix describing the probability of
  # observing a crossover in interval l with model parameters yl and pl
  R_1 <- (1/2)*Dk_y1
  R_2 <- (1/2)*Dk_y2
  R_3 <- (1/2)*Dk_y3
  R_4 <- (1/2)*Dk_y4
  
  # CASE OF k = 0
  #-------------------------------------------------------------------------
  # For loops giving probability of "k" number of crossovers in interval "l"
  #   1) assuming "s" is the number of "C" precursor events in interval "l" 
  #      and events are produced by poisson process with rate parameter "y"
  #   2) assuming "C" precursor events mature to "k" number of "Cx" events
  #      with "s-k" number of "Co" events in the counting process "Cx(Co)^m"
  #   3) assuming the "Cx(Co)^m" counting process is stationary in that the
  #      first "C" equally likely to mature into "Cx" or any one of the "Co"
  #   4) defining "i" as the identity of the first "C" in interval "l" that
  #      numerically is position in series CxCoCoCo plus 1 (i for Cx = 1)
  #   5) defining "j" as the identity of the last "C" in interval "l" that
  #      numerically is p minus position in series CoCoCoCx (j for Cx = 0)
  #
  # Check internal inconsistency of Zhao et al. 1995 main text and appendix
  # R script uses pk-p+i+j definition from Zhao et al. 1995 appendix proof
  #
  # Lower triangular matrix D0 in Zhao et al. 1995 was not defined anywhere
  # R script uses implied definition in Zhao et al. 1995 appendix theorem 1
  # D0_yl is p(l-1) by p(l) matrix of probability no crossover in interval l
  #
  # Matrix multiplication using interval and parameter notation requires 
  # First interval treated separately for valid matrix and multiplication
  # First interval is square matrix by virtue of stationary assumption but
  # p[l-1] not defined (would have been number of columns in 0th interval)
  
  # Nested for loops for the first interval
  interval_matrix <- as.matrix(k0class[[1]])
  for (i in 1:p[1]) {
    for (j in (p[1]-1):0) {
      if ((0-p[1]+i+j)>-1) {
        interval_matrix[i,(p[1]-j)] <- ((exp(-y[1]))*(y[1]^(0-p[1]+i+j)))/(factorial(0-p[1]+i+j))
      } 
      else {
        interval_matrix[i,(p[1]-j)]<-0 
      } 
    } 
  }
  assign(paste0("D0_y", 1), interval_matrix)
  
  # Nested for loops for all other intervals
  for (l in 2:intervals) {
    interval_matrix <- as.matrix(k0class[[l]])
      for (i in 1:p[l-1]) {
        for (j in (p[l]-1):0) {
          if ((0-p[l]+i+j)>-1) {
            interval_matrix[i,p[l]-j] <- ((exp(-y[l]))*(y[l]^(0-p[l]+i+j)))/(factorial(0-p[l]+i+j))
          } 
          else {
            interval_matrix[i,(p[l]-j)]<-0 
          } 
        } 
      }
    assign(paste0("D0_y", l), interval_matrix)
  }
  
  # N_intervals is a p(l-1) by p(l) matrix describing the probability of
  # Not observing a crossover in interval l with the parameters yl and pl
  N_1 <- D0_y1+((1/2)*Dk_y1)
  N_2 <- D0_y2+((1/2)*Dk_y2)
  N_3 <- D0_y3+((1/2)*Dk_y3)
  N_4 <- D0_y4+((1/2)*Dk_y4)
  
  # Probability of observing no crossover events
  Pr11111<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%summation_vector)*((v[1]*v[2]*v[3]*v[4]*v[5])/2) #Pr(y,cv,v,f,car)
  Pr00000<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%N_4%*%summation_vector)*((1.00*1.00*1.00*1.00*1.00)/2) #Pr(+,+,+,+,+)
      
  # Probability of observing only a single crossover in interval y-cv
  Pr00001<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%summation_vector)*((1.00*1.00*1.00*1.00*v[5])/2) #Pr(y,+,+,+,+)
  Pr11110<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%N_3%*%R_4%*%summation_vector)*((v[1]*v[2]*v[3]*v[4]*1.00)/2) #Pr(+,cv,v,f,car)
  # Probability of observing only a single crossover in interval cv-v
  Pr00011<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%summation_vector)*((1.00*1.00*1.00*v[4]*v[5])/2) #Pr(y,cv,+,+,+)
  Pr11100<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%N_4%*%summation_vector)*((v[1]*v[2]*v[3]*1.00*1.00)/2) #Pr(+,+,v,f,car)
  # Probability of observing only a single crossover in interval v-f
  Pr00111<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%summation_vector)*((1.00*1.00*v[3]*v[4]*v[5])/2) #Pr(y,cv,v,+,+)
  Pr11000<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%N_4%*%summation_vector)*((v[1]*v[2]*1.00*1.00*1.00)/2) #Pr(+,+,+,f,car)
  # Probability of observing only a single crossover in interval f-car
  Pr01111<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%summation_vector)*((1.00*v[2]*v[3]*v[4]*v[5])/2) #Pr(y,cv,v,f,+)
  Pr10000<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%N_4%*%summation_vector)*((v[1]*1.00*1.00*1.00*1.00)/2) #Pr(+,+,+,+,car)
  
  # Probability of observing double crossover in intervals y-cv and cv-v
  Pr11101<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%summation_vector)*((v[1]*v[2]*v[3]*1.00*v[5])/2) #Pr(y,+,v,f,car)
  Pr00010<-((1/p[1])%*%stationary_vector%*%N_1%*%N_2%*%R_3%*%R_4%*%summation_vector)*((1.00*1.00*1.00*v[4]*1.00)/2) #Pr(+,cv,+,+,+)
  # Probability of observing double crossover in intervals y-cv and v-f
  Pr11001<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%summation_vector)*((v[1]*v[2]*1.00*1.00*v[5])/2) #Pr(y,+,+,f,car)
  Pr00110<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%N_3%*%R_4%*%summation_vector)*((1.00*1.00*v[3]*v[4]*1.00)/2) #Pr(+,cv,v,+,+)
  # Probability of observing double crossover in intervals y-cv and f-car
  Pr10001<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%summation_vector)*((v[1]*1.00*1.00*1.00*v[5])/2) #Pr(y,+,+,+,car)
  Pr01110<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%N_3%*%R_4%*%summation_vector)*((1.00*v[2]*v[3]*v[4]*1.00)/2) #Pr(+,cv,v,f,+)
  # Probability of observing double crossover in intervals cv-v and v-f
  Pr11011<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%summation_vector)*((v[1]*v[2]*1.00*v[4]*v[5])/2) #Pr(y,cv,+,f,car)
  Pr00100<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%N_4%*%summation_vector)*((1.00*1.00*v[3]*1.00*1.00)/2) #Pr(+,+,v,+,+)
  # Probability of observing double crossover in intervals cv-v and f-car
  Pr10011<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%summation_vector)*((v[1]*1.00*1.00*v[4]*v[5])/2) #Pr(y,cv,+,+,car)
  Pr01100<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%N_4%*%summation_vector)*((1.00*v[2]*v[3]*1.00*1.00)/2) #Pr(+,+,v,f,+)
  # Probability of observing double crossover in intervals v-f and f-car
  Pr10111<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%summation_vector)*((v[1]*1.00*v[3]*v[4]*v[5])/2) #Pr(y,cv,v,+,car)
  Pr01000<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%N_4%*%summation_vector)*((1.00*v[2]*1.00*1.00*1.00)/2) #Pr(+,+,+,f,+)
  
  # Probability of observing triple crossover in intervals y-cv cv-v v-f
  Pr00101<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%summation_vector)*((1.00*1.00*v[3]*1.00*v[5])/2) #Pr(y,+,v,+,+)
  Pr11010<-((1/p[1])%*%stationary_vector%*%N_1%*%R_2%*%R_3%*%R_4%*%summation_vector)*((v[1]*v[2]*1.00*v[4]*1.00)/2) #Pr(+,cv,+,f,car)
  # Probability of observing triple crossover in intervals y-cv cv-v f-car
  Pr01101<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%summation_vector)*((1.00*v[2]*v[3]*1.00*v[5])/2) #Pr(y,+,v,f,+)
  Pr10010<-((1/p[1])%*%stationary_vector%*%R_1%*%N_2%*%R_3%*%R_4%*%summation_vector)*((v[1]*1.00*1.00*v[4]*1.00)/2) #Pr(+,cv,+,+,car)
  # Probability of observing triple crossover in intervals cv-v v-f f-car
  Pr01001<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%summation_vector)*((1.00*v[2]*1.00*1.00*v[5])/2) #Pr(y,+,+,f,+)
  Pr10110<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%N_3%*%R_4%*%summation_vector)*((v[1]*1.00*v[3]*v[4]*1.00)/2) #Pr(+,cv,v,+,car)
  # Probability of observing triple crossover in intervals cv-v v-f f-car
  Pr01011<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%summation_vector)*((1.00*v[2]*1.00*v[4]*v[5])/2) #Pr(y,cv,+,f,+)
  Pr10100<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%N_4%*%summation_vector)*((v[1]*1.00*v[3]*1.00*1.00)/2) #Pr(+,+,v,+,car)
  
  # Probability of observing quadruple crossover in intervals y-cv cv-v v-f f-car
  Pr10101<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%summation_vector)*((v[1]*1.00*v[3]*1.00*v[5])/2) #Pr(y,+,v,+,car)
  Pr01010<-((1/p[1])%*%stationary_vector%*%R_1%*%R_2%*%R_3%*%R_4%*%summation_vector)*((1.00*v[2]*1.00*v[4]*1.00)/2) #Pr(+,cv,+,f,+)
  
  # Probability of zygotic lethality due only to differential marker viability effects
  Pr_dead<-1-sum(Pr11111,
                 Pr00000,
                 Pr00001,
                 Pr11110,
                 Pr00011,
                 Pr11100,
                 Pr00111,
                 Pr11000,
                 Pr01111,
                 Pr10000,
                 Pr11101,
                 Pr00010,
                 Pr11001,
                 Pr00110,
                 Pr10001,
                 Pr01110,
                 Pr11011,
                 Pr00100,
                 Pr10011,
                 Pr01100,
                 Pr10111,
                 Pr01000,
                 Pr00101,
                 Pr11010,
                 Pr01101,
                 Pr10010,
                 Pr01001,
                 Pr10110,
                 Pr01011,
                 Pr10100,
                 Pr10101,
                 Pr01010)
  
  # Vector of expected frequencies of each recombination pattern
  expected_freqs <- c(Pr11111,
                      Pr00000,
                      Pr00001,
                      Pr11110,
                      Pr00011,
                      Pr11100,
                      Pr00111,
                      Pr11000,
                      Pr01111,
                      Pr10000,
                      Pr11101,
                      Pr00010,
                      Pr11001,
                      Pr00110,
                      Pr10001,
                      Pr01110,
                      Pr11011,
                      Pr00100,
                      Pr10011,
                      Pr01100,
                      Pr10111,
                      Pr01000,
                      Pr00101,
                      Pr11010,
                      Pr01101,
                      Pr10010,
                      Pr01001,
                      Pr10110,
                      Pr01011,
                      Pr10100,
                      Pr10101,
                      Pr01010,
                      Pr_dead)
  
  # raw_log_likelihoods for each of the individual recombination classes
  # log_likelihood is the sum of the individual recombination classes
  
  # Function returns negative log-likelihood value because optimization
  # with the Nelder-Mead method in optim() is, by default, minimization 
  raw_log_likelihoods <- observed_count*log(expected_freqs)
  log_likelihood <- sum(raw_log_likelihoods)
  return(-log_likelihood)
}

# Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
# Log likelihoods calculated for each combination of y & p model parameters
MLE_Uniform_CxCoM <- matrix(nrow=0, ncol=((4*intervals)+4))
colnames(MLE_Uniform_CxCoM)<-c("MLE_y1","MLE_y2","MLE_y3","MLE_y4","p1","p2","p3","p4","x1","x2","x3","x4","v1","v2","v3","v4","v5","lnL","Iterations","Convergence")

# For loop Uniform.CxCoM.Likelihood for all possible combinations of p
# Output model parameters and likelihoods in matrix MLE_Uniform_CxCoM
for (u in 1:max_p_value) {
  
  # Define parameters y and p for single run of optimization algorithm
  parameters <- c(y_uniform[u,],v_uniform[u,])
  p <- p_uniform[u,]
  
  # Optimize CxCoM model parameter "y" with Nelder-Mead method
  Uniform_CxCoM_Table <- optim(par = parameters,
                               fn = Uniform.CxCoM.Likelihood,
                               gr = NULL,
                               method = "Nelder-Mead", 
                               control = list(maxit = 10000,
                                              alpha = 2.0,
                                              beta = 0.25,
                                              gamma = 4.0))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_outputs <- c(Uniform_CxCoM_Table$par[1],
                   Uniform_CxCoM_Table$par[2],
                   Uniform_CxCoM_Table$par[3],
                   Uniform_CxCoM_Table$par[4],
                   p[1],
                   p[2],
                   p[3],
                   p[4],
                   Uniform_CxCoM_Table$par[1]/(2*p[1]),
                   Uniform_CxCoM_Table$par[2]/(2*p[2]),
                   Uniform_CxCoM_Table$par[3]/(2*p[3]),
                   Uniform_CxCoM_Table$par[4]/(2*p[4]),
                   Uniform_CxCoM_Table$par[5],
                   Uniform_CxCoM_Table$par[6],
                   Uniform_CxCoM_Table$par[7],
                   Uniform_CxCoM_Table$par[8],
                   Uniform_CxCoM_Table$par[9],
                   Uniform_CxCoM_Table$value,
                   Uniform_CxCoM_Table$counts[1],
                   Uniform_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_Uniform_CxCoM <- rbind(MLE_Uniform_CxCoM, MLE_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
output<-as.data.frame(MLE_Uniform_CxCoM)
df2 <- output[order(output$lnL),]
df2[1,]

#write.csv(df2,"~/Desktop/test_observed_output.csv",row.names = TRUE)


