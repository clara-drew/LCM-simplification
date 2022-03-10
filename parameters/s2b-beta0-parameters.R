# Parameters for scenario 2b with beta=0

maxT <- 3:8 # number of test results
n <- 100 # sample size
its <- 1000 # number of iterations
beta <- 0 # 1-specificity
theta <- 0.3 # prevalence of shedders

# Define alpha for each number of tests (alpha_list)
sigma <- 1
a0 <- 2
a2 <- -7/(maxT-1)
alpha_list <- list()
omega0 <- 1
omega1 <- 0.5
for(j in 1:length(maxT)){
  vt <- a0+a2[j]*(0:(maxT[j]-1))
  alpha_list[[j]] <- invlogit(omega0+omega1*vt)
}

# Function to simulate sensitivity for each number of tests
alpha_fun <- alpha_logistic_unident

folder <- "s2b-beta0" # folder name to save results
system(paste("mkdir -p ", folder, sep="")) # Create folder if it does not already exist
ncores <- 23 # number of cores to use for parallel bootstrapping (set to 1 if not computing in parallel)
