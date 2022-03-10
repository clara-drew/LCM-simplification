# Parameters for scenario 3 with beta=0.05

maxT <- 3:8 # number of test results
n <- 100 # sample size
its <- 1000 # number of iterations
beta <- 0.05 # 1-specificity
theta <- 0.3 # prevalence of shedders

# Define alpha for each number of tests (alpha_list)
sigma <- NULL
a0 <- NULL
a2 <- NULL
alpha_list <- list()
for(j in 1:length(maxT)){
  alpha_list[[j]] <- alpha_polynomial(maxT)[1,]
}

# Function to simulate sensitivity for each number of tests
alpha_fun <- alpha_polynomial

folder <- "s3-beta05" # folder name to save results
system(paste("mkdir -p ", folder, sep="")) # Create folder if it does not already exist
ncores <- 23 # number of cores to use for parallel bootstrapping (set to 1 if not computing in parallel)
