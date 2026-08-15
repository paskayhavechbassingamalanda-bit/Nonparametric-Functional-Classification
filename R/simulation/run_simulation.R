###############################################################
# Generate a simulated sample
#
# Master's Thesis:
# Nonparametric Classification for Spatial Functional Data
# under Random Censoring
#
# Author: Paska Bassinga Malanda
###############################################################

# Load the data generation function
source("R/data_generation.R")

# ------------------------------------------------------------
# Simulation parameters
# ------------------------------------------------------------

n1 <- 10
n2 <- 10
a <- 0.1
sigma <- 0.5

# ------------------------------------------------------------
# Generate one sample
# ------------------------------------------------------------

Dn <- data2(
  n1 = n1,
  n2 = n2,
  a = a,
  sigma = sigma
)

# ------------------------------------------------------------
# Dimensions of the simulated data
# ------------------------------------------------------------

dim(Dn)

# Display the first observations
head(Dn)
