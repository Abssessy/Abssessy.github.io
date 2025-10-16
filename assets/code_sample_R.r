# ========================
# INTERNATIONAL SANCTIONS AND INTERNATIONAL STUDENT FLOWS
# Code Sample
# ========================

# ------------------------
# 1. LOAD LIBRARIES
# ------------------------
library(igraph)
library(network)
library(ergm)
library(btergm)
library(xergm)
library(sna)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(doParallel)
library(matrixStats)
library(data.table)
library(parallel)
library(VCERGM)     
library(splines)    
library(gridExtra) 
library(statnet.common)
library(Rcpp)      
library(RcppArmadillo)
library(MASS)       
library(sdpt3r)   
library(tidyverse)
library(doParallel)

# ------------------------
# 2. DATA IMPORT FUNCTION (with error handling)
# ------------------------
load_network_data <- function(years, file_prefix) {
  data_list <- list()
  for (year in years) {
    filename <- paste0(file_prefix, year, ".csv")
    if (!file.exists(filename)) {
      warning(paste("File not found:", filename))
      next
    }
    data_list[[as.character(year)]] <- as.matrix(read.csv(filename, header = TRUE, row.names = 1))
  }
  return(data_list)
}

years <- 2000:2020
student_networks <- load_network_data(years, "stud")
sanction_networks <- load_network_data(years, "sanc")
financial_networks <- load_network_data(years, "finan")
trade_networks <- load_network_data(years, "trade")
travel_networks <- load_network_data(years, "travel")

# Load time-invariant networks
language <- as.matrix(read.csv("lan.csv", header = TRUE, row.names = 1))
distance <- as.matrix(read.csv("dist.csv", header = TRUE, row.names = 1)) * 0.0001
colonial <- as.matrix(read.csv("colo.csv", header = TRUE, row.names = 1))

# Load node attributes (assuming same years)
node_attributes <- load_network_data(years, "")

# ------------------------
# 3. NETWORK CONSTRUCTION
# ------------------------
create_binary_network <- function(matrix) {
  threshold <- quantile(matrix, probs = 0.6, na.rm = TRUE)
  binary_net <- ifelse(matrix > threshold, 1, 0)
  diag(binary_net) <- 0
  return(binary_net)
}

binary_student_networks <- lapply(student_networks, create_binary_network)

country_names <- c("AGO","ALB","ARE","ARG","ARM","ATG","AUS","AUT","AZE",
                   "BDI","BEL","BEN","BFA","BGD","BGR","BHR","BHS","BIH",
                   "BLR","BLZ","BMU","BOL","BRA","BRB","BRN","BTN","BWA",
                   "CAF","CAN","CHE","CHL","CHN","CIV","CMR","COG","COL",
                   "COM","CPV","CRI","CYP","CZE","DEU","DJI","DMA","DNK",
                   "DOM","DZA","ECU","EGY","ERI","ESP","EST","ETH","FIN",
                   "FJI","FRA","FSM","GAB","GBR","GEO","GHA","GIN","GMB",
                   "GNB","GNQ","GRC","GRD","GTM","GUY","HKG","HND","HRV",
                   "HTI","HUN","IDN","IND","IRL","IRN","IRQ","ISL","ISR",
                   "ITA","JAM","JOR","JPN","KAZ","KEN","KGZ","KHM","KOR",
                   "KWT","LAO","LBN","LBR","LBY","LCA","LKA","LSO","LTU",
                   "LUX","LVA","MAC","MAR","MDA","MDG","MDV","MEX","MKD",
                   "MLI","MLT","MMR","MNG","MOZ","MRT","MUS","MWI","MYS",
                   "NAM","NER","NGA","NIC","NLD","NOR","NPL","NZL","OMN",
                   "PAK","PAN","PER","PHL","PNG","POL","PRI","PRT","PRY",
                   "QAT","RUS","RWA","SAU","SDN","SEN","SGP","SLB","SLE",
                   "SLV","SUR","SVK","SVN","SWE","SWZ","SYC","TCD","TGO",
                   "THA","TJK","TKM","TTO","TUN","TUR","TZA","UGA","UKR",
                   "URY","USA","UZB","VCT","VEN","VNM","VUT","WSM","YEM",
                   "ZAF","ZMB","ZWE")

for (year in years) {
  rownames(binary_student_networks[[as.character(year)]]) <- country_names
  colnames(binary_student_networks[[as.character(year)]]) <- country_names
}

# ------------------------
# 4. NETWORK VISUALIZATION (Parameterized years)
# ------------------------
plot_network <- function(net_matrix, year) {
  net <- network(net_matrix, directed = TRUE)
  graph_tbl <- as_tbl_graph(net) %>% 
    mutate(deg = centrality_degree(mode = 'in'),
           group = group_infomap())
  
  ggraph(graph_tbl, layout = 'kk') + 
    geom_edge_fan(color = "lightblue", show.legend = FALSE) + 
    geom_node_point(aes(fill = factor(group)), shape = 21) + 
    geom_node_text(aes(filter = deg > 12, label = name), size = 2) +
    scale_color_discrete() +
    scale_edge_width(range = c(0.2, 10)) +
    guides(size = "none", fill = "none") +
    theme_graph() +
    ggtitle(paste("International Student Network", year))
}

viz_years <- c(2000, 2005, 2010, 2015, 2020)
for (viz_year in viz_years) {
  print(plot_network(binary_student_networks[[as.character(viz_year)]], viz_year))
}

# ------------------------
# 5. NETWORK DESCRIPTIVES
# ------------------------
calculate_network_stats <- function(net_matrix) {
  net_graph <- graph.adjacency(net_matrix, mode = "directed")
  return(list(
    mean_degree = mean(rowSums(net_matrix)),
    sd_degree = sd(rowSums(net_matrix)),
    eigen_centrality = mean(eigen_centrality(net_graph)$vector),
    sd_eigen_centrality = sd(eigen_centrality(net_graph)$vector),
    transitivity = transitivity(net_graph, type = "global"),
    diameter = diameter(net_graph)
  ))
}

network_stats <- lapply(binary_student_networks[as.character(viz_years)], calculate_network_stats)

# TABLE 4: Dynamic Trend Similarity Test
trend <- read.csv("trend_analysis.csv", header = TRUE, row.names = 1)

# Function to calculate all similarity metrics for a pair of columns
calc_similarity <- function(col1, col2, col_names) {
  list(
    pearson = cor.test(col1, col2, method = "pearson"),
    mape = MLmetrics::MAPE(col1, col2),
    dtw = dtw::dtw(col1, col2, step = asymmetricP1, keep = TRUE)
  )
}

# Calculate similarities for all three metrics
metrics <- c("Link_number", "Transitivity", "Reciprocity")
similarity_results <- map(1:3, ~ calc_similarity(trend[,.x], trend[,.x+3], metrics[.x]))

# Optional: Plot DTW alignments
walk(similarity_results, ~ dtwPlotTwoWay(.x$dtw))

# QAP Regression Analysis
library(qap)
library(sna)

# Function to create QAP covariates for a given year
create_qap_covariates <- function(year_data, year_suffix) {
  list(
    Phom = as.matrix(dist(year_data[,1])),
    Pshom = as.matrix(dist(year_data[,3])),
    Prec = matrix(year_data[,1], nrow(year_data), nrow(year_data), byrow = TRUE),
    Psrec = matrix(year_data[,3], nrow(year_data), nrow(year_data), byrow = TRUE),
    eqrec = matrix(year_data[,8], nrow(year_data), nrow(year_data), byrow = TRUE),
    s = get(paste0("s", year_suffix)),  # Assuming these are pre-loaded
    dist = dist,
    colo = colo,
    lan = lan
  )
}

# Run QAP analysis for all years
years <- c(2000, 2003, 2005, 2016, 2020)
year_data_list <- list(a2000, a2003, a2005, a2016, a2020)
network_list <- list(bi2000, bi2003, bi2005, bi2016, bi2020)

qap_results <- map2(year_data_list, years, ~ {
  covariates <- create_qap_covariates(.x, .y)
  netlm(network_list[[which(years == .y)]], covariates, nullhyp = "qap", reps = 100)
})

names(qap_results) <- paste0("qap", years)

# ------------------------
# 6. ECONOMETRIC ANALYSIS
# ------------------------
# 6.1 ERGM Analysis
prepare_ergm_data <- function(year) {
  net_year <- binary_student_networks[[as.character(year)]]
  attr_year <- node_attributes[[as.character(year)]]
  sanction_year <- sanction_networks[[as.character(year)]]
  
  network_obj <- network(net_year, 
                         vertex.attr = attr_year,
                         vertex.attrnames = c('pgdp','va','ps','ge','rq','rl','coc','eq','sp','ss','sq'),
                         directed = TRUE)
  
  return(list(network = network_obj, 
              sanctions = sanction_year,
              distance = distance,
              language = language,
              colonial = colonial))
}

ergm_results <- list()
for (ergm_year in c(2000, 2005, 2010, 2015, 2020)) {
  data <- prepare_ergm_data(ergm_year)
  
  set.seed(234)
  ergm_fit <- ergm(data$network ~ 
                     edgecov(data$sanctions) +
                     edges + 
                     dgwesp(2, fixed = TRUE, type = "RTP") + 
                     mutual +
                     absdiff("pgdp") + 
                     absdiff("ps") +
                     nodeocov("pgdp") + 
                     nodeocov("ps") + 
                     nodeocov("eq") +
                     nodeicov("pgdp") + 
                     nodeicov("ps") + 
                     nodeicov("eq") +
                     edgecov(data$distance) + 
                     edgecov(data$colonial) + 
                     edgecov(data$language),
                   control = control.ergm(MCMC.samplesize = 5000, 
                                          MCMC.burnin = 10000, 
                                          MCMLE.maxit = 10),
                   verbose = TRUE)
  
  ergm_results[[as.character(ergm_year)]] <- summary(ergm_fit)
}

# 6.2 TERGM Analysis
# Prepare panel data
panel_networks <- binary_student_networks
panel_sanctions <- sanction_networks

set.seed(123)
tergm_fit <- btergm(panel_networks ~ 
                      edges + 
                      dgwesp(2, fixed = TRUE, type = "OTP") + 
                      mutual +
                      absdiff("pgdp") + 
                      absdiff("ps") +
                      nodeocov("pgdp") + 
                      nodeocov("ps") + 
                      nodeocov("eq") +
                      nodeicov("pgdp") + 
                      nodeicov("ps") + 
                      nodeicov("eq") +
                      edgecov(panel_sanctions) + 
                      edgecov(distance) + 
                      edgecov(colonial) + 
                      edgecov(language) + 
                      timecov(panel_sanctions) +
                      memory(type = "stability"),
                    R = 100)

# Extract coefficients for counterfactual analysis
tergm_coef <- coef(tergm_fit)

# ------------------------
# 6.3. SIMPLE REGRESSION ANALYSIS
# ------------------------
# Compare two time periods (2000 vs 2003 as example)
first <- bi2000
third <- bi2003

N <- dim(first)[1] # Number of agents 

# Create regression design matrix
X <- cbind(
  as.vector(matrix(1, 2*N, 1)),
  c(as.vector(matrix(0, N, 1)), as.vector(matrix(1, N, 1)))
)

# Create response variables for three network measures
graphFirst <- graph.adjacency(first, mode = "directed")
graphThird <- graph.adjacency(third, mode = "directed")

Y1 <- c(rowSums(first), rowSums(third))  # Degree
Y2 <- c(
  as.vector(eigen_centrality(graphFirst)$vector),
  as.vector(eigen_centrality(graphThird)$vector)
)  # Eigenvector centrality

Y3 <- c(
  transitivity(graphFirst, type = "local"),
  transitivity(graphThird, type = "local")
)  # Local clustering

# Regression for degree
beta1 <- ginv(t(X) %*% X) %*% (t(X) %*% Y1)
variance1 <- as.double(t(Y1 - X %*% beta1) %*% (Y1 - X %*% beta1)) * ginv(t(X) %*% X) / (2*N - 2)

# Regression for eigenvector centrality
beta2 <- ginv(t(X) %*% X) %*% (t(X) %*% Y2)
variance2 <- as.double(t(Y2 - X %*% beta2) %*% (Y2 - X %*% beta2)) * ginv(t(X) %*% X) / (2*N - 2)

# Regression for clustering (handle NaN values)
X3 <- X[is.nan(Y3) == 0, ]
Y3 <- Y3[is.nan(Y3) == 0]
beta3 <- ginv(t(X3) %*% X3) %*% (t(X3) %*% Y3)
variance3 <- as.double(t(Y3 - X3 %*% beta3) %*% (Y3 - X3 %*% beta3)) * ginv(t(X3) %*% X3) / (2*N - 2)

# Calculate t-statistics and p-values
tstat11 <- beta1[1]^2 / (variance1[1, 1])
pval11 <- 1 - pchisq(tstat11, 1)

tstat12 <- beta1[2]^2 / (variance1[2, 2])
pval12 <- 1 - pchisq(tstat12, 1)

tstat21 <- beta2[1]^2 / (variance2[1, 1])
pval21 <- 1 - pchisq(tstat21, 1)

tstat22 <- beta2[2]^2 / (variance2[2, 2])
pval22 <- 1 - pchisq(tstat22, 1)

tstat31 <- beta3[1]^2 / (variance3[1, 1])
pval31 <- 1 - pchisq(tstat31, 1)

tstat32 <- beta3[2]^2 / (variance3[2, 2])
pval32 <- 1 - pchisq(tstat32, 1)

# Create Table 2
table2 <- data.frame(
  Variable = rep(c("Intercept", "Time"), 3),
  Measure = rep(c("Degree", "Eigenvector", "Clustering"), each = 2),
  Coefficient = c(beta1[1], beta1[2], beta2[1], beta2[2], beta3[1], beta3[2]),
  P_Value = c(pval11, pval12, pval21, pval22, pval31, pval32)
)

print("Table 2: Regression Results")
print(table2)

# ------------------------
# 6.4. RANDOMIZATION TEST USING THE 2-2 NORM (Auerbach, 2022)
# ------------------------
# Setup parallel processing
registerDoParallel(cores = detectCores() - 1)

# Use 2000 and 2003 as example years
first <- bi2000
third <- bi2003
N <- dim(first)[1] # Number of agents

# Summary statistics for the actual data
graphFirst <- graph.adjacency(first)
graphThird <- graph.adjacency(third)

# Degree statistics
totDeg <- abs(mean(first) - mean(third))
diffDeg <- mean((rowSums(first) - rowSums(third))^2)

# Eigenvector centrality statistics
diffCent <- mean((
  as.vector(eigen_centrality(graph.adjacency(first))$vector) -
    as.vector(eigen_centrality(graph.adjacency(third))$vector)
)^2)

# Clustering statistics
diffTran <- abs(
  transitivity(graphThird, type = "global") - 
    transitivity(graphFirst, type = "global")
)

# Diameter statistics
diffDiam <- abs(diameter(graphThird) - diameter(graphFirst))

# All statistics together
stats <- c(totDeg, diffDeg, diffCent, diffTran, diffDiam)

# Randomization tests using summary statistics
permuteStats <- foreach(r = 1:100, .combine = 'rbind', 
                        .packages = c("igraph", "foreach", "MASS", "Matrix", "sdpt3r")) %dopar% {
  # Create random sign matrix
  v <- matrix(sign(runif(N^2, -1, 1)), N)
  v <- upper.tri(v, diag = FALSE) * v + t(upper.tri(v, diag = FALSE) * v) # Make symmetric
  
  # Permute networks
  permuteFirst <- pmax(first * (v == 1), third * (v == -1))
  permuteThird <- pmax(third * (v == 1), first * (v == -1))
  
  # Calculate statistics for permuted data
  permuteTotDeg <- abs(mean(permuteFirst) - mean(permuteThird))
  permuteDiffDeg <- mean((rowSums(permuteFirst) - rowSums(permuteThird))^2)
  
  permuteDiffCent <- mean((
    as.vector(eigen_centrality(graph.adjacency(permuteFirst))$vector) -
      as.vector(eigen_centrality(graph.adjacency(permuteThird))$vector)
  )^2)
  
  pemuteDiffTran <- abs(
    transitivity(graph.adjacency(permuteThird), type = "global") - 
      transitivity(graph.adjacency(permuteFirst), type = "global")
  )
  
  permuteDiffDiam <- abs(
    diameter(graph.adjacency(permuteThird)) - 
      diameter(graph.adjacency(permuteFirst))
  )
  
  c(permuteTotDeg, permuteDiffDeg, permuteDiffCent, pemuteDiffTran, permuteDiffDiam)
}

# Add actual statistics to permutation results
permuteStats <- rbind(permuteStats, stats)

# Calculate quantiles and p-values
quantiles <- t(colQuantiles(permuteStats, probs = c(.5, .75, .9, .95, .97, .99)))
p_values <- c(
  mean(permuteStats[, 1] >= stats[1]),
  mean(permuteStats[, 2] >= stats[2]),
  mean(permuteStats[, 3] >= stats[3]),
  mean(permuteStats[, 4] >= stats[4]),
  mean(permuteStats[, 5] >= stats[5])
)

# Create Table 3 part 1
table3_part1 <- data.frame(
  Statistic = c("Total Degree", "Degree Difference", "Eigenvector Difference", 
                "Clustering Difference", "Diameter Difference"),
  Value = stats,
  P_Value = p_values
)

print("Table 3 Part 1: Randomization Test Results")
print(table3_part1)
print("Quantiles of Reference Distribution:")
print(quantiles)

# 2-2 norm (spectral norm) test
spectralNorm <- max(svd(first - third)$d)

# Randomization test for spectral norm
permuteSpectral <- foreach(r = 1:100, .combine = 'rbind', 
                           .packages = c("igraph", "foreach", "MASS", "Matrix", "sdpt3r")) %dopar% {
  # Create random sign matrix
  v <- matrix(sign(runif(N^2, -1, 1)), N)
  v <- upper.tri(v, diag = FALSE) * v + t(upper.tri(v, diag = FALSE) * v) # Make symmetric
  
  # Permute networks
  permuteFirst <- pmax(first * (v == 1), third * (v == -1))
  permuteThird <- pmax(third * (v == 1), first * (v == -1))
  
  # Calculate spectral norm for permuted data
  max(svd(permuteFirst - permuteThird)$d)
}

# Add actual spectral norm to permutation results
permuteSpectral <- rbind(permuteSpectral, spectralNorm)

# Calculate results
spectral_quantiles <- cbind(colQuantiles(permuteSpectral, probs = c(.5, .75, .9, .95, .97, .99)))
spectral_pvalue <- mean(spectralNorm <= permuteSpectral)

# Create Table 3 part 2
table3_part2 <- data.frame(
  Statistic = "Spectral Norm (2-2 Norm)",
  Value = spectralNorm,
  P_Value = spectral_pvalue
)

print("Table 3 Part 2: Spectral Norm Test")
print(table3_part2)
print("Quantiles of Reference Distribution:")
print(spectral_quantiles)

# Stop parallel cluster
stopImplicitCluster()

# ------------------------
# 7. HETEROGENEITY ANALYSIS
# ------------------------
# 7.1 Different sanction types
tergm_trade <- btergm(panel_networks ~ 
                        edges + mutual + 
                        edgecov(trade_networks) + 
                        edgecov(panel_sanctions) +
                        edgecov(distance) + 
                        edgecov(colonial) + 
                        edgecov(language),
                      R = 100)

tergm_financial <- btergm(panel_networks ~ 
                            edges + mutual + 
                            edgecov(financial_networks) + 
                            edgecov(panel_sanctions) +
                            edgecov(distance) + 
                            edgecov(colonial) + 
                            edgecov(language),
                          R = 100)

# 7.2 Interaction effects
create_interaction_effects <- function(sanction_mat, moderator_mat) {
  interaction <- sanction_mat * moderator_mat
  return(scale(interaction, center = TRUE, scale = TRUE))
}

# GDP distance interactions
gdp_interactions <- list()
for (year in years) {
  gdp_interactions[[as.character(year)]] <- 
    create_interaction_effects(sanction_networks[[as.character(year)]], 
                               distance)
}

# ------------------------
# 8. COUNTERFACTUAL ANALYSIS
# ------------------------
# Use estimated coefficients from TERGM for counterfactual analysis
par <- list(
  beta = list(
    sanction   = tergm_coef["edgecov.panel_sanctions[[i]]"],
    language   = tergm_coef["edgecov.language"],
    colony     = tergm_coef["edgecov.colonial"],
    distance   = tergm_coef["edgecov.distance"],
    constant   = tergm_coef["edges"],
    gwesp      = tergm_coef["dgwesp.OTP.2"],
    mutual     = tergm_coef["mutual"],
    diff_gdp   = tergm_coef["absdiff.pgdp"],
    diff_trade = 0,  # Placeholder - not in original model
    gdp_sender = c(tergm_coef["nodeocov.pgdp"], tergm_coef["nodeocov.ps"], tergm_coef["nodeocov.eq"]),
    gdp_target = c(tergm_coef["nodeicov.pgdp"], tergm_coef["nodeicov.ps"], tergm_coef["nodeicov.eq"])
  ),
  gw_decay = 0.25,   # geometrically-weighted decay
  thresh   = 0.20,   # binary-adjacency threshold (quantile)
  max_esp  = 3L,     # max shared partners
  n_nodes  = length(country_names),
  iter_max = 15L     # max counter-factual iterations
)

# Helper: GWESP change statistic
compute_gwesp <- function(adj, i, j, decay, max_esp) {
  adj1 <- adj0 <- adj
  adj1[i, j] <- 1; adj0[i, j] <- 0
  ch <- summary(adj1 ~ esp(0:max_esp)) - summary(adj0 ~ esp(0:max_esp))
  ch <- as.numeric(ch)[-1]               # drop esp=0 term
  sum(exp(decay) * (1 - (1 - exp(-decay))^(seq_along(ch))) * ch)
}

# Single-edge probability (logistic)
link_prob <- function(i, j, adj, par, X) {
  mutual  <- adj[j, i]
  gwesp   <- compute_gwesp(adj, i, j, par$gw_decay, par$max_esp)
  san     <- X$san[i, j]
  lang    <- X$lang[i, j]
  colony  <- X$col[i, j]
  dist    <- X$dist[i, j]

  sender_gdp <- X$gdp[i, ]; target_gdp <- X$gdp[j, ]
  diff_gdp   <- abs(sender_gdp[1] - target_gdp[1])
  diff_trade <- abs(sender_gdp[3] - target_gdp[3])

  idx <- 1:3
  linPred <- c(
    san    * par$beta$sanction,
    lang   * par$beta$language,
    colony * par$beta$colony,
    dist   * par$beta$distance,
    par$beta$constant,
    gwesp  * par$beta$gwesp,
    mutual * par$beta$mutual,
    diff_gdp   * par$beta$diff_gdp,
    diff_trade * par$beta$diff_trade,
    sum(sender_gdp * par$beta$gdp_sender[idx]),
    sum(target_gdp * par$beta$gdp_target[idx])
  )
  plogis(sum(linPred))
}

# Full network probability matrix
prob_matrix <- function(adj, par, X) {
  n <- nrow(adj)
  out <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {  # Skip self-loops
        out[i, j] <- link_prob(i, j, adj, par, X)
      }
    }
  }
  out
}

# Counterfactual iteration
iterate_cf <- function(adj_init, prob_init, treat_i, treat_j, par, X) {
  adj  <- adj_init
  prob <- prob_init
  converged <- FALSE
  k <- 0
  while (!converged && k < par$iter_max) {
    k <- k + 1
    # flip treated edge
    adj[treat_i, treat_j] <- 1
    prob[treat_i, treat_j] <- link_prob(treat_i, treat_j, adj, par, X)

    # recompute entire network
    prob_new <- prob_matrix(adj, par, X)
    adj_new  <- 1 * (prob_new > quantile(prob_new, par$thresh, na.rm = TRUE))

    converged <- identical(adj_new, adj)
    adj <- adj_new; prob <- prob_new
  }
  list(adj = adj, prob = prob, iter = k, converged = converged)
}

# Prepare data for counterfactual analysis
X <- list(
  san  = as.matrix(sanction_networks[["2000"]]),
  lang = as.matrix(language),
  col  = as.matrix(colonial),
  dist = as.matrix(distance),
  gdp  = node_attributes[["2000"]]
)

# Baseline network for 2000
adj0  <- binary_student_networks[["2000"]]
prob0 <- prob_matrix(adj0, par, X)

# Example: Remove sanction on CHN→USA (indices based on country_names)
treat_i <- which(country_names == "CHN")
treat_j <- which(country_names == "USA")

# Run counterfactual analysis
res <- iterate_cf(adj0, prob0, treat_i, treat_j, par, X)

# Summary of results
cat(sprintf("Converged in %d iterations\n", res$iter))
short_direct <- res$prob[treat_i, treat_j] - prob0[treat_i, treat_j]
cat(sprintf("Short-run direct effect: %.4f\n", short_direct))

delta_prob <- res$prob - prob0
delta_prob[treat_i, treat_j] <- NA
spillover <- mean(delta_prob, na.rm = TRUE)
cat(sprintf("Average spill-over effect: %.4f\n", spillover))

# ------------------------
# 9. ADVANCED NETWORK ANALYSIS
# ------------------------
# 9.1 VCERGM (Varying-Coefficient ERGM) Analysis (Lee et al., 2020)
# Prepare network list for VCERGM
xnet <- binary_student_networks

# Model setup
object_formula <- Net ~ edges + mutual + dgwesp(2, fixed = TRUE, type = "RTP")
stat_names <- c("edges", "mutual", "gwesp")
directed <- TRUE

# Degree of spline and number of knots for basis expansion
degree_spline <- 3
interior_knot <- 10

# Cross-sectional ERGM estimation
ergmest <- cross_sectional_ergm(
  object = object_formula, 
  network = xnet, 
  directed = directed, 
  degree.spline = degree_spline, 
  interior.knot = interior_knot
)

crossERGM_est <- ergmest$phi.hat # Cross-sectional ERGMs: hat(phi(t))
adhoc_est <- ergmest$phi.hat.smooth # Ad hoc 2-step procedure: hat(phi(t))

# VCERGM estimation
vcergmest <- estimate_vcergm(
  object = object_formula, 
  network = xnet, 
  degree.spline = degree_spline, 
  interior.knot = interior_knot, 
  directed = directed, 
  constant = FALSE
)

vcergm_est <- vcergmest$phi.hat # VCERGM: hat(phi(t))

# Format results for comparison
vcergm_est <- as.matrix(vcergm_est)
crossERGM_est <- as.matrix(crossERGM_est)
rownames(vcergm_est) <- c("vc_edges", "vc_mutual", "vc_gwesp")
colnames(vcergm_est) <- paste0("N", 1:ncol(vcergm_est))

rownames(crossERGM_est) <- c("cross_edges", "cross_mutual", "cross_gwesp")
colnames(crossERGM_est) <- paste0("N", 1:ncol(crossERGM_est))

# Combine results for comparison
combined_results <- rbind(vcergm_est, crossERGM_est)

# 9.2 PSTERGM (Panel Stochastic TERGM) Analysis (Kei et al., 2023)

# Find maximum value across all networks for PSTERGM
max_val <- max(sapply(student_networks, max, na.rm = TRUE))

# Create network list
y_list <- student_networks

# For demonstration, we'll create a placeholder for node attributes
node_attr <- matrix(0, nrow = length(country_names), ncol = ncol(node_attributes[["2000"]]))
colnames(node_attr) <- colnames(node_attributes[["2000"]])

# Placeholder for PSTERGM parameter estimation
# In practice, this would require compiling and running C++ code
cat("PSTERGM analysis would require compiling C++ code from:\n")
cat("https://github.com/allenkei/PSTERGM\n")
cat("This is a placeholder for the PSTERGM analysis.\n")

# The following would be the R code to call the compiled C++ functions:
eta <- rep(0, 8)
eta <- partial_stepping_temporal(20, 100, n, y_list, eta, node_attr, pi0 = 0.2, m = max_val + 1)
eta <- newton_raphson_temporal(10, 100, 2*n, y_list, eta, node_attr, pi0 = 0.2, m = max_val + 1)
write.csv(eta, 'eta.csv', row.names = FALSE)
se_eta <- SE_temporal(100, 2*n, y_list, eta, node_attr, pi0 = 0.2, m = max_val + 1)
write.csv(se_eta, "se_eta.csv", row.names = FALSE)

# ------------------------
# 10. RESULTS VISUALIZATION AND COMPARISON
# ------------------------
# Compare VCERGM and cross-sectional ERGM results
comparison_df <- data.frame(
  Method = rep(c("VCERGM", "Cross-sectional"), each = ncol(combined_results)/2),
  Time = rep(2000:2020, 2),
  Edges = c(combined_results["vc_edges", ], combined_results["cross_edges", ]),
  Mutual = c(combined_results["vc_mutual", ], combined_results["cross_mutual", ]),
  GWESP = c(combined_results["vc_gwesp", ], combined_results["cross_gwesp", ])
)

# Plot comparison of edges parameter over time
edges_plot <- ggplot(comparison_df, aes(x = Time, y = Edges, color = Method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  ggtitle("Comparison of Edges Parameter Estimates") +
  ylab("Parameter Estimate")

# Plot comparison of mutual parameter over time
mutual_plot <- ggplot(comparison_df, aes(x = Time, y = Mutual, color = Method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  ggtitle("Comparison of Mutual Parameter Estimates") +
  ylab("Parameter Estimate")

# Plot comparison of GWESP parameter over time
gwesp_plot <- ggplot(comparison_df, aes(x = Time, y = GWESP, color = Method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  ggtitle("Comparison of GWESP Parameter Estimates") +
  ylab("Parameter Estimate")

# Display plots
print(edges_plot)
print(mutual_plot)
print(gwesp_plot)

# ========================
# END OF CODE SAMPLE
# ========================