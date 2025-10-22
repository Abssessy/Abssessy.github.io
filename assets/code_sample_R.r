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
library(statnet)
library(openxlsx)

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

# Load and prepare network data
prepare_network_data <- function(years, country_names) {
  networks <- list()
  attributes <- list()
  covariates <- list()
  
  for (year in years) {
    # Load network data
    networks[[as.character(year)]] <- get_network_data(year)
    
    # Load node attributes
    attributes[[as.character(year)]] <- get_node_attributes(year)
    
    # Load covariate networks
    covariates[[as.character(year)]] <- list(
      sanction = get_sanction_network(year),
      language = get_language_network(),
      colonial = get_colonial_network(),
      distance = get_distance_network()
    )
  }
  
  return(list(
    networks = networks,
    attributes = attributes,
    covariates = covariates
  ))
}

# Estimate TERGM model
estimate_tergm_model <- function(networks, covariates, node_attributes) {
  # Prepare btergm formula
  formula <- prepare_tergm_formula(networks, covariates, node_attributes)
  
  # Estimate model
  set.seed(123)
  tergm_fit <- btergm(
    formula$formula,
    R = 100  # Number of bootstrap samples
  )
  
  return(tergm_fit)
}

# Setup counterfactual simulation parameters
setup_counterfactual_params <- function(tergm_coef, config = list()) {
  par <- list(
    beta = list(
      sanction   = tergm_coef["edgecov.sanction"],
      language   = tergm_coef["edgecov.language"],
      colony     = tergm_coef["edgecov.colonial"],
      distance   = tergm_coef["edgecov.distance"],
      constant   = tergm_coef["edges"],
      gwesp      = tergm_coef["gwesp"],
      mutual     = tergm_coef["mutual"],
      diff_gdp   = tergm_coef["absdiff.gdp"],
      diff_trade = tergm_coef["absdiff.trade"],
      gdp_sender = c(tergm_coef["nodeocov.gdp"], tergm_coef["nodeocov.size"]),
      gdp_target = c(tergm_coef["nodeicov.gdp"], tergm_coef["nodeicov.size"])
    ),
    gw_decay = config$gw_decay %||% 0.25,
    thresh   = config$thresh %||% 0.20,
    max_esp  = config$max_esp %||% 3L,
    n_nodes  = config$n_nodes,
    iter_max = config$iter_max %||% 15L
  )
  return(par)
}

# Calculate GWESP change statistic
compute_gwesp <- function(adj, i, j, decay, max_esp) {
  adj1 <- adj0 <- adj
  adj1[i, j] <- 1; adj0[i, j] <- 0
  
  # Calculate shared partner changes
  esp1 <- summary(adj1 ~ esp(0:max_esp))
  esp0 <- summary(adj0 ~ esp(0:max_esp))
  ch <- esp1 - esp0
  ch <- as.numeric(ch)[-1]  # Remove esp=0 term
  
  # Calculate geometric weights
  weights <- decay * (1 - (1 - exp(-decay))^(1:length(ch)))
  sum(weights * ch)
}

# Calculate single link probability
link_prob <- function(i, j, adj, par, X) {
  # Network structural features
  mutual  <- adj[j, i]
  gwesp   <- compute_gwesp(adj, i, j, par$gw_decay, par$max_esp)
  
  # Covariates
  san     <- X$sanction[i, j]
  lang    <- X$language[i, j]
  colony  <- X$colonial[i, j]
  dist    <- X$distance[i, j]
  
  # Node attribute differences
  sender_attr <- X$attributes[i, ]
  target_attr <- X$attributes[j, ]
  diff_gdp   <- abs(sender_attr["gdp"] - target_attr["gdp"])
  diff_trade <- abs(sender_attr["trade"] - target_attr["trade"])
  
  # Linear predictor
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
    sum(sender_attr * par$beta$gdp_sender),
    sum(target_attr * par$beta$gdp_target)
  )
  
  # Logistic transformation
  plogis(sum(linPred, na.rm = TRUE))
}

#' Calculate full network probability matrix
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

# Counterfactual iteration simulation
iterate_counterfactual <- function(adj_init, prob_init, treat_edges, par, X) {
  adj <- adj_init
  prob <- prob_init
  convergence_history <- list()
  
  for (k in 1:par$iter_max) {
    # Apply treatment effects
    for (edge in treat_edges) {
      i <- edge$i; j <- edge$j
      adj[i, j] <- edge$new_value
      prob[i, j] <- link_prob(i, j, adj, par, X)
    }
    
    # Recompute entire network
    prob_new <- prob_matrix(adj, par, X)
    threshold <- quantile(prob_new, par$thresh, na.rm = TRUE)
    adj_new <- 1 * (prob_new > threshold)
    
    # Check convergence
    converged <- identical(adj_new, adj)
    convergence_history[[k]] <- list(
      iteration = k,
      converged = converged,
      edges_changed = sum(adj_new != adj),
      mean_prob = mean(prob_new)
    )
    
    if (converged) break
    adj <- adj_new
    prob <- prob_new
  }
  
  return(list(
    adj = adj,
    prob = prob,
    iterations = k,
    converged = converged,
    history = convergence_history
  ))
}

# Run multi-period counterfactual analysis
run_multiperiod_counterfactual <- function(data, par, treat_spec, years) {
  results <- list()
  
  for (year in years) {
    cat("Processing year:", year, "\n")
    
    # Prepare current year data
    X <- list(
      sanction = data$covariates[[as.character(year)]]$sanction,
      language = data$covariates[[as.character(year)]]$language,
      colonial = data$covariates[[as.character(year)]]$colonial,
      distance = data$covariates[[as.character(year)]]$distance,
      attributes = data$attributes[[as.character(year)]]
    )
    
    # Initial network
    adj0 <- data$networks[[as.character(year)]]
    prob0 <- prob_matrix(adj0, par, X)
    
    # Run counterfactual simulation
    results[[as.character(year)]] <- iterate_counterfactual(
      adj0, prob0, treat_spec, par, X
    )
  }
  
  return(results)
}

# Calculate counterfactual effects
calculate_counterfactual_effects <- function(results_baseline, results_counterfactual) {
  effects <- list()
  
  for (year in names(results_baseline)) {
    # Direct effects
    direct_effect <- results_counterfactual[[year]]$prob - results_baseline[[year]]$prob
    
    # Indirect effects (spillover effects)
    indirect_effect <- mean(direct_effect, na.rm = TRUE)
    
    # Network statistics changes
    stats_baseline <- calculate_network_stats(results_baseline[[year]]$adj)
    stats_counterfactual <- calculate_network_stats(results_counterfactual[[year]]$adj)
    
    effects[[year]] <- list(
      direct_effects = direct_effect,
      indirect_effect = indirect_effect,
      network_stats_change = stats_counterfactual - stats_baseline,
      convergence = results_counterfactual[[year]]$converged,
      iterations = results_counterfactual[[year]]$iterations
    )
  }
  
  return(effects)
}

# Calculate network statistics
calculate_network_stats <- function(adj) {
  net <- network(adj, directed = FALSE)
  c(
    edges = network.edgecount(net),
    triangles = summary(net ~ triangles),
    density = network.density(net),
    reciprocity = if(is.directed(net)) summary(net ~ mutual) else NA
  )
}

# Save results to Excel
save_results_to_excel <- function(results, effects, filename) {
  wb <- createWorkbook()
  
  # Save networks for each year
  for (year in names(results)) {
    addWorksheet(wb, paste("Network", year))
    writeData(wb, paste("Network", year), results[[year]]$adj, rowNames = TRUE)
    
    addWorksheet(wb, paste("Probability", year))
    writeData(wb, paste("Probability", year), results[[year]]$prob, rowNames = TRUE)
  }
  
  # Save effect analysis
  addWorksheet(wb, "Effects Summary")
  effects_summary <- do.call(rbind, lapply(effects, function(x) {
    data.frame(
      indirect_effect = x$indirect_effect,
      edges_change = x$network_stats_change["edges"],
      triangles_change = x$network_stats_change["triangles"],
      converged = x$convergence,
      iterations = x$iterations
    )
  }))
  writeData(wb, "Effects Summary", effects_summary, rowNames = TRUE)
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# Main execution function
main_counterfactual_analysis <- function(config) {
  cat("Starting counterfactual simulation analysis...\n")
  
  # Prepare data
  cat("1. Preparing data...\n")
  data <- prepare_network_data(config$years, config$country_names)
  
  # Estimate model
  cat("2. Estimating TERGM model...\n")
  tergm_fit <- estimate_tergm_model(
    data$networks, 
    data$covariates, 
    data$attributes
  )
  
  # Setup counterfactual parameters
  cat("3. Setting up counterfactual parameters...\n")
  par <- setup_counterfactual_params(
    coef(tergm_fit),
    list(
      n_nodes = length(config$country_names),
      gw_decay = config$gw_decay,
      thresh = config$thresh,
      iter_max = config$iter_max
    )
  )
  
  # Define treatment effects
  treat_spec <- lapply(config$treatment_edges, function(edge) {
    list(
      i = which(config$country_names == edge$from),
      j = which(config$country_names == edge$to),
      new_value = edge$new_value
    )
  })
  
  # Run baseline analysis
  cat("4. Running baseline analysis...\n")
  baseline_results <- run_multiperiod_counterfactual(
    data, par, list(), config$years  # Empty treatment list = baseline
  )
  
  # Run counterfactual analysis
  cat("5. Running counterfactual analysis...\n")
  counterfactual_results <- run_multiperiod_counterfactual(
    data, par, treat_spec, config$years
  )
  
  # Calculate effects
  cat("6. Calculating counterfactual effects...\n")
  effects <- calculate_counterfactual_effects(baseline_results, counterfactual_results)
  
  # Save results
  cat("7. Saving results...\n")
  save_results_to_excel(counterfactual_results, effects, config$output_file)
  
  return(list(
    baseline = baseline_results,
    counterfactual = counterfactual_results,
    effects = effects,
    model = tergm_fit
  ))
}

# Configuration parameters
config <- list(
  years = 2000:2019,
  country_names = c("CHN", "USA", "RUS", "DEU", "FRA", "GBR", "JPN", "IND"),
  treatment_edges = list(
    list(from = "CHN", to = "USA", new_value = 0),  # Remove sanction
    list(from = "RUS", to = "USA", new_value = 0)   # Remove sanction
  ),
  gw_decay = 0.25,
  thresh = 0.20,
  iter_max = 20,
  output_file = "counterfactual_results.xlsx"
)

# Execute analysis
results <- main_counterfactual_analysis(config)

# Visualize counterfactual effects
visualize_counterfactual_effects <- function(results) {
  # Extract effect data
  effect_data <- do.call(rbind, lapply(names(results$effects), function(year) {
    data.frame(
      year = as.numeric(year),
      indirect_effect = results$effects[[year]]$indirect_effect,
      edges_change = results$effects[[year]]$network_stats_change["edges"]
    )
  })
  
  # Plot indirect effect trends
  plot(effect_data$year, effect_data$indirect_effect, 
       type = "b", col = "blue", lwd = 2,
       xlab = "Year", ylab = "Indirect Effect",
       main = "Counterfactual Indirect Effects Over Time")
  
  # Plot network size changes
  plot(effect_data$year, effect_data$edges_change,
       type = "b", col = "red", lwd = 2,
       xlab = "Year", ylab = "Change in Number of Edges",
       main = "Network Size Changes Under Counterfactual")
}

# Generate visualizations
visualize_counterfactual_effects(results)

cat("Counterfactual simulation analysis completed!\n")

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