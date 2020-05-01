#########################################################
### Simulation of social interactions with infections ###
### Supplementary to: "Efficient social distancing    ###  
### strategies to reduce the spread of COVID"         ###
### Authors: Marion Hoffman, Per Block                ###
#########################################################

library(igraph)
library(ggplot2)
library(Matrix)
source("functions_SN_covid_spread.R")


# Steps of the script:
# 1- generate network
# 2- Pre-calculate statistics
# 3- Calculate the entropy of the model with repetition strategy
# 4- Estimate parameters for the models with homophily and triadic strategies
# 5- Simulate random models and the three strategies in parallel


# 1- generate network

nAct <- 1000

net_sim <- socialCirclesNet(nActors = nAct, 
                            groupSize = 8, densityInGroups = 0.9,
                            minContactsGeo = 4, multMaxContactsGeo = 3, densityGeo = 0.3,
                            randomTiesPerPerson = 0.5,
                            minContactsHomo = 4, multMaxContactsHomo = 3, homoStrength = 2,
                            orderAndUnique = T)

edgelist <- net_sim[[1]]
graph <- graph_from_edgelist(edgelist, directed = F) # igra
network <- get.adjacency(graph)
attributes <- net_sim[[2]]

# 2- Pre-calculate statistics
triads <- calculate_triads(graph,network)
homophilies <- calculate_homophilies(network,attributes)


# 3- Calculate the entropy of the model with repetition strategy
params_R <- c(0,0,2.5,2) # params: 1: triad, 2: homophily, 3: repetition, 4: window repetition
num_simulations <- 50000
thining <- 2000
burnin <- 1000
entropies_base <- calculate_entropy_timedependent(network,
                                                  triads,
                                                  homophilies,
                                                  param_triad = 0,
                                                  param_homophily = 0,
                                                  param_repetition = 2.5,
                                                  window_repetition = 2,
                                                  num_simulations = num_simulations,
                                                  thining = thining, 
                                                  burnin = burnin, 
                                                  seed = 1)
entropy_base <- mean(entropies_base$entropies)
lastcontacts_base <- entropies_base$lastcontacts


# 4- Estimate parameters for the models with homophily and triadic strategies
pT <- optim(1, optim_entropy_T, 
            network=network, 
            triads=triads, 
            homophilies=homophilies, 
            entropy_base=entropy_base, 
            lower=0, upper=40, method="L-BFGS-B")
params_T <- c(pT$par,0,0,0)
pH <- optim(1, optim_entropy_H, 
            network=network, 
            triads=triads, 
            homophilies=homophilies, 
            entropy_base=entropy_base, 
            lower=0, upper=40, method="L-BFGS-B")
params_H <- c(0,pH$par,0,0)


# 5- Simulate random models and the three strategies in parallel
percentage_keepcontact <- 0.5
percentage_infection <- 0.8
window_exposure <- 1 * nAct
window_recovery <- 4 * nAct
starting_nodes <- sample(1:nAct, 1) # second number determines how many infected seeds
seed <- 56
burnin <- 10000

res_null <- simulate(network, triads, homophilies, percentage_keepcontact = 1, percentage_infection, 
                     c(0,0,0,0), window_exposure, window_recovery, starting_nodes, 
                     burnin = NULL, lastcontacts = NULL, seed, max_steps = Inf)

plot(unlist(res_null$infections))


res_random <- simulate(network, triads, homophilies, percentage_keepcontact, percentage_infection, 
                       c(0,0,0,0), window_exposure, window_recovery, starting_nodes, 
                       burnin = NULL, lastcontacts = NULL, seed, max_steps = Inf)

plot(unlist(res_random$infections))


res_T <- simulate(network, triads, homophilies, percentage_keepcontact, percentage_infection, 
                  params_T, window_exposure, window_recovery, starting_nodes, 
                  burnin = NULL, lastcontacts = NULL, seed, max_steps = Inf)

plot(unlist(res_T$infections))


res_H <- simulate(network, triads, homophilies, percentage_keepcontact, percentage_infection, 
                  params_H, window_exposure, window_recovery, starting_nodes, 
                  burnin = NULL, lastcontacts = NULL, seed, max_steps = Inf)

plot(unlist(res_H$infections))


res_R <- simulate(network, triads, homophilies, percentage_keepcontact, percentage_infection, 
                  params_R, window_exposure, window_recovery, starting_nodes, 
                  burnin = burnin, lastcontacts = lastcontacts_base, seed, max_steps = Inf)

plot(unlist(res_R$infections))


