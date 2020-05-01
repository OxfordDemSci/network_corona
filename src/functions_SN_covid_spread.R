#########################################################
### Simulation of social interactions with infections ###
### Supplementary to: "Efficient social distancing    ###  
### strategies to reduce the spread of COVID"         ###
### Authors: Marion Hoffman, Per Block                ###
#########################################################


############ PART 1 ################
#
# generation of social networks
#
####################################
#                                  #
#  Make some networks that look    #
#  like real ones                  #
#                                  #
#  - distribute people in 2D space #
#  - make people have ties with    #
#  others very close               #
#  - add a few random ties         #
#  - add another var and add homo  #
#  philous ties                    #
#  - give people random group      #
#  membership                      #
#  - connect people within groups  #
#                                  #
####################################


socialCirclesNet <- function(nActors = 1000, 
                             groupSize = 0, densityInGroups = 0.9,
                             minContactsGeo = 0, multMaxContactsGeo = 3, densityGeo = 0.3,
                             randomTiesPerPerson = 0,
                             minContactsHomo = 0, multMaxContactsHomo = 3, homoStrength = 5,
                             orderAndUnique = T){
  
  ##### generate groups in which people are highly connected #####
  if(groupSize>0){
    nGroups <- round(nActors/groupSize)
    membership <- sample(1:nGroups, nActors, replace = T)
    group_edges <- which((((outer(membership, membership, "==")*1) * (runif(nActors^2) < densityInGroups)) == 1), arr.ind = T)
  }
  
  ##### generate geo based edges #####
  # Distribute people on 2D space
  
  geo_x <- runif(nActors) * 100
  geo_y <- runif(nActors) * 100
  
  if(minContactsGeo > 0){
    # draw number of contacts that everybody gets on geo
    
    geo_friends <- floor(runif(nActors, min = minContactsGeo, max = (minContactsGeo*multMaxContactsGeo + 1)))
    
    geo_data <- data.frame(x = geo_x, y = geo_y, altN = geo_friends)
    
    geo_edges <- lapply(1:nActors, function(i) geo.to.net(geo_data, i, densityGeo))
  }
  
  ##### generate homophilous edges #####
  
  coVar <- runif(nActors)
  
  if(minContactsHomo > 0){
    
    # draw number of contacts that everybody gets on homo
    
    homo_friends <- floor(runif(nActors, min = minContactsHomo, max = (minContactsHomo*multMaxContactsHomo + 1)))
    
    homo_data <- data.frame(coVar = coVar, altN = homo_friends)
    
    homo_edges <- lapply(1:nActors, function(i) homo.to.net(homo_data, i, homoStrength))
  }
  
  if(randomTiesPerPerson > 0){
    ##### generate random long range ties #####
    
    random_edges <- cbind(sample(1:nActors, round(nActors * randomTiesPerPerson), replace = T), 
                          sample(1:nActors, round(nActors * randomTiesPerPerson), replace = T))
  }
  
  all_edges <- c(0,0)
  
  if(groupSize > 0){
    all_edges <- rbind(all_edges, group_edges)
  }
  
  if(randomTiesPerPerson > 0){
    all_edges <- rbind(all_edges, random_edges)
  }
  
  if(minContactsHomo > 0){
    all_edges <- rbind(all_edges, Reduce(rbind, homo_edges))
  }
  
  if(minContactsGeo > 0){
    all_edges <- rbind(all_edges, Reduce(rbind, geo_edges))
  }
  
  all_edges <- all_edges[-1,]
  
  row.names(all_edges) <- NULL
  colnames(all_edges) <- NULL
  
  if(orderAndUnique){
    all_edges <- all_edges[(all_edges[,1] != all_edges[,2]),]
    all_edges <- t(apply(all_edges, 1, sort))
    all_edges <- all_edges[order(all_edges[,1], all_edges[,2]),]
    all_edges <- unique(all_edges)
  }
  
  return(list(all_edges, coVar, geo_x))
}

##### functions used to make these networks #####

# function to take friends from geo data

geo.to.net <- function(dataframe, i, densityGeo){
  distanceToI <- sqrt((dataframe$x[i] - dataframe$x)^2 + (dataframe$y[i] - dataframe$y)^2)
  
  maxDist <- sort(distanceToI)[dataframe$altN[i]/densityGeo + 1]
  
  alters <- which(distanceToI <= maxDist)[!which(distanceToI <= maxDist) %in% i]
  
  alters <- sample(alters, dataframe$altN[i], replace = F)
  
  return(cbind(i, alters))
}


# function to take friends from homo data

homo.to.net <- function(dataframe, i, homoStrength){
  
  simToI <- (1 - abs(dataframe$coVar[i] - dataframe$coVar))^homoStrength
  
  simToI[i] <- 0
  
  sampleProb <- simToI / sum(simToI)
  
  alters <- sample(1:nrow(dataframe), dataframe$altN[i], replace = F, prob = sampleProb)
  
  return(cbind(i, alters))
}


############ PART 2 ################
#
# functions used as statistics in
# simulations
#
####################################

# function that returns the number of triads an edge is embedded in
calculate_triads <- function(graph,network){
  
  cocit <- cocitation(graph)
  triads <- cocit * network
  return(triads)
}

# function that returns the homophily of edges
calculate_homophilies <- function(network, attributes){
  
  range <- max(attributes) - min(attributes)
  
  sim <- 1 - (abs(outer(attributes, attributes, "-")) / range)
  
  homophilies <- sim * network
  
  return(homophilies)
}

calculate_homophilies_m <- function(network, attributes){
  
  diff1 <- abs(outer(attributes[,1], attributes[,1], "-"))
  diff2 <- abs(outer(attributes[,2], attributes[,2], "-"))
  
  diff <- sqrt((diff1 * diff1) + (diff2 * diff2))
  
  range <- max(diff) - min(diff)
  
  sim <- 1 - (diff / range)
  
  homophilies <- sim * network
  
  return(homophilies)
}

# pre-calculation of the link functions for the triad and homophily effects
calculate_fixed_probas <- function(network, triads, homophilies, param_triad, param_homophily) {
  
  probas <- param_triad*triads + param_homophily*homophilies
  probas <- exp(probas)
  probas[network == 0] <- 0
  return(probas)
  
}

############ PART 3 ################
#
# function to simulate infection
# spread
#
####################################


# Function to simulate contacts/infections
simulate_infection_curve <- function(network, 
                                     triads = NULL, 
                                     homophilies = NULL, 
                                     percentage_keepcontact = 1, # default: no contact reduction
                                     percentage_infection = 1, # default: deterministic infection
                                     param_repetition, 
                                     param_triad,
                                     param_homophily,
                                     window_repetition = 7, 
                                     window_exposed = 0, # default: no exposure time
                                     window_recovery = Inf, # default: no recovery
                                     starting_nodes,
                                     burnin = NULL,
                                     lastcontacts = NULL,
                                     seed = 1, 
                                     max_steps = Inf) # default: run the simulation until everyone is healthy or recovered
{
  
  set.seed(seed)
  
  n <- nrow(network)
  if(is.null(triads)) triads <- sparseMatrix(i=1:n,j=1:n)
  if(is.null(homophilies)) homophilies <- sparseMatrix(i=1:n,j=1:n)
  
  # 1 calculate probabilities beforehand, possibly start keeping track of last contacts
  probas_fixed <- calculate_fixed_probas(network, triads, homophilies, param_triad, param_homophily)
  probas <- probas_fixed / matrix(rep(rowSums(probas_fixed),n),nrow=n,ncol=n)
  if(param_repetition != 0 && is.null(lastcontacts)) lastcontacts <- vector(mode = "list", length = n)
  if(is.null(burnin)) burnin <- 0
  
  contacts <- list()
  infections <- list()
  step <- 1
  all_infected <- FALSE
  infected <- rep(-1,n)
  # Here
  infected[starting_nodes] <- window_exposed
  cpt_burnin <- 0
  
  while(!all_infected && step <= max_steps){
    
    keepcontact <- runif(1) <= percentage_keepcontact
    
    if(keepcontact) {
      
      # 2 pick a random node
      i <- sample(1:n,1)
      
      # 3 pick a neighbor
      neighbors <- which(network[i,] > 0)
      ps <- probas[i,neighbors]
      if(length(neighbors) > 1){
        j <- sample(neighbors, prob = ps, size=1)
      } else if(length(neighbors) == 1){
        j <- neighbors
      } else {
        next
      }
      contacts[[step]] <- c(i,j)
      
      # 4 infect and update the recovery
      if(cpt_burnin >= burnin) {
        newinfected <- infected
        
        # If i or j got infected and passed the exposure time and the other one was healthy, infect with a certain probability
        infect <- sample(c(F,T),1,prob=c(1-percentage_infection,percentage_infection))
        if(infected[i] >= window_exposed && infected[j] == -1 && infect) newinfected[j] <- 0
        if(infected[i] == -1 && infected[j] >= window_exposed && infect) newinfected[i] <- 0
        
        infected <- newinfected
      }
      
      # 5 update network of repetition
      if(param_repetition != 0) {
        lci <- lastcontacts[[i]]
        lci <- c(j,lci)
        lci <- lci[1:min(window_repetition,length(lci))]
        lastcontacts[[i]] <- lci
        
        lcj <- lastcontacts[[j]]
        lcj <- c(i,lcj)
        lcj <- lcj[1:min(window_repetition,length(lcj))]
        lastcontacts[[j]] <- lcj
      }
      
      # 6 update probas
      if(param_repetition != 0) {
        ps <- probas_fixed[i,neighbors]
        t <- table(lci)
        newterms <- t[match(neighbors, names(t))]
        newterms[is.na(newterms)] <- 0
        ps <- ps * exp(param_repetition*newterms)
        probas[i,neighbors] <- as.vector(ps / sum(ps))
        
        neighbors <- which(network[j,] > 0)
        ps <- probas_fixed[j,neighbors]
        t <- table(lcj)
        newterms <- t[match(neighbors, names(t))]
        newterms[is.na(newterms)] <- 0
        ps <- ps * exp(param_repetition*newterms)
        probas[j,neighbors] <- as.vector(ps / sum(ps))
      }
      
    } 
    
    # 7 update the infection counters
    if(cpt_burnin >= burnin) {
      infected[infected >= 0] <- infected[infected >= 0] + 1
      infections[[step]] <- sum(infected >= 0)
    }
    
    # 8 update recoveries
    if(cpt_burnin >= burnin) {
      infected[infected >= window_exposed + window_recovery] <- -2
    } 
    
    # check if half or all network is infected
    all_infected <- (sum(infected < 0) == n)
    step <- step + 1
    cpt_burnin <- cpt_burnin+1
  }
  
  return(list(contacts = contacts,
              infections = infections,
              number_steps_max = step-1,
              infected = infected))
  
}

# Umbrella function
simulate <- function(network, 
                     triads, 
                     homophilies, 
                     percentage_keepcontact,
                     percentage_infection,
                     params, 
                     window_exposure,
                     window_recovery,
                     starting_nodes,
                     burnin,
                     lastcontacts,
                     seed,
                     max_steps){
  res <- simulate_infection_curve(network, 
                                  triads, 
                                  homophilies, 
                                  percentage_keepcontact,
                                  percentage_infection,
                                  param_repetition = params[3], 
                                  param_triad = params[1],
                                  param_homophily = params[2],
                                  window_repetition = params[4], 
                                  window_exposure,
                                  window_recovery,
                                  starting_nodes,
                                  burnin,
                                  lastcontacts,
                                  seed,
                                  max_steps)
  return(res)
}


############ PART 4 ################
#
# functions related to entropy and
# estaimting proper parameter sizes
#
####################################

# Calculate entropy (R_H) with no time-dependent statistics
calculate_entropy_fixed <- function(network,
                                    triads,
                                    homophilies,
                                    param_triad,
                                    param_homophily) {
  
  n <- nrow(network)
  degrees <- rowSums(network)
  isolates <- degrees == 0
  if(is.null(triads)) triads <- sparseMatrix(i=1:n,j=1:n)
  if(is.null(homophilies)) homophilies <- sparseMatrix(i=1:n,j=1:n)
  
  probas <- calculate_fixed_probas(network, triads, homophilies, param_triad, param_homophily)
  probas <- probas / matrix(rep(rowSums(probas),n),nrow=n,ncol=n)
  fakeprobas <- probas
  fakeprobas[fakeprobas == 0] <- 1
  fakedegrees <- degrees
  fakedegrees[fakedegrees == 1] <- 2
  
  allRH <- 1 + rowSums(probas[!isolates,] * log2(fakeprobas[!isolates,])) / log2(fakedegrees[!isolates])
  RH <- mean(allRH)
  return(RH)
}

# Calculate entropy (R_H) with time-dependent statistics (repetition)
calculate_entropy_timedependent <- function(network,
                                            triads,
                                            homophilies,
                                            param_triad,
                                            param_homophily,
                                            param_repetition,
                                            window_repetition,
                                            num_simulations,
                                            burnin,
                                            thining,
                                            seed) {
  
  set.seed(seed)
  
  n <- nrow(network)
  degrees <- rowSums(network)
  isolates <- degrees == 0
  
  if(is.null(triads)) triads <- sparseMatrix(i=1:n,j=1:n)
  if(is.null(homophilies)) homophilies <- sparseMatrix(i=1:n,j=1:n)
  probas_fixed <- calculate_fixed_probas(network, triads, homophilies, param_triad, param_homophily)
  probas <- probas_fixed / matrix(rep(rowSums(probas_fixed),n),nrow=n,ncol=n)
  
  # initiate the last contacts
  lastcontacts <- vector(mode = "list", length = n)
  for(i in 1:n){
    neighbors <- which(network[i,] > 0)
    if(length(neighbors) > 1){
      lci <- c()
      for(k in 1:window_repetition){
        j <- sample(neighbors, size=1)
        lci <- c(lci,j)
      }
      lastcontacts[[i]] <- lci
    } else if(length(neighbors) == 1){
      lci <- rep(neighbors, window_repetition)
      lastcontacts[[i]] <- lci
    }
    ps <- probas_fixed[i,neighbors]
    t <- table(lci)
    newterms <- t[match(neighbors, names(t))]
    newterms[is.na(newterms)] <- 0
    ps <- ps * exp(param_repetition*newterms)
    probas[i,neighbors] <- as.vector(ps / sum(ps))
  }
  
  
  # simulate some steps and store all probas
  step <- 1
  cpt_burnin <- 0
  cpt_thining <- 0
  averageRH <- 0
  ntotal <- 0
  all <- c()
  
  while(step <= num_simulations){
    #print(step)
    
    # 2 pick a random node
    i <- sample(1:n,1)
    
    # 3 pick a neighbor
    neighbors <- which(network[i,] > 0)
    ps <- probas[i,neighbors]
    if(length(neighbors) > 1){
      j <- sample(neighbors, prob = ps, size=1)
    } else if(length(neighbors) == 1){
      j <- neighbors
    } else {
      next
    }
    
    # 6 update network of repetition
    lci <- lastcontacts[[i]]
    lci <- c(j,lci)
    lci <- lci[1:min(window_repetition,length(lci))]
    lastcontacts[[i]] <- lci
    
    lcj <- lastcontacts[[j]]
    lcj <- c(i,lcj)
    lcj <- lcj[1:min(window_repetition,length(lcj))]
    lastcontacts[[j]] <- lcj
    
    # 7 update probas
    ps <- probas_fixed[i,neighbors]
    t <- table(lci)
    newterms <- t[match(neighbors, names(t))]
    newterms[is.na(newterms)] <- 0
    ps <- ps * exp(param_repetition*newterms)
    probas[i,neighbors] <- as.vector(ps / sum(ps))
    
    neighbors <- which(network[j,] > 0)
    ps <- probas_fixed[j,neighbors]
    t <- table(lcj)
    newterms <- t[match(neighbors, names(t))]
    newterms[is.na(newterms)] <- 0
    ps <- ps * exp(param_repetition*newterms)
    probas[j,neighbors] <- as.vector(ps / sum(ps))
    
    
    # Entropy
    cpt_burnin <- cpt_burnin + 1
    if(cpt_burnin > burnin) cpt_thining <- cpt_thining + 1
    if(cpt_burnin > burnin && cpt_thining == thining){
      fakeprobas <- probas
      fakeprobas[fakeprobas == 0] <- 1
      fakedegrees <- degrees
      fakedegrees[fakedegrees == 1] <- 2
      allRH <- 1 + rowSums(probas[!isolates,] * log2(fakeprobas[!isolates,])) / log2(fakedegrees[!isolates])
      print(mean(allRH))
      averageRH <- averageRH +mean(allRH)
      
      ntotal <- ntotal + 1
      cpt_thining <- 0
      all <- c(all,mean(allRH))
    }
    
    step <- step + 1
  }
  
  #plot(all)
  
  return(list(entropies = all,
              lastcontacts = lastcontacts))
}


# Umbrella functions for the optimization
optim_entropy_T <- function(p1,network,triads,homophilies,entropy_base){
  newentropy <- calculate_entropy_fixed(network,
                                        triads,
                                        homophilies,
                                        param_triad = p1,
                                        param_homophily = 0)
  return(abs(newentropy - entropy_base))
}
optim_entropy_H <- function(p2,network,triads,homophilies,entropy_base){
  newentropy <- calculate_entropy_fixed(network,
                                        triads,
                                        homophilies,
                                        param_triad = 0,
                                        param_homophily = p2)
  return(abs(newentropy - entropy_base))
}


