###  load the required packages  ###
library(data.table)
library(dplyr)
library(igraph)
library(Matrix)
library(ggplot2)

# read the g6 data
graph6.file <- 'graph9c.g6'
g6 <- readLines(con = graph6.file)

# initial parameters
R0 <- 2.8 # basic reproduction number
Time <- 10000 # total time steps (the process ends earlier)
DeltaT <- 0.1 # time steps increment
I0 <- 0.1 # proportion of initial individuals that are seeded as symptomatic infected
prop_Is <- 0.6 # proportion of infected individuals that are symptomatic
prop_D <- 0.0104 # proportion of a symptomatic infected to die
threshold <- 0.2 # proportion of the population that specifies if the outbreak size is big enough to be considered as valid
dur_Is <- 5/DeltaT # mean time for symptomatic infectious to recover
dur_Ia <- 2.9/DeltaT # mean time for asymptomatic infectious to recover


simu <- function(network){
  
  # network and its characteristics
  g <- network
  N <- length(V(g)) # number of nodes in the network
  deg <- degree(g) # degree of each node in the network
  d.m <- mean(deg) # average degree of the network
  
  #  beta <- beta_est_explicit(dur_Is=dur_Is, dur_Ia=dur_Ia, p=prop_Is, mixing_matrix=adj.mat, R0=R0)
  beta <- R0/(prop_Is*dur_Is + (1-prop_Is)*dur_Ia)
  
  # assign random times from gamma distributions
  sigma <- rgamma(n=N, shape=1.3521, scale=3.7719)/DeltaT # incubation period
  gamma_s <- rgamma(n=N, shape=4, scale=1.25)/DeltaT # recovery period for symptomatic
  gamma_a <- rgamma(n=N, shape=1, scale=1/0.35)/DeltaT  # recovery period for asymptomatic
  d_s <- rgamma(n=N, shape=4.9383, scale=3.6045)/DeltaT # time to death
  ind_death <- rbinom(n = N, size = 1, prob = prop_D)
  ind_recov <- 1 - ind_death
  
  #  SEIR simulation
  V1 <- data.table(id = c(1:N), S = 1, E = 0, Ia = 0, Is = 0,
                   R = 0, D = 0, rand = 0, tmp = 0, Ias = 0, 
                   firstE = 0, firstI = 0, firstR = 0, firstD = 0) %>% as.matrix()
  
  #  set initial spreader(s)
  spreaders <- sample(1:nrow(V1), round(I0*N))
  V1[spreaders, 'Is'] <- 1
  V1[spreaders, 'S'] <- 0
  V1[spreaders, 'firstI'] <- 1
  
  Res <- V1 %>% as.data.table %>% 
    .[, .(s = sum(S), e = sum(E), ia = sum(Ia), is = sum(Is), r = sum(R), d = sum(D),
          TOT = sum(S+E+Ia+Is+R+D))] %>% 
    .[, `:=`(t = 1)]
  
  valid_sim <- FALSE
  
  for(t in 2:Time){
    
    # TRANSMITTION PROCESSES
    
    #  exposed: S -> E
    #  link contagion
    E(g)$weight <- runif(length(E(g))) < 1-exp(-beta)
    V1[, 'Ias'] <- (V1[, 'Ia'] + V1[, 'Is']) >= 1  # infectious nodes
    if (is.null(E(g)$weight)==TRUE){
      A <- as.matrix(0, ncol=N, nrow=N)
    }else{
      A <- as_adjacency_matrix(g, attr = "weight")  # adjacency matrix of link infection
    }  
    
    tmp0 <- (A %*% V1[, 'Ias']) >= 1
    tmp1 <- tmp0 * V1[, 'S']
    V1[, 'tmp'] <- tmp1[, 1]
    V1[, 'E'] <- V1[, 'E'] + V1[, 'tmp']
    V1[, 'S'] <- V1[, 'S'] - V1[, 'tmp'] 
    V1[, 'Ias'] <- 0
    V1[, 'tmp'] <- 0
    V1 <- V1 %>% as.data.table %>% 
      .[, firstE := ifelse(firstE == 0 & E == 1, t, firstE)] %>% 
      as.matrix() #  first time exposed
    
    # infect: E -> I
    V1[, 'rand'] <- ifelse(V1[,'firstE'] == 0, 0, (t-V1[,'firstE']) >= sigma)
    V1[, 'tmp'] <- V1[, 'rand'] * V1[, 'E'] #  select nodes to change state
    idx.inf <- which(V1[, 'tmp'] == 1)
    if (length(idx.inf) >= 1){
      idx.Is <- rbinom(length(idx.inf),1,prop_Is) * idx.inf # each infected is symptomatic with probability prop_Is
      if(sum(idx.Is) == 0){idx.Is <- c()}
      idx.Ia <- idx.inf[! idx.inf %in% idx.Is]
      tmp.Is <- V1[, 'tmp']
      tmp.Is[idx.Ia] <- 0
      tmp.Ia <- V1[, 'tmp']
      tmp.Ia[idx.Is] <- 0
      if(!is.null(tmp.Ia)) {V1[, 'Ia'] <- V1[, 'Ia'] + tmp.Ia}
      if(!is.null(tmp.Is)) {V1[, 'Is'] <- V1[, 'Is'] + tmp.Is}
      V1[, 'E'] <- V1[, 'E'] - V1[, 'tmp']
      V1[, 'rand'] <- 0
      V1[, 'tmp'] <- 0
    }  
    V1 <- V1 %>% as.data.table %>% 
      .[, firstI := ifelse(firstI == 0 & (Ia == 1 | Is == 1), t, firstI)] %>% 
      as.matrix() #  first time infected
    
    # recover: Ia -> R
    V1[, 'rand'] <- ifelse(V1[,'firstI'] == 0, 0, (t-V1[,'firstI']) >= gamma_a)
    V1[, 'R'] <- V1[, 'R'] + V1[, 'Ia'] * V1[, 'rand']
    V1[, 'Ia'] <- V1[, 'Ia'] - V1[, 'Ia'] * V1[, 'rand']
    V1[, 'rand'] <- 0
    V1 <- V1 %>% as.data.table %>% 
      .[, firstR := ifelse(firstR == 0 & R == 1, t, firstR)] %>% 
      as.matrix() #  first time recovered
    
    # death: Is -> D
    V1[, 'rand'] <- ind_death * ifelse(V1[,'firstI'] == 0, 0, (t-V1[,'firstI']) >= d_s)
    V1[, 'D'] <- V1[, 'D'] + V1[, 'Is'] * V1[, 'rand']
    V1[, 'Is'] <- V1[, 'Is'] - V1[, 'Is'] * V1[, 'rand']
    V1[, 'rand'] <- 0
    V1 <- V1 %>% as.data.table %>% 
      .[, firstD := ifelse(firstD == 0 & D == 1, t, firstD)] %>% 
      as.matrix() #  first time dead
    
    # recover: Is -> R
    V1[, 'rand'] <- ind_recov * ifelse(V1[,'firstI'] == 0, 0, (t-V1[,'firstI']) >= gamma_s)
    V1[, 'R'] <- V1[, 'R'] + V1[, 'Is'] * V1[, 'rand']
    V1[, 'Is'] <- V1[, 'Is'] - V1[, 'Is'] * V1[, 'rand']
    V1[, 'rand'] <- 0
    V1 <- V1 %>% as.data.table %>% 
      .[, firstR := ifelse(firstR == 0 & R == 1, t, firstR)] %>% 
      as.matrix() #  first time recovered
    
    # summary stats
    V2 <- V1 %>% as.data.table %>% 
      .[, .(s = sum(S), e = sum(E), ia = sum(Ia), is = sum(Is), r = sum(R), d = sum(D), 
            TOT = sum(S+E+Ia+Is+R+D))] %>% 
      .[, `:=`(t = t)]
    Res <- Res %>% rbind(V2)
    
    #  terminate process if the epidemic disappears
    if((V2[1, 'e'] + V2[1, 'is'] + V2[1, 'ia']) == 0){
      break
    }
    
  }  
  
  # indicator for valid simulation
  if( sum(tail(Res,1)[,c('r','d')]) > round(threshold*N) ){
    valid_sim <- TRUE
  }
  
  return(list(Res, valid_sim))
  
}  


# run 1000 realisations for each network
for (netID in i:length(g6)){
  
  network <- rgraph6::igraph_from_text(g6[[netID]])
  
  res <- list()
  
  for (ii in c(1:1000)){  # ii is for simulation
    validity_msg <- FALSE
    while(validity_msg == FALSE){
      out <- simu(network[[1]])
      validity_msg <- out[[2]]
    }  
    res[[ii]] <- list(out)
    print(c(net=netID, sim=ii))
  }
  
  save(res, file = paste0('simOutcomes/net_', netID, '.RData'))
  
}
