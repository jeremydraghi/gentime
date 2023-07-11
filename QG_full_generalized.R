library(doParallel)
library("flock")

floatMatch = function(x,y)
{
  return(abs(x - y) < 1e-9)
}

# Breaks integer into cleaned-up bit vector of length L, most sig. bit first
bits = function(x)
{
  return( as.integer(rev(intToBits(x)[1:L])) )
}

# Map maladaptation onto traits
tMap = function(z, env, effectOnDeath, selectionStrength)
{
  if(pFunction == "quad")
  {
    penalty = 1 + (z - env)^2 / selectionStrength
    return(c(1 / penalty^(1 - effectOnDeath), penalty^effectOnDeath))
  }
  if(pFunction == "abs")
  {
    penalty = abs(z - env) / selectionStrength
    return(c(1 / exp((1 - effectOnDeath)*penalty), exp(effectOnDeath * penalty)))
  }
}

pFunction = "abs"

#Environmental variation
VE = 1
sigmaE = sqrt(VE)

# Window size to evaluate population recovery
coda = 20

# Number of time steps before environmental change begins
burnIn = 0

# Total scale of environmental change
upperLimit = 10

# Minimum number of observations of both extinction and survival
minSamples = 100

# Number of loci per block
L = 8
# Number of blocks, each of L loci
blocks = 125
# Total genotype space size per block
size = 2^L
# Total number of loci
loci = L * blocks

# Convenient vectors
indexer = 1:blocks
bAndD = c(2,3)

# Pre-compute lookup tables for mutations and free recombination within a block.
powers = 2^((L-1):0)
bitMap = matrix(0, nrow=size, ncol=L)
for(i in 1:size)
{
  bitMap[i,] = bits(i-1)
}

recMaps = array(0, c(size, size, size))
mutMaps = matrix(0, ncol=L, nrow=size)
for(i in 1:size)
{
  vec = bitMap[i,]
  for(j in 1:L)
  {
    mutant = vec
    mutant[j] = (mutant[j] + 1) %% 2
    mutMaps[i,j] = sum(mutant * powers) + 1
  }
  for(j in 1:size)
  {
    parents = rbind(bitMap[i,], bitMap[j,])
    for(k in 1:size)
    {
      val = sum(parents[cbind(bitMap[k,]+1, 1:L)] * powers) + 1
      recMaps[i,j,k] = val
    }
  }
}

setwd("/Users/draghi/My Drive/birth_versus_death/QG_outputs")
project = "blue"
setwd(paste0("./",project))
tableFile = paste0(project, "_table.txt")
dt = read.table(tableFile, header=TRUE)
replicates = dim(dt)[1]
startingN = dt$N[1]

#deathLevels = c(0.05, 0.1, 0.2, 0.3, 0.5)
deathLevels = c(0.1)
# omega--fractional effect of maladaptation on death rates
betaLevels = seq(from = 0, to = 1, length.out=9)
phiLevels = seq(from = 0, to = 1, length.out=9)
# Per haploid genome mutation rate
muLevels = c(0.01)
parameterNames = c("d0", "beta", "phi", "mu")
nPara = length(parameterNames)
vals = list(deathLevels, betaLevels, phiLevels, muLevels)
levels = unlist(lapply(vals, length))
nTreatments = prod(levels)
treats = data.frame(matrix(ncol=nPara, nrow=nTreatments))
colnames(treats) = parameterNames
for(i in 1:nPara)
{
  previous = 1
  if(i > 1) previous = prod(levels[1:(i-1)])
  factor = levels[i] * previous
  treats[,i] = rep(rep(vals[[i]], times = previous), each = nTreatments / factor)
}

treats = treats[which(treats$beta + treats$phi == 1),]
nTreatments = nrow(treats)

mode = "evolving"

seeds = sample((1):1e8, nTreatments, replace=FALSE)

hubC = makeForkCluster(9)
registerDoParallel(hubC)

temp = foreach(r=1:nTreatments) %dopar%
{
  # Load parameters from treats and tableFile
  set.seed(seeds[r])
  mu = treats$mu[r]
  d0 = treats$d0[r]
  
  NAt5 = 0.1
  fXAt5 = (1 - NAt5) / 0.1
  if(pFunction == "abs")
  {
  	selectionStrength = 5 / log(fXAt5)
  }
  if(pFunction == "quad")
  {
    selectionStrength = 5^2 / (fXAt5 - 1)
  }
  
  phi = treats$phi[r]
  beta = treats$beta[r]
  dCoefficient = d0^phi
  M = ceiling(startingN / (1 - d0))
  pop = array(1, c(M, blocks, 2))
  child = matrix(0, nrow=2, ncol=blocks)
  traits = matrix(0, ncol=3, nrow=M)
  outfile = paste0(project,"_resC_", pFunction, ".txt")
  if(mode == "counting") countfile = paste0(project,"_counts_", pFunction, ".txt")
  
  SEARCHING = TRUE
  while(SEARCHING)
  {
    if(file.exists(outfile))
    {
      data = read.table(outfile, header=FALSE, col.names=c(parameterNames, "rate", "outcome"))
      data = data[which(floatMatch(data$d0, treats$d0[r]) & floatMatch(data$beta, treats$beta[r]) & floatMatch(data$phi, treats$phi[r]) & floatMatch(data$mu, treats$mu[r])),]
      survived = data$outcome > 0
      if(sum(survived) < 1 | sum(survived == 0) < 1)
      {
        if(dim(data)[1] == 0)
        {
          optRate = 0.05
          if(mode == "jumping") optRate = 10
        } else {
          if(sum(survived) < 1)
          {
            optRate = min(data$rate) * 0.5
          } else {
            optRate = max(data$rate) * 2
          }
        }
      } else {
        mod = suppressWarnings(glm(survived ~ data$rate, family=binomial(link = "logit")))
        guess = as.numeric(-mod$coefficients[1]/mod$coefficients[2])
        
        if(mode == "counting")
        {
          optRate = guess
          if(file.exists(countfile) == TRUE)
          {
            counts = read.table(countfile, header=FALSE, col.names=c(parameterNames, "rate", "N", "z", "bg", "dg"))
            counts = counts[which(floatMatch(counts$beta, treats$beta[r]) & floatMatch(counts$phi, treats$phi[r])),]
            if(sum(counts$N > 0) >= minSamples) SEARCHING = FALSE
          }
        } else {
          
          if(guess < 0) guess = 0
          optRate = 0
          hits = min(c(sum(survived == 1), sum(survived == 0)))
          while(optRate <= 0.002) optRate = rnorm(1, guess, guess/(0.5*sqrt(hits)))
          if(sum(survived) >= minSamples & sum(survived == 0) >= minSamples) SEARCHING = FALSE
        }
      }
    } else {
      optRate = 0.05
      if(mode == "jumping") optRate = 10
    }
    id = sample(1:replicates, 1)
    # Load population and allele effects
    dp = read.table(paste0(project,"_pop_", formatC(id, width=4, format="d", flag="0"), ".txt"), header=FALSE)
    de = read.table(paste0(project,"_effects_", formatC(id, width=4, format="d", flag="0"), ".txt"), header=FALSE)
    effects = matrix(c(de[,1], de[,2]), ncol=2, nrow=loci)
    startingPop = array(1, c(startingN, blocks, 2))
    startingPop[,,1] = as.matrix(dp[,2:(blocks+1)])
    startingPop[,,2] = as.matrix(dp[,(blocks+2):(2*blocks+1)])
    startingTraits = dp[,1]
    
    # Map phenotype effects
    phenoMaps = matrix(0, ncol=blocks, nrow=size)
    for(i in 1:size)
    {
      vec = bitMap[i,]
      for(j in 1:blocks)
      {
        phenoMaps[i,j] = sum(effects[cbind((1 + (j-1)*L):(j*L), vec+1)])
      }
    }
    
    N = startingN
    pop[1:N,,] = startingPop
    traits[1:N,1] = startingTraits
    env = 0

    for(i in 1:N)
    {
      traits[i,bAndD] = tMap(traits[i,1], env, beta, selectionStrength)
    }
    testPeriod = burnIn + upperLimit / optRate + coda
    
    countRunning = FALSE
    if(mode == "counting") countRunning = TRUE
    deathN = N
    birthN = N
    deathCount = 0
    birthCount = 0
    deathG = 0
    birthG = 0
    
    t = 0
    nextT = 1
    RUNNING = TRUE
    totalBirth = sum(traits[1:N,2])
    totalDeath = sum(traits[1:N,3])
    
    popSizes = NULL
    # Optional plotting variables
    # meanZs = NULL
    # envs = NULL
    
    while(RUNNING)
    {
      if(t >= nextT)
      {
        if(nextT > testPeriod)
        {
          ts = nextT - (coda-1):0
          pops = popSizes[(nextT-coda+1):nextT]
          slope = lm(pops ~ ts)$coefficients[2]
          if(slope > 0 & N > startingN / 10) RUNNING = FALSE
        }
        
        if(mode == "jumping")
        {
          if(nextT == burnIn)
          {
            env = optRate
            for(i in 1:N)
            {
              traits[i,bAndD] = tMap(traits[i,1], env, beta, selectionStrength)
            }
          }
        } else if(nextT >= burnIn & env < upperLimit) {
          env = env + optRate
          if(env > upperLimit) env = upperLimit
          for(i in 1:N)
          {
            traits[i,bAndD] = tMap(traits[i,1], env, beta, selectionStrength)
          }
        }
        
        if(env >= upperLimit & countRunning == TRUE)
        {
          countRunning = FALSE
          deathG = deathG + deathCount / deathN
          birthG = birthG + birthCount / birthN
          meanZ = mean(traits[1:N,1])
          
        }
        
        popSizes = c(popSizes, N)
        # meanZs = c(meanZs, mean(traits[1:N, 1]))
        # envs = c(envs, env)
        nextT = nextT + 1
        totalBirth = sum(traits[1:N,2])
        totalDeath = sum(traits[1:N,3])
      }
      
      effectiveBirth = dCoefficient * totalBirth  * (1 - N/M)^(1-phi)
      effectiveDeath = dCoefficient * totalDeath * d0 / (1 - N/M)^phi

      gTotal = effectiveBirth + effectiveDeath
      t = t + rexp(1, gTotal)
      event = rbinom(1, 1, effectiveBirth / gTotal)
      if(event == 0)
      {
        # Death event
        dead = sample(1:N, 1, prob=traits[1:N,3])
        totalBirth = totalBirth - traits[dead,2]
        totalDeath = totalDeath - traits[dead,3]
        # Replace dead individual with Nth individual
        if(dead < N)
        {
          pop[dead,,] = pop[N,,]
          traits[dead,] = traits[N,]
        }
        N = N - 1
        if(countRunning == TRUE)
        {
          deathCount = deathCount + 1
          if(deathCount == deathN)
          {
            deathG = deathG + 1
            deathCount = 0
            deathN = N
          }
        }
      } else {
        # Birth event
        parents = sample(1:N, 2, replace=FALSE, prob=traits[1:N,2])
        N = N + 1
        if(countRunning == TRUE)
        {
          birthCount = birthCount + 1
          if(birthCount == birthN)
          {
            birthG = birthG + 1
            birthCount = 0
            birthN = N
          }
        }
        recs = sample(1:size, 2 * blocks, replace=TRUE)
        for(j in 1:2)
        {
          pop[N,,j] = recMaps[cbind(pop[parents[j],,1], pop[parents[j],,2], recs[(1+blocks * (j-1)):(blocks * j)])]
        }
        
        nMut = rpois(1, mu)
        if(nMut > 0)
        {
          sites = sample(1:blocks, nMut, replace=TRUE)
          chromosomes = rbinom(nMut, 1, 0.5) + 1
          for(j in 1:nMut)
          {
            pop[N,sites[j], chromosomes[j]] = mutMaps[pop[N,sites[j], chromosomes[j]], sample(1:L, 1)]
          }
        }
        
        traits[N,1] = sum(phenoMaps[rbind(cbind(pop[N,,1], indexer), cbind(pop[N,,2], indexer))]) + rnorm(1, 0, sigmaE)
        traits[N,bAndD] = tMap(traits[N,1], env, beta, selectionStrength)
        totalBirth = totalBirth + traits[N,2]
        totalDeath = totalDeath + traits[N,3]
      }
      if(N <= 1)
      {
        RUNNING = FALSE
        if(N == 1) N = 0
      }
    }
    
    fl = lock(outfile)
    if(mode == "counting")
    {
      write.table(cbind(treats[r,], optRate, N, meanZ, birthG, deathG), countfile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    } else {
      write.table(cbind(treats[r,], optRate, N), outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    }
    unlock(fl)
  }
}  
stopCluster(hubC)

#write.table(cbind(1:length(envs), envs, meanZs, popSizes), "survived_example_03.txt", col.names = c("ts", "env", "meanZ", "N"), quote=FALSE, row.names=FALSE)



