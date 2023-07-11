# Breaks integer into cleaned-up bit vector of length L, most sig. bit first
bits = function(x)
{
  return( as.integer(rev(intToBits(x)[1:L])) )
}

setwd("/Users/draghi/My Drive/birth_versus_death/QG_outputs/")

#Environmental variation
VE = 1
sigmaE = sqrt(VE)

# Number of loci per block
L = 8
# Number of blocks, each of L loci
blocks = 125
# Total genotype space size per block
size = 2^L
# Total number of loci
loci = L * blocks

# sd of Gaussian allele effect
alleleEffect = 0.035

# Convenient vectors
indexer = 1:blocks

# Number of samples at the end of "ancestral" runs.
hertSamples = 1e6

# How many populations to make
replicates = 100

targetVA = 0.35
tries = 1e6
pFixed = 0.5

# Population size at well-adapted equilibrium
N = 10000

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

project = "cerulean"
if(!dir.exists(project)) dir.create(project)
setwd(paste0("./",project))
tableFile = paste0(project, "_table.txt")
if(file.exists(tableFile))
{
  dt = read.table(tableFile, header=TRUE)
  lastID = max(dt$id)
} else {
  write.table(paste("id", "seed", "N", "VA", sep=" "), tableFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  lastID = 0
}

seeds = sample((1):1e8, replicates, replace=FALSE)

for(r in 1:replicates)
{
  set.seed(seeds[r])
  print(r)
  pop = array(1, c(N, blocks, 2))
  traits = rep(0, N)
  
  # Draw effects and map phenotypes
  effects = matrix(rnorm(2*loci, 0, alleleEffect), ncol=2, nrow=loci)
  AA = 2*effects[,1]
  AB = (effects[,1] + effects[,2])
  BB = 2*effects[,2]
  
  # Draw initial frequencies
  fixed = rbinom(1, loci, pFixed)
  probs = rep(0, loci)
  probs = rbinom(loci, 1, 0.5)
  probs[sample(1:loci, loci - fixed, replace=FALSE)] = runif(loci - fixed, 0, 1)
  
  means = AA * probs^2 + AB * 2 * probs * (1 - probs) + BB * (1-probs)^2
  meanZ = sum(means)
  VA = sum((AA - means)^2 * probs^2 + (AB - means)^2 * 2 * probs * (1 - probs) + (BB - means)^2 * (1-probs)^2)
  dist = (meanZ)^2 + 10 * (VA - targetVA)^2
  for(tr in 1:tries)
  {
    proposed = probs
    locus = sample(1:loci, 1)
    if(rbinom(1, 1, 1 - pFixed) == 1)
    {
      proposed[locus] = runif(1,0,1)
    } else {
      proposed[locus] = sample(c(0,1), 1)
    }
    means = AA * proposed^2 + AB * 2 * proposed * (1 - proposed) + BB * (1-proposed)^2
    newMeanZ = sum(means)
    newVA = sum((AA - means)^2 * proposed^2 + (AB - means)^2 * 2 * proposed * (1 - proposed) + (BB - means)^2 * (1-proposed)^2)
    newDist = (newMeanZ)^2 + 10 * (newVA - targetVA)^2
    if(newDist < dist)
    {
      #print(c(VA, newVA, meanZ, newMeanZ))
      VA = newVA
      meanZ = newMeanZ
      probs = proposed
      dist = newDist
    }
  }
  
  for(i in 1:blocks)
  {
    probMat = matrix(0, ncol=2, nrow=L)
    probMat[,1] = probs[(1 + (i-1)*L):(i*L)]
    probMat[,2] = 1 - probMat[,1]
    pickProbs = rep(0, size)
    for(j in 1:size)
    {
      pickProbs[j] = prod(probMat[cbind(1:L, bitMap[j,]+1)])
    }
    pop[1:N,i,] = sample(1:size, 2*N, replace=TRUE, prob=pickProbs)
  }
  phenoMaps = matrix(0, ncol=blocks, nrow=size)
  for(i in 1:size)
  {
    vec = bitMap[i,]
    for(j in 1:blocks)
    {
      phenoMaps[i,j] = sum(effects[cbind((1 + (j-1)*L):(j*L), vec+1)])
    }
  }
  for(i in 1:N)
  {
    traits[i] = sum(phenoMaps[rbind(cbind(pop[i,,1], indexer), cbind(pop[i,,2], indexer))]) + rnorm(1, 0, sigmaE)
  }
  child = matrix(0, nrow=2, ncol=blocks)
  midparent = rep(0, hertSamples)
  offspring = rep(0, hertSamples)
  allParents = matrix(sample(1:N, 2* hertSamples, replace=TRUE), ncol=2)
  deviates = rnorm(hertSamples, 0, sigmaE)
  for(i in 1:hertSamples)
  {
    parents = allParents[i,]
    while(parents[1] == parents[2]) parents[2] = sample(1:N, 1)
    midparent[i] = mean(traits[parents])
    recs = sample(1:size, 2 * blocks, replace=TRUE)
    for(j in 1:2)
    {
      child[j,] = recMaps[cbind(pop[parents[j],,1], pop[parents[j],,2], recs[(1+blocks * (j-1)):(blocks * j)])]
    }
    offspring[i] = sum(phenoMaps[rbind(cbind(child[1,], indexer), cbind(child[2,], indexer))]) + deviates[i]
  }
  mod = lm(offspring ~ midparent) 
  measuredVA = as.numeric(mod$coefficients[2] * var(offspring))
  
  write.table(cbind(r, seeds[r], N, measuredVA), tableFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  
  write.table(cbind(traits, pop[,,1], pop[,,2]), paste0(project,"_pop_", formatC(r, width=4, format="d", flag="0"), ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
  write.table(effects, paste0(project,"_effects_", formatC(r, width=4, format="d", flag="0"), ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
}  



