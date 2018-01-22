load("genome1.rdata")
load("genome2.rdata")
load("genome3.rdata")

# Ex5 Ex6
genome1.subset = genome1[1:1000,]
cov.set = apply(genome1.subset[,2:5], 1, sum)
ref.subset = reference[1:1000]
genome1.subset = cbind(genome1.subset,cov.set)

cov = apply(genome1[,2:5], 1, sum)
genome1 = cbind(genome1,cov)
mean.coverage = mean(genome1[,6],2) # 20
max.coverage = max(genome1[,6],2)

plot(genome1.subset[1:1000,1],genome1.subset[1:1000,6],"l",col="blue", 
     main="Coverage for first 1000 position", 
     xlab = "Position", ylab = "Coverage")

# Ex7

genome1.subset$matches <- apply(genome1.subset, 1, function(row) row[ref.subset[row[1]]])

# Ex8

calculate.pvalue = function(genome,perror){
  genome.length = nrow(genome)
  pvalue = vector(length=genome.length, mode = "double")
  for (i in 1:genome.length) {
    if (genome[i,6] == 0) {
      pvalue[i] = 1
    } else {
      pvalue[i] = 
        binom.test(genome[i,6]-genome[i,7], genome[i,6], p = perror, 
                   alternative = "greater")$p.value
    }
  }
  return(pvalue)
  
}

pval = calculate.pvalue(genome1.subset,0.05)


# Ex9
genome1.subset = cbind(genome1.subset,pval)
genome1.subset.order=order(genome1.subset[,8], decreasing=FALSE)
genome1.subset.sorted = genome1.subset[genome1.subset.order,]

# combine into one function
sort.pval = function(genome){
  coverage = apply(genome[,2:5], 1, sum)
  genome = cbind(genome,coverage)
  genome$matches <- apply(genome, 1, function(row) row[reference[row[1]]])
  pval = calculate.pvalue(genome,0.05)
  genome = cbind(genome,pval)
  genome.order=order(genome[,8], decreasing=FALSE)
  genome.pval.sorted = genome[genome.order,]
  return(genome.pval.sorted)
}








