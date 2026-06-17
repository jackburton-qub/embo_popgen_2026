### Transition probabilities under the Wright-Fisher model

# Q1)**. Let’s warm up by refreshing the Wright-Fisher model.

# a)** What is the average number of offspring per individual?


# Since the population size is constant, each individuals has on average 1
# offspring.


# b)** What is the probability that an individual has no offspring at
# all? Does that depend on the population size?

# This is happening when each of the 2N individuals of the next generation 
# picks another individual as its parent. Hence,
# P(nooffspring) = (1-1/2N)^(2N), which
# does depend on the population size but becomes approximately 0.37 if
# 2*N* is large.


# Q2)** Consider a population of size 2N = 100 and a current allele
# frequency of f = 10/100 = 0.1.

# a)** What is the probability that the allele frequency remains
# unchanged after one generation of random mating?


# Under the Wright-Fisher Model, transition probabilities are binomial.
# Hence, P(n|2N,f) = 2Nnf^n(1-f)^{2N}). In R, we can get that
# probability using the function `dbinom()`:

``` r
dbinom(x=10, size = 100, prob = 0.1)
```
# x = number of alleles in the next generation
# size = 2N
# probability = current allele frequency
# answer = 0.1318653


# b)** What is the probability that the allele is lost after one
# generation of random mating?


Again using `dbinom()`:

``` r
dbinom(x=0, size = 100, prob = 0.1)
```
# x = number of alleles in the next generation
# size = 2N
# probability = current allele frequency
#answer = 2.65614e-05



# c)** What is the probability that the allele is fixed after one
# generation of random mating?

Again using `dbinom()`:

``` r
dbinom(x=100, size = 100, prob = 0.1)
```
# x = number of alleles in the next generation
# size = 2N
# probability = current allele frequency
# answer = 1e-100


# d)** What is the probability that the allele frequency decreases after
# one generation of random mating?


# To calculate the probability of a decrease, we need to integrate the
# transition probability from 0 to n-1. In R, we can use the function
#`pbinom()`:

``` r
pbinom(9, size = 100, prob = 0.1)
```
# siz
#answer = 0.4512902


# d)** What is the probability that the allele frequency increases after
# one generation of random mating?

# We can calculate the probability of increase analogously as the integral
# from 201 to 1000, or as 1 minus the integral from 0 to 200:

``` r
1 - pbinom(10, size = 100, prob = 0.1)
```

# answer = 0.4168445

# Note that the probabilities for an increase and decrease are not equal,
# i.e. this distribution is not symmetric.


# Q3)** Plot the probability distribution on the allele frequency for a
# current allele frequency *f* = 0.1 and 2*N* = 10, 2*N* = 100,
# 2*N* = 1000 and 2*N* = 10000. What is the effect of population size?

# plotting just 2N = 10000

p <- 0:10000
d <- dbinom(p, size = 10000, prob = 0.1)

plot(p, d, type = 'b')

# plotting just 2N = 1000

p <- 0:1000
d <- dbinom(p, size = 1000, prob = 0.1)
plot (p, d, type = 'b')

##etc

# to write a function that runs all these effective population sizes at some time
# first write a par function for plotting grid style

par(mfrow=c(2,2))

# then write the function, first categorising twoN as the possible effective population
# sizes 10, 100, 1000, 10000. Then using the same code as above to carry out the dbinom and plotting

for(twoN in c(10, 100, 1000, 10000)){
  p <- 0:twoN # 
  d <- dbinom(p, size = twoN, prob = 0.1)
  plot(p, d, type = 'b')
}
# this function plots the probability distribution of allele frequencies for different population sizes. 
# As the population size increases, the distribution becomes more peaked around the initial allele frequency, 
# indicating less fluctuation in allele frequencies due to genetic drift.

### Simulating neutral allele frequency trajectories

# Q4)** Write a function `simulateWF()` to simulate allele trajectories
# under the Wright-Fisher model. Your function should take as input i) the
# population size 2*N*, ii) the initial allele frequency *f* and iii) the
# number *G* of generations to simulate. It should then return the allele
# frequency (between 0 and 1) for each generation as a vector.

# We can use the function `rbinom()` to simulate trajectories:

# to write a function using simulateWF(), we need to include population size 2N (e.g. 100)
# initial allele frequency f (e.g. 0.1) and number of generations G (e.g. 1000)

simulateWF <- function(twoN, f, G){
  p <- numeric(G + 1) # create a vector of length G+1 to store allele frequencies
  p[1] <- f # set the initial allele frequency
  for(i in 1:G){ # loop over generations
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN # simulate the next generation's allele frequency using binomial sampling
  }
  return(p) # return the vector of allele frequencies
}


# Q5)** Use your function `simulateWF()` to simulate 1000 trajectories
# with 2*N* = 100 and *f* = 0.1 for *G* = 1000 and plot them in one plot.
# In how many cases was the allele lost? Does this match the expectation?
# Repeat for different population sizes and initial allele frequencies.
# What is the effect of the population size?

trajectories <- replicate(1000, simulateWF(twoN = 100, f = 0.1, G = 1000)) 
# this provides a matrix of 1000 columns, each column being a trajectory of allele frequencies over 1000 generations
# it uses the simulateWF function by filling in the necessary parameters

plot(0, type = 'n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
# this provides an empty plot with the correct x and y limits to plot the trajectories

invisible(apply(trajectories, 2, lines, type='l'))
# this provides the lines for each trajectory in the plot, using the apply function 
# to apply the lines function to each column of the trajectories matrix

print(paste("Allele was lost in", sum(trajectories[1000,] == 0), "/", ncol(trajectories), "cases."))
# this provides the number of cases where the allele was lost 
# (i.e. the final allele frequency is 0) out of the total number of trajectories simulated

# [1] "Allele was lost in 907 / 1000 cases."

# The allele is expected to go to fixation with the probability of its
# initial allele frequency *f*, and hence to be lost with probability
# 1 − *f*.

# If you use larger population sizes, you should observe less fluctuation
# and eventually many sites will remain polymorphic even after *G* = 1000
# generations. However, if run for enough generations, the same fraction
# of alleles will be lost.


# Q6)** Use your function `simulateWF()` to study fixation probability
# of a new mutation under different population sizes (ensure *G* is large
# enough). Does the fixation probability depend on the population size?
# How does this affect the substitution rate?


# To study the fixation probability of a new mutation, set
# f= 1/2N. or 1 / twoN

trajectories <- replicate(1000, simulateWF(twoN = 100, f = 1/100, G = 1000))
# this changes population size to 100 and initial allele frequency to 1/100, 
# simulating 1000 trajectories over 1000 generations
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
# same plot function as previous 
invisible(apply(trajectories, 2, lines, type='l'))
# same application of trajectories as previous code

print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

## [1] "Allele was fixed in 11 / 1000 cases."

# Doing another example with an effective population size of 1000 (twoN = 1000) 
# and initial allele frequency of 1/1000:

trajectories <- replicate(1000, simulateWF(twoN = 1000, f = 1/1000, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

# [1] "Allele was fixed in 0 / 1000 cases."

# As population size increases, the probability of the mutated allele being fixed decreases

# The fixation probability of any mutation is given by its current
# frequency *f* and is hence markedly different between small and large
# populations. The substitution rate, however, is not since there are
# 2*N**μ* mutations occurring per generation, and hence many more in large
# than small populations.


### Simulating selection and genetic drift
<br/>

### Simulating selection and genetic drift

**Q7)** Write a function `simulateWFWithSelection()` that simulates
allele trajectories under both genetic drift and viability selection.
Similar to your function `simulateWF()` you wrote above, it should take
as input i) the population size 2*N*, ii) the initial allele frequency
*f*, iii) the number *G* of generations to simulate and iv) also a
vector *v* of viabilities for the genotypes AA, Aa and aa. It should
then return the allele frequency (between 0 and 1) for each generation
as a vector. In each generation, your function should apply selection to
alter the allele frequency, and then use binomial sampling to simulate
genetic drift in that modified allele frequency.

<details>
<summary>
Solutions
</summary>

We can use the function `rbinom()` to simulate trajectories:

``` r
simulateWFWithSelection <- function(twoN, f, G, v){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    # selection
    fA <- p[i]
    fa <- 1-fA
    fPrime <- (v[1]*fA*fA + v[2]*fA*fa)/(v[1]*fA*fA + v[2]*2*fA*fa + v[3]*fa*fa);
    # drift
    p[i+1] <- rbinom(1, size = twoN, prob = fPrime) / twoN
  }
  return(p)
}
```

</details>

<br/>

**Q8)** Use your function `simulateWFWithSelection()` to simulate 1000
trajectories with 2*N* = 100, and *f* = 0.1 and genic selection with
viabilities *v* = (1, 1 − *s*, (1 − *s*)<sup>2</sup>) for *G* = 1000 and
plot them in one plot. Start with *s* = 0.01. How strong does the
selection coefficient *s* have to be to see differences to the neutral
case simulated above (or simulated with `simulateWFWithSelection()` when
setting *v* = (1, 1, 1))? How is this affected by the population size
2*N*?

<details>
<summary>
Solutions
</summary>

``` r
s <- 0.01
trajectories <- replicate(100, simulateWFWithSelection(twoN = 1000, f = 0.1, G = 1000, v=c(1,1-s,(1-s)^2)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was fixed in 86 / 100 cases."

A weak selection coefficient of *s* = 0.01 leads to only very marginal
difference to the neutral case (allele is lost in about 80% rather than
90% of the cases). With *s* = 0.1, the allele goes to fixation in the
majority of cases. When using 2*N* = 1000, the allele goes to fixation
in about 90% of the cases even with *s* = 0.01.

</details>

<br/>

**Q9)** Set 2*N* = 10<sup>6</sup> and dominant *v* = (1, 1, 1 − *s*)
with *s* = 0.05. How often does the allele go to fixation within
*G* = 1000? And in case of a smaller population size? How does it look
like for the recessive case?

<details>
<summary>
Solutions
</summary>

``` r
s <- 0.05
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10^6, f = 0.1, G = 1000, v=c(1,1,1-s)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was fixed in 0 / 100 cases."

Despite a large population 2*N* = 10<sup>6</sup> and strong selection
*s* = 0.05, the allele is essentially never fixed. This is because there
is no selection benefiting the heterozygous Aa over the homozygous AA
genotype and drift is too weak to push the allele to fixation. When
using a smaller population size such as 2*N* = 10<sup>3</sup> the allele
is essentially always fixed because selection is strong enough to push
it quickly to high frequencies, and drift strong enough to fix it by
random fluctuations.

In the recessive case *v* = *c*(1, 1 − *s*, 1 − *s*) with *s* = 0.05 and
2*N* = 10<sup>6</sup>, the allele is always fixed. When using a smaller
population size (e.g. 2*N* = 1<sup>3</sup>), the allele is occasionally
lost early on because of drift (random fluctuations) before selection
could push the allele to a high enough frequency, also because there is
no selection benefiting the heterozygous Aa over the homozygous aa
genotype.
</details>

<br/>
