# Infr
Infr provides a set of functions for performing Bayesian Inference on systems which can be described as Markov Chains. Given a model, a data series and a prior distribution it can compute estimates for the marginal and posterior probabilities for the parameters.

## Dependencies

* [lisp-stat](https://lisp-stat.dev)
* [iterate](https://iterate.common-lisp.dev)

## Usage
Infr operates on Markov Chain models which map data points to random variables. A model should be functions of the form `(f beta x) -> distribution`, where distribution is a lisp-stat distribution object.

Load the system with asdf.
```(asdf:load-system :infr)```

For testing purposes one can generate a vector containing a realisation of a process defined by the model function, f:
```(defvar y (generate-markov-chain f n :xi (x0 x1 ...)))```.

To compute the log posterior probability density for beta=(2.0) given the data sequence y and the model f, with a normal prior:
```(log-posterior f (vector (r-normal 0d0 0.5d0)) (vector 2d0) y)```.

## Methodology
Infr explots the fact the the joint probability of a realization of a markov chain is equal to the prodduct of the probability of each step.
𝑃 (𝑥𝑛, 𝑥𝑛−1, …𝑥0) = 𝑃 (𝑥𝑛|𝑥𝑛−1)𝑃 (𝑥𝑛−1 | 𝑥𝑛−2)…𝑃 (𝑥1|𝑥0)𝑃 (𝑥0)
This reduces computing the Bayesian Likelihood to a simple product.

The marginal likelihood is estimated by a one-step particle filter method. A sample of particles is drawn from the prior distribution and used to estimate the integral.

The method is explained fully in this [article](http://homocarbonis.xyz/infr/infr.html).

## Reference
### generate-markov-chain
generate-markov-chain f n &key (xi #(0d0))
=> vector
* f - Function (f beta x) => distribution
* n - integer
* xi - vector

Generate a realization of n elements of the Markov chain defined by f, beginning with xi.

### estimate-parameters
estimate-parameters f prior y &key (offset 2) (sample-count 100)
=> vector
* f - model function
* prior - distribution vector
* y - data vector
* offset - integer - Offset from start of y to allow for the 'memory' of the model.
* Number of samples taken from the prior distribution to estimate beta.

Find a Monte Carlo estimate of the parameter of a model.  We draw a large
number of samples from 'prior. Then for each sampled value of beta we find
P(y|model,beta). Finally we take the mean of the samples using P(y|model,beta)
as weights, to find an expected value for beta.

### log-marginal
log-marginal (model prior y &key (offset 2) (sample-count 100)
=> double-float

Estimate the log of the marginal probability for y. We take a random sample of beta and use it to estimate the integral of P(y|beta)P(beta) over beta.
The result can be used to find the Bayes factor to compare different models.

### log-posterior
log-posterior (model prior beta y &key (offset 2) (sample-count 100)
=> double-float
* beta - parameter vector

Calculate the posterior probability density of beta: log P(beta| model, y)
