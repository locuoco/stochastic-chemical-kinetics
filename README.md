# Stochastic chemical kinetics
Stochastic chemical kinetics using Gillespie algorithm (also stochastic simulation algorithm or SSA) and chemical master equation (CME), with application to enzyme kinetics.

Bibliography:
- D. T. Gillespie, *Exact Stochastic Simulation of Coupled Chemical Reactions*, The Journal of Physical Chemistry, Vol. 81, No. 25, 1977.
- L. A. Segel, M. Slemrod, *The quasi-steady-state assumption: a case study in perturbation*, SIAM Rev. 1989; 31(3):446-477.
- J. Borghans, R. de Boer, L. Segel, *Extending the quasi-steady state approximation by changing variables*, Bull. Math. Biol. 58 (1996) 43-63.
- R. A. Copeland, *ENZYMES: A Practical Introduction to Structure, Mechanism, and Data Analysis*, second edition, Wiley-VCH, 2000.
- A. R. Tzafriri, *Michaelis-Menten kinetics at high enzyme concentrations*, Bull. Math. Biol. 2003; 65(6):1111-1129.
- D. T. Gillespie, *Stochastic Simulation of Chemical Kinetics*, Annu. Rev. Phys. Chem. 2007, 58:35-55.
- M. G. Pedersen, A. M. Bersani, E. Bersani, G. Cortese, *The total quasi-steady-state approximation for complex enzyme reactions*, Mathematics and Computers in Simulation 79 (2008) 1010-1019.
- Y. Cao, D. C. Samuels, *Discrete Stochastic Simulation Methods for Chemically Reacting Systems*, Methods Enzymol. 2009; 454: 115-140.
- K. R. Sanft, D. T. Gillespie, L. R. Petzold, *Legitimacy of the stochastic Michaelis-Menten approximation*, IET Syst. Biol., 2011, Vol. 5, Iss. 1, pp. 58-69.
- V. Sunkara, *On the Properties of the Reaction Counts Chemical Master Equation*, Entropy 2019, 21(6), 607.
- J. K. Kim, J. J. Tyson, *Misuse of the Michaelis-Menten rate law for protein interaction networks and its remedy*, PLoS Computational Biology 16(10), 2020.
- T.-H. Ahn, Y. Cao, L. T. Watson, *Stochastic Simulation Algorithms for Chemical Reactions*, Virginia Polytechnic Institute and State University, Blacksburg, Virginia.

## The code

...TODO...

## Chemical master equation

The chemical master equation (CME) is given by

```math
\frac{\partial p(\mathbf{x},t)}{\partial t} = \sum_{n=1}^{N_r} a_n \left(\mathbf{x} - \mathbf{\nu}_n\right) p\left(\mathbf{x} - \mathbf{\nu}_n, t\right) - \sum_{n=1}^{N_r} a_n (\mathbf{x}) p(\mathbf{x}, t),
```

where
$\mathbf{x}\in\mathbb{N}^{N_s}$
is the total number of molecules for each chemical species,
$N_r$
is the number of reaction channels and
$N_s$
is the number of chemical species. Furthermore,
$\nu_n$
is the difference of population counts involved in a single reaction (i.e., the stoichiometric vector),
$a_n (\mathbf{x})$
are the propensity functions (the transition rates for a particular reaction
$n$
).

Sunkara (2019) found an equivalent formulation of the CME, where the probability of certain species counts
$\mathbf{x}\in\mathbb{N}^{N_s}$
is substituted with the probability of certain reaction counts
$\mathbf{r}\in\mathbb{N}^{N_r}$
from a certain initial state
$\mathbf{x}_0$
:

```math
\frac{\partial p(\mathbf{r},t)}{\partial t} = \sum_{n=1}^{N_r} a_n \left(\mathbf{r} - \mathbf{1}_n\right) p\left(\mathbf{r} - \mathbf{1}_n, t\right) - \sum_{n=1}^{N_r} a_n (\mathbf{r}) p(\mathbf{r}, t),
```

where
$\mathbf{1}_n$
is the Kronecker delta. This formulation may be more convenient since it leads to much simpler dynamics: in particular, since the reaction counts can only increase, this gives a simple feedforward propagation of probabilities for the associated Markov chain. We only need a map that from the reaction counts gives the species counts, since these latter variables are needed to compute the propensity functions:

```math
\mathbf{x} = \mathbf{x}_0 + \sum_{n=1}^{N_r} \mathbf{\nu}_n r_n.
```

Note that this map is not in general bijective. Given a solution to the reaction count CME
$p(\mathbf{r},t)$
, we obtain the probabilities associated to each state at time
$t$
through

```math
p(\mathbf{x},t) = \sum_{\mathbf{r}\in\Gamma_\mathbf{x}} p(\mathbf{r},t),
```

where
$`\Gamma_\mathbf{x} := \left\{\mathbf{r} \middle| \mathbf{x} = \mathbf{x}_0 + \Sigma_{n=1}^{N_r} \mathbf{\nu}_n r_n \right\}`$
is the set of reaction counts that gives the same population counts
$\mathbf{x}$
. One problem with this formulation is that the reaction count itself is unbounded: for practical purposes, one can choose to limit the total number of reactions, as long as this choice does not affect the accuracy of the result (this choice will depend on the required length of the simulation).

## Stochastic simulation algorithm

A way to solve the CME is through a Monte-Carlo method, i.e., performing several stochastic realizations, and then calculating the relevant statistics. Gillespie (1977) found such algorithm, which is based on the reaction probability density function:

```math
\rho(\tau,j;\mathbf{x},t) = a_j (\mathbf{x}) \exp \left(-\tau a_0 (\mathbf{x})\right),
```

where
$`a_0 (\mathbf{x}) = \Sigma_{n=1}^{N_r} a_n (\mathbf{x})`$
, while
$\tau$
is the time between a reaction and the next one and
$j$
is the reaction channel index. When calculating the marginal probability density functions of the two variables,
$\tau$
is shown to be an exponentially distributed random variable, while
$j$
is a statistically independent integer random variable with point probabilities
$`a_j (\mathbf{x}) / a_0 (\mathbf{x})`$
. Thus, the Monte-Carlo method is to draw two uniform (pseudo)random numbers
$r_1, r_2 \in U(0,1)$
and then compute:

```math
\tau = - \frac{1}{a_0 (\mathbf{x})} \ln r_1,
```

while
$j$
is given by the smallest integer satisfying

```math
\sum_{n=1}^j a_n (\mathbf{x}) > r_2 a_0 (\mathbf{x})
```

and, finally, carry out the reaction by replacing

```math
t \leftarrow t + \tau, \qquad \mathbf{x} \leftarrow \mathbf{x} + \mathbf{\nu}_j.
```

## Enzyme kinetics

...TODO...
