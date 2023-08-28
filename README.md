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
\frac{\partial p(\mathbf{x},t)}{\partial t} = \sum_{n=1}^{N_r} a_n (\mathbf{x} - \mathbf{\nu}_n) p(\mathbf{x} - \mathbf{\nu}_n, t) - \sum_{n=1}^{N_r} a_n (\mathbf{x}) p(\mathbf{x}, t),
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
\frac{\partial p(\mathbf{r},t)}{\partial t} = \sum_{n=1}^{N_r} a_n (\mathbf{r} - \mathbf{1}_n) p(\mathbf{r} - \mathbf{1}_n, t) - \sum_{n=1}^{N_r} a_n (\mathbf{r}) p(\mathbf{r}, t),
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

Both these formulations can be casted into a system of ordinary differential equations with many different variables depending on the maximum allowed population numbers or reaction numbers. For small enough systems, these equations can be simply solved using a standard integration method, for example a Runge-Kutta method. However, we are quickly caught by the curse of dimensionality for systems with many chemical species: a lot of computer memory may be required to perform a simulation.

## Stochastic simulation algorithm

An alternative way to solve the CME is through a Monte-Carlo method, i.e., performing several stochastic realizations, and then calculating the relevant statistics. Gillespie (1977) found such algorithm, which is based on the reaction probability density function:

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
and then compute

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

Note that, since we require a large amount of realizations to obtain statistically significant results, stochastic simulation algorithm can be slower than direct integration of the CME for systems with few chemical species or small populations. In practice, this algorithm is useful when solving the CME direcly would be prohibitively expensive in terms of computational time.

## Single-substrate enzyme-catalyzed reaction

A single-substrate enzyme-catalyzed reaction can be described by the following chain of reactions

```math
E + S \rightleftarrows C \rightarrow P + S
```

which correspond, due to the mass action law, to the following system of ordinary differential equations (ODEs):

```math
\frac{d[S]}{dt} = -k_f [S] [E] + k_b [C],
```
```math
\frac{d[E]}{dt} = -k_f [S] [E] + k_b [C] + k_\textrm{cat} [C],
```
```math
\frac{d[C]}{dt} =  k_f [S] [E] - k_b [C] - k_\textrm{cat} [C],
```
```math
\frac{d[P]}{dt} = k_\textrm{cat} [C],
```

where
$[S]$
is the substrate molar concentration,
$[E]$
is the enzyme molar concentration,
$[C]$
is the enzyme-substrate complex molar concentration and
$[P]$
is the product molar concentration. From the conservation of the total enzyme concentration, we have
$`[E_T] := [E] + [C]`$
such that
$`\frac{d}{dt}[E_T]=0`$
. From the conservation of the total substrate and product concentration, we have
$`[S_T] := [S] + [C] + [P]`$
such that
$`\frac{d}{dt}[S_T]=0`$
. The existence of these constants means that there are only two independent ODEs to solve, corresponding to two independent variables, for example
$[C]$
and
$[P]$
:

```math
\frac{d[C]}{dt} =  k_f ([S_T] - [C] - [P]) ([E_T] - [C]) - k_b [C] - k_\textrm{cat} [C],
```
```math
\frac{d[P]}{dt} = k_\textrm{cat} [C],
```

The number of independent variables can be further reduced to one using approximations.

### Quasi-steady state approximation
The quasi-steady state approximation (QSSA, or sQSSA, meaning *standard QSSA*) is obtained under the assumption that
$[C]$
does not change appreciably before
$[S]$
varies appreciably. Mathematically speaking, we require that

```math
\frac{d[C]}{dt} \approx 0,
```

from which we obtain a closed expression for the enzyme-substrate complex concentration
$[C]$
as a function of the substrate concentration
$[S]$
:

```math
[C] = \frac{[E_T] [S]}{[S] + K_M},
```

where
$`K_M = (k_b + k_\textrm{cat})/k_f`$
is the Michaelis-Menten constant. Substituting in the products rate equation, we obtain

```math
\frac{d[P]}{dt} = \frac{k_\textrm{cat} [E_T] [S]}{[S] + K_M},
```

which is the Michaelis-Menten rate law. Segel and Slemrod (1989) showed that this approximation is valid for
$`[E] \ll [S] + [P] + K_M`$
, which also implies that
$[C]$
is negligible with respect to the total substrate-product concentration:
$[C] \ll [S] + [P]$
. This enables us to make another approximation:
$`[S_T] \approx [S] + [P]`$
, so we do not need to know
$[C]$
to compute the time evolution of concentration of substrates and products.

### Total quasi-steady state approximation
Introduced by Borghans et al. (1996), this approximation is analogous to QSSA, with the difference that we perform a change of variable before applying the condition
$d[C]/dt \approx 0$, i.e.:

```math
[\hat{S}] = [S] + [C].
```

Following the same steps as before, we obtain the following expression for
$[C]$
:

```math
[C] = \frac{1}{2} \left([E_T] + [\hat{S}] + K_M - \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}\right).
```

From the validity condition given by Tzafriri (2003), one can show that this approximation is reasonably accurate even in the worst case, so the range of validity is much wider for tQSSA than standard QSSA. Substituting in the products rate equation, we obtain

```math
\frac{d[P]}{dt} = \frac{k_\textrm{cat}}{2} \left([E_T] + [\hat{S}] + K_M - \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}\right),
```

which can also be rewritten so that it is more computationally stable for small values of
$[\hat{S}]$
:

```math
\frac{d[P]}{dt} = \frac{2 k_\textrm{cat} [E_T] [\hat{S}]}{[E_T] + [\hat{S}] + K_M + \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}}.
```

Note that
$`[S_T] = [\hat{S}] + [P]`$
holds exactly.

### Stochastic enzyme kinetics
...TODO...
