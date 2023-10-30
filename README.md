# Stochastic chemical kinetics
Stochastic chemical kinetics library using Gillespie algorithm (also stochastic simulation algorithm or SSA) and chemical master equation (CME), with application to enzyme kinetics. What follows are the graphs obtained using the library. See the [experiments](https://github.com/locuoco/stochastic-chemical-kinetics/tree/main/experiments) folder for the source code of the Python scripts.

cme_averages.py | cme_completion_times.py| cme_example.py | cme_stationary_dist.py
----------------|------------------|----------------|----------------
![figure_cme_averages](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/118d850f-9a7e-4c26-87a2-b8f531226197) | ![figure_cme_completion_times](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/311aef21-25b4-4b61-9d63-996d5929472d) | <a href="https://user-images.githubusercontent.com/80710565/269441710-972ac818-80db-4719-86e4-fcbb7708edc0.mp4"><img src="https://user-images.githubusercontent.com/80710565/269446124-b4fe5e4d-b0df-48fe-bbaa-8b1b51e41347.png"/></a> | ![figure_cme_stationary_dist](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/8c7ba109-9d05-4835-9517-21018bbccf92)

gillespie_averages.py | gillespie_completion_times.py | gillespie_example.py | gillespie_stationary_dist.py
----------------|------------------|----------------|----------------
![figure_gillespie_averages](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/a342e6de-ca53-4da5-a6ab-c233f5379926) | ![figure_gillespie_completion_times](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/0b8b8b7c-84e0-4f99-ad96-00bfd2aedc73) | ![figure_gillespie_example](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/1cf16610-f668-4a63-8ced-1450efbf46aa) | ![figure_gillespie_stationary_dist](https://github.com/locuoco/stochastic-chemical-kinetics/assets/80710565/700d8e35-e7bf-4cfa-bc0b-f0b650457a41)

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

### Code structure

The code is written as a library in C++20 with Python bindings using [Pybind11](https://github.com/pybind/pybind11). The repository is structured in the following way:

* `experiments`: directory containing example code using the library and experiments for University projects.
* `include/sck`: directory containing the C++ header files to include in implementation files.
  * `cme.hpp`: it includes classes for generic CME equation integration and applications to enzyme kinetics.
  * `gillespie.hpp`: it includes classes for generic Gillespie algorithm and applications to enzyme kinetics.
  * `runge_kutta.hpp`: explicit Runge-Kutta methods used for integration of the CME equation.
  * `tensor.hpp`: classes, aliases and data structures for vectors, matrices and tensors with some helper functions.
* `pybind`: directory containing C++ implementation files that binds the code inside the `include` directory. It also contains the Windows dynamic-link libraries which can be directly imported in Python scripts (on Linux, you will need to recompile them).
  * `cme_pybind.cpp`: Python bindings for CME.
  * `gillespie_pybind.cpp`: Python bindings for Gillespie algorithm.
  * `runge_kutta_pybind.cpp`: Python bindings for Runge-Kutta methods.
* `tests`: Tests that verify the correctness of the current implementation.

### Dependencies

The C++ header files have no dependencies other than a C++20-complying compiler. To compile the Python bindings, [Pybind11](https://github.com/pybind/pybind11) must be installed on the current machine. For more information on how to compile these bindings, see the individual files inside the `pybind` folder. Rename the library paths accordingly, if needed.

The resulting Python library will require NumPy 1.7.0 or any later version.

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
). The propensity functions can be seen as the elements of the transition-rate matrix
$A_{ij}$
for a general master equation:

```math
\frac{d\mathbf{P}(t)}{dt} = A \mathbf{P}(t).
```

Sunkara (2019) found an equivalent formulation of the CME, where the probability of certain species counts
$\mathbf{x}\in\mathbb{N}^{N_s}$
is substituted with the probability of certain reaction counts
$\mathbf{r}\in\mathbb{N}^{N_r}$
from a certain initial state
$\mathbf{x}_0:$

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
$p(\mathbf{r},t),$
we obtain the probabilities associated to each state at time
$t$
through

```math
p(\mathbf{x},t) = \sum_{\mathbf{r}\in\Gamma_\mathbf{x}} p(\mathbf{r},t),
```

where
$`\Gamma_\mathbf{x} := \left\{\mathbf{r} \middle| \mathbf{x} = \mathbf{x}_0 + \Sigma_{n=1}^{N_r} \mathbf{\nu}_n r_n \right\}`$
is the set of reaction counts that gives the same population counts
$\mathbf{x}.$
One problem with this formulation is that the reaction count itself is unbounded: for practical purposes, one can choose to limit the total number of reactions, as long as this choice does not affect the accuracy of the result (this choice will depend on the required length of the simulation).

Both these formulations can be casted into a system of ordinary differential equations with many different variables depending on the maximum allowed population numbers or reaction numbers. For small enough systems, these equations can be simply solved using a standard integration method, for example a Runge-Kutta method. However, we are quickly caught by the curse of dimensionality for systems with many chemical species: a lot of computer memory may be required to perform a simulation.

## Stochastic simulation algorithm

An alternative way to solve the CME is through a Monte-Carlo method, i.e., performing several stochastic realizations, and then calculating the relevant statistics. Gillespie (1977) found such algorithm, which is based on the reaction probability density function:

```math
\rho(\tau,j;\mathbf{x},t) = a_j (\mathbf{x}) \exp \left(-\tau a_0 (\mathbf{x})\right),
```

where
$`a_0 (\mathbf{x}) = \Sigma_{n=1}^{N_r} a_n (\mathbf{x}),`$
while
$\tau$
is the time between a reaction and the next one and
$j$
is the reaction channel index. When calculating the marginal probability density functions of the two variables,
$\tau$
is shown to be an exponentially distributed random variable, while
$j$
is a statistically independent integer random variable with point probabilities
$`a_j (\mathbf{x}) / a_0 (\mathbf{x}).`$
Thus, the Monte-Carlo method is to draw two uniform (pseudo)random numbers
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

Note that, since we require a large amount of realizations to obtain statistically significant results, stochastic simulation algorithm can be slower than direct integration of the CME for systems with few chemical species or small populations. In practice, this algorithm is useful when solving the CME directly would be prohibitively expensive in terms of required memory or computational time.

## Single-substrate enzyme-catalyzed reaction

A single-substrate enzyme-catalyzed reaction can be described by the following chain of reactions

```math
E + S \overset{k_f}{\underset{k_b}\rightleftarrows} C \xrightarrow{k_\textrm{cat}} E + P
```

which correspond, due to the mass action law, to the following system of ordinary differential equations (ODEs):

```math
\frac{d[S]}{dt} = -k_f [E] [S] + k_b [C],
```
```math
\frac{d[E]}{dt} = -k_f [E] [S] + k_b [C] + k_\textrm{cat} [C],
```
```math
\frac{d[C]}{dt} =  k_f [E] [S] - k_b [C] - k_\textrm{cat} [C],
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
$`\frac{d}{dt}[E_T]=0.`$
From the conservation of the total substrate and product concentration, we have
$`[S_T] := [S] + [C] + [P]`$
such that
$`\frac{d}{dt}[S_T]=0.`$
The existence of these constants implies that there are only two independent ODEs to solve, corresponding to two independent variables, for example
$[C]$
and
$[P]:$

```math
\frac{d[C]}{dt} =  k_f ([E_T] - [C]) ([S_T] - [C] - [P]) - k_b [C] - k_\textrm{cat} [C],
```
```math
\frac{d[P]}{dt} = k_\textrm{cat} [C].
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
$[S]:$

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
$`[E] \ll [S] + [P] + K_M,`$
which also implies that
$[C]$
is negligible with respect to the total substrate-product concentration:
$[C] \ll [S] + [P].$
This enables us to make another approximation:
$`[S_T] \approx [S] + [P],`$
so we do not need to know
$[C]$
to compute the time evolution of concentration of substrates and products.

### Total quasi-steady state approximation
Introduced by Borghans et al. (1996), this approximation is analogous to QSSA, with the difference that we perform a change of variable before applying the condition
$d[C]/dt \approx 0$, i.e.:

```math
[\hat{S}] = [S] + [C].
```

Following the same steps as before, we obtain the following expression for
$[C]:$

```math
[C] = \frac{1}{2} \left([E_T] + [\hat{S}] + K_M - \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}\right).
```

From the validity condition given by Tzafriri (2003), one can show that this approximation is reasonably accurate even in the worst case, so the range of validity is much wider for tQSSA than standard QSSA. Substituting in the products rate equation, we obtain

```math
\frac{d[P]}{dt} = \frac{k_\textrm{cat}}{2} \left([E_T] + [\hat{S}] + K_M - \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}\right),
```

which can also be rewritten so that it is more computationally stable for small values of
$[\hat{S}]:$

```math
\frac{d[P]}{dt} = \frac{2 k_\textrm{cat} [E_T] [\hat{S}]}{[E_T] + [\hat{S}] + K_M + \sqrt{([E_T] + [\hat{S}] + K_M)^2 - 4 [E_T] [\hat{S}]}}.
```

Note that
$`[S_T] = [\hat{S}] + [P]`$
holds exactly.

### Stochastic enzyme kinetics
The propensity functions for a single-substrate enzyme-catalyzed reaction are associated to each reaction channel:

```math
E + S \xrightarrow{\kappa_f} C, \qquad C \xrightarrow{\kappa_b} E + S, \qquad C \xrightarrow{\kappa_\textrm{cat}} E + P.
```

So, the three propensity functions are given by

```math
a_f (E, S) = \kappa_f E S, \qquad a_b (C) = \kappa_b C, \qquad a_\textrm{cat} (C) = \kappa_\textrm{cat} C,
```

where
$E,S,C,P$
refer to the numbers of molecules for the different chemical species rather than their concentrations. The relation between a population number
$X$
and the relative molar concentration
$[X]$
is

```math
[X] = \frac{X}{N_A \Omega},
```

where
$`N_A`$
is the Avogadro constant and
$\Omega$
is the total volume of the system. The constants
$`\kappa_i`$
are related to the constants
$`k_i`$
used in the ODE treatment of the reaction rates (Cao and Samuels, 2009):

```math
\kappa_f = \frac{k_f}{N_A \Omega}, \qquad \kappa_b = k_b, \qquad \kappa_\textrm{cat} = k_\textrm{cat}.
```

Remembering that we have only two independent variables thanks to the conservation of
$`S_T = S + C + P`$
and
$`E_T = E + C,`$
the CME gives

```math
\frac{\partial p(C,P,t)}{\partial t} = a_f (C-1, P) p(C-1, P, t) + a_b (C+1) p(C+1, P, t) + a_\textrm{cat} (C+1) p(C+1, P-1, t) - a_0 (C, P) p(C, P, t),
```

where
$`a_f (C, P) = \kappa_f (E_T - C) (S_T - C - P)`$
and
$`a_0 (C, P) = a_f (C, P) + a_b (C) + a_\textrm{cat} (C).`$
The Gillespie algorithm for this system is directly obtained from the propensity functions. Note that we can perform the tQSSA in the stochastic formulation as in the deterministic one by using only one propensity function
$a (\hat{S})$
(Kim and Tyson, 2020), i.e.:

```math
a (\hat{S}) = \frac{2 \kappa_\textrm{cat} E_T \hat{S}}{E_T + \hat{S} + \kappa_M + \sqrt{(E_T + \hat{S} + \kappa_M)^2 - 4 E_T \hat{S}}},
```

where
$`\hat{S} = S + C = S_T - P`$
and
$`\kappa_M = \frac{\kappa_b + \kappa_\textrm{cat}}{\kappa_f},`$
related to the usual Michaelis-Menten constant through

```math
\kappa_M = K_M N_A \Omega.
```

This can be viewed as having only one monomolecular reaction described by

```math
\hat{S} \rightarrow P.
```

## Goldbeter-Koshland switch

The Goldbeter-Koshland (GK) switch consists of a substrate-product pair
$`S,S_P`$
that is interconverted by two enzymes
$E,D:$

```math
E + S \overset{k_{fe}}{\underset{k_{be}}\rightleftarrows} C \xrightarrow{k_e} E + S_P
```
```math
D + S_P \overset{k_{fd}}{\underset{k_{bd}}\rightleftarrows} C_P \xrightarrow{k_d} D + S
```

These reactions are commonly used to describe phosphorylation and dephosphorylation processes, where
$S$
is phosphorylated by kinase
$E$
and
$`S_P`$
is dephosphorylated by phosphatase
$D.$
From the law of mass action we obtain these six differential equations:

```math
\frac{d[S]}{dt} = -k_{fe} [E] [S] + k_{be} [C] + k_d [C_P],
```
```math
\frac{d[S_P]}{dt} = -k_{fd} [D] [S_P] + k_{bd} [C_P] + k_e [C],
```
```math
\frac{d[E]}{dt} = -k_{fe} [E] [S] + k_{be} [C] + k_e [C],
```
```math
\frac{d[D]}{dt} = -k_{fd} [D] [S_P] + k_{bd} [C_P] + k_d [C_P],
```
```math
\frac{d[C]}{dt} =  k_{fe} [E] [S] - k_{be} [C] - k_e [C],
```
```math
\frac{d[C_P]}{dt} =  k_{fd} [D] [S_P] - k_{bd} [C_P] - k_d [C_P].
```

Like for the single-substrate case, we have conserved quantities: the total substrate concentration
$`[S_T] := [S] + [S_P] + [C] + [C_P],`$
the total kinase concentration
$`[E_T] := [E] + [C]`$
and the total phosphatase concentration
$`[D_T] := [D] + [C_P].`$
Having these three conserved quantities, only three of the above differential equations are independent. For example, we can choose

```math
\frac{d[S_P]}{dt} = -k_{fd} ([D_T] - [C_P]) [S_P] + k_{bd} [C_P] + k_e [C],
```
```math
\frac{d[C]}{dt} =  k_{fe} ([E_T] - [C]) ([S_T] - [S_P] - [C] - [C_P]) - k_{be} [C] - k_e [C],
```
```math
\frac{d[C_P]}{dt} =  k_{fd} ([D_T] - [C_P]) [S_P] - k_{bd} [C_P] - k_d [C_P].
```

The QSSA can be obtained by assuming

```math
\frac{d[C]}{dt} \approx 0, \qquad \frac{d[C_P]}{dt} \approx 0,
```

from which we obtain

```math
[C] = \frac{[E_T][S]}{K_{ME} + [S]}, \qquad [C_P] = \frac{[D_T][S_P]}{K_{MD} + [S_P]},
```

where

```math
K_{ME} = \frac{k_{be} + k_e}{k_{fe}}, \qquad K_{MD} = \frac{k_{bd} + k_d}{k_{fd}}.
```

Substituting into the equation for
$`d[S_P]/dt:`$

```math
\frac{d[S_P]}{dt} = \frac{k_e [E_T][S]}{K_{ME} + [S]} - \frac{k_d [D_T][S_P]}{K_{MD} + [S_P]}.
```

Note that
$`[S_T] \approx [S] + [S_P].`$
The tQSSA is obtained by making two variable changes:

```math
[\hat{S}] := [S] + [C], \qquad [\hat{S}_P] = [S_P] + [C_P].
```

Then, applying the quasi-steady state condition as before, we obtain the following expression for the phosphorylated substrate concentration:

```math
\frac{d[S_P]}{dt} = \frac{2 k_e [E_T][\hat{S}]}{[E_T] + [\hat{S}] + K_{ME} + \sqrt{\left([E_T] + [\hat{S}] + K_{ME}\right)^2 - 4 [E_T] [\hat{S}]}}
- \frac{2 k_d [D_T][\hat{S}_P]}{[D_T] + [\hat{S}_P] + K_{MD} + \sqrt{\left([D_T] + [\hat{S}_P] + K_{MD}\right)^2 - 4 [D_T] [\hat{S}_P]}}.
```

Note that
$`[S_T] = [\hat{S}] + [\hat{S}_P]`$
holds exactly. For the stochastic formulation of the GK switch, we recognize six propensity functions related to the six reaction channels

```math
E + S \xrightarrow{\kappa_{fe}} C, \qquad C \xrightarrow{\kappa_{be}} E + S, \qquad C \xrightarrow{\kappa_e} E + S_P,
```
```math
D + S_P \xrightarrow{\kappa_{fd}} C_P, \qquad C_P \xrightarrow{\kappa_{bd}} D + S_P, \qquad C_P \xrightarrow{\kappa_d} D + S.
```

Choosing the coordinates
$`S_P, C, C_P`$
as the independent variables, the propensity functions are given by

```math
a_{fe} (S_P, C, C_P) = \kappa_{fe} (E_T - C) (S_T - S_P - C - C_P), \qquad a_{be} (C) = \kappa_{be} C, \qquad a_e (C) = \kappa_e C,
```
```math
a_{fd} (S_P, C_P) = \kappa_{fd} (D_T - C_P) S_P, \qquad a_{bd} (C_P) = \kappa_{bd} C_P, \qquad a_d (C_P) = \kappa_d C_P.
```

For stochastic tQSSA, we have only two reactions corresponding to

```math
\hat{S} \rightleftarrows \hat{S}_P,
```

for which the propensity functions are given by

```math
a_e (\hat{S}) = \frac{2 k_e E_T \hat{S}}{E_T + \hat{S} + K_{ME} + \sqrt{\left(E_T + \hat{S} + K_{ME}\right)^2 - 4 E_T \hat{S}}}, \qquad
a_d (\hat{S}_P) = \frac{2 k_d D_T \hat{S}_P}{D_T + \hat{S}_P + K_{MD} + \sqrt{\left(D_T + \hat{S}_P + K_{MD}\right)^2 - 4 D_T \hat{S}_P}},
```

with only one independent variable since
$`\hat{S},\hat{S}_P`$
are related by
$`\hat{S} = S_T - \hat{S}_P.`$
