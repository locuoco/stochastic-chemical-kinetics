// Stochastic enzyme kinetics: Gillespie algorithm

#ifndef SEK_GILLESPIE

#include <random> // uniform_real_distribution, mt19937_64
#include <valarray>
#include <concepts> // floating_point
#include <functional> // function
#include <stdexcept> // domain_error
#include <string>
#include <cmath> // log

template <std::size_t N_s, std::size_t N_r, std::floating_point T = double>
requires (N_r > 0 && N_s > 0)
// T: underlying floating point type
// N_s: number of substances (chemical species)
// N_r: number of reaction channels
class gillespie
{
	std::mt19937_64 gen; // random number generator (64-bit Mersenne Twister)
	std::uniform_real_distribution<T> u_dist = std::uniform_real_distribution<T>(0, 1);

protected:

	std::array<std::valarray<int>, N_r> nu; // stoichiometric vector

public:

	std::valarray<int> x = std::valarray<int>(N_s); // population numbers
	T t = 0; // time

	T total_propensity() const
	// return the total propensity, i.e., the sum of all propensity functions
	// calculated at the current population numbers x.
	{
		T a_tot = 0;

		for (std::size_t i = 0; i < N_r; ++i)
			a_tot += a(i);

		return a_tot;
	}

	void step()
	// a single step of the stochastic simulation algorithm (Gillespie)
	{
		T r1 = u_dist(gen);
		T r2 = u_dist(gen);

		T a_tot = total_propensity();

		if (a_tot == 0)
			return; // no reaction is possible

		T tau = -std::log(r1)/a_tot;

		std::size_t j;

		for (j = 0; j < N_r-1; ++j)
			if (a(j) > r2*a_tot)
				break;

		t += tau;
		x += nu[j];
	}

	void simulate(std::size_t n)
	// simulate for n steps
	{
		for (std::size_t i = 0; i < n; ++i)
			step();
	}

	virtual T a(std::size_t i) const = 0;
	// propensity functions
	//	i: reaction channel index
};

enum enzyme_kinetics_substance : std::size_t {eks_C, eks_P, eks_N};
// enzyme kinetics substances

enum enzyme_kinetics_reaction_channel : std::size_t {ekrc_f, ekrc_b, ekrc_cat, ekrc_N};
// enzyme kinetics reaction channels

template <std::floating_point T = double>
class ekinetics_gillespie : public gillespie<eks_N, ekrc_N, T>
{
	using base = gillespie<eks_N, ekrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	std::array<T, ekrc_N> kappa;
	int ET, ST;

	ekinetics_gillespie(const std::array<T, ekrc_N>& kappa, int ET, int ST)
		: kappa(kappa), ET(ET), ST(ST)
	// default constructor
	{
		x = {0, 0};
		nu[ekrc_f] = {1, 0};
		nu[ekrc_b] = {-1, 0};
		nu[ekrc_cat] = {-1, 1};
	}

	T a(std::size_t i) const override
	{
		if (x[eks_C] > ET || x[eks_C] + x[eks_P] > ST)
			throw std::domain_error(std::string("Current state ")
				+ std::to_string(x[0]) + ", " + std::to_string(x[1]) + " is incopatible with constants of motion.");
		switch (i)
		{
			case ekrc_f:
				return kappa[ekrc_f] * ((ET - x[eks_C]) * (ST - x[eks_C] - x[eks_P]));
			case ekrc_b:
				return kappa[ekrc_b] * x[eks_C];
			case ekrc_cat:
				return kappa[ekrc_cat] * x[eks_C];
			default:
				throw std::out_of_range("Reaction channel index out of bounds");
		}
	}
};

#endif // SEK_GILLESPIE
















