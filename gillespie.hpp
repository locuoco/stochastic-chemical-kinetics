// Stochastic enzyme kinetics: Gillespie algorithm

#ifndef SEK_GILLESPIE

#include <random> // uniform_real_distribution, mt19937_64
#include <valarray>
#include <concepts> // floating_point
#include <functional> // function
#include <stdexcept> // domain_error
#include <string>
#include <vector>
#include <cmath> // log, sqrt

template <std::floating_point T = double>
struct stoch_state
{
	std::valarray<long long> x;
	T t;
};

template <std::size_t N_s, std::size_t N_r, std::floating_point T = double>
requires (N_r > 0 && N_s > 0)
// T: underlying floating point type
// N_s: number of substances (chemical species)
// N_r: number of reaction channels
class gillespie
// Gillespie general algorithm
{
	std::mt19937_64 gen; // random number generator (64-bit Mersenne Twister)
	std::uniform_real_distribution<T> u_dist = std::uniform_real_distribution<T>(0, 1);

protected:

	std::array<std::valarray<long long>, N_r> nu; // stoichiometric vector

public:

	std::valarray<long long> x = std::valarray<long long>(N_s); // population numbers
	T t = 0; // time

	gillespie() noexcept
	// default constructor
	{
		x = 0;
	}

	T total_propensity() const
	// return the total propensity, i.e., the sum of all propensity functions
	// calculated at the current population numbers x.
	{
		T a_tot = 0;

		for (std::size_t i = 0; i < N_r; ++i)
			a_tot += a(i);

		return a_tot;
	}

	bool step(T t_final = 0)
	// a single step of the stochastic simulation algorithm (Gillespie)
	// return whether a reaction has been performed successfully before time t_final or not
	{
		using std::log;

		T r1 = u_dist(gen);
		T r2 = u_dist(gen);

		T a_tot = total_propensity();

		if (a_tot == 0)
			return false; // no reaction is possible

		T tau = -log(r1)/a_tot;

		if (t + tau > t_final && t_final > 0)
			return false; // reaction would be performed after t_final

		std::size_t j;

		T a_accum = 0;
		for (j = 0; j < N_r-1; ++j)
		{
			a_accum += a(j);
			if (a_accum > r2*a_tot)
				break;
		}

		t += tau;
		x += nu[j];

		return true;
	}

	void simulate(std::size_t n, T t_final = 0)
	// simulate for n steps or until t >= t_final
	// set t_final to 0 or negative number for infinity
	// if the total propensity gets to zero, the simulation will be terminated
	{
		for (std::size_t i = 0; i < n && (t <= t_final || t_final <= 0); ++i)
			if (!step(t_final))
				break;
	}

	void simulate(std::vector<stoch_state<T>>& states, std::size_t n, T t_final = 0, bool b_initial = true)
	// simulate for n steps or until t >= t_final, and save the states inside a list (final state is always included)
	// set t_final to 0 or negative number for infinity
	// set b_initial to false if the initial state should not be included (default is true)
	// if the total propensity gets to zero, the simulation will be terminated
	{
		if (b_initial)
			states.push_back({x, t});
		for (std::size_t i = 0; i < n && (t <= t_final || t_final <= 0); ++i)
		{
			if (!step(t_final))
				break;
			states.push_back({x, t});
		}
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
// Gillespie algorithm applied to enzyme kinetics
{
	using base = gillespie<eks_N, ekrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	std::array<T, ekrc_N> kappa;
	long long ET, ST;

	ekinetics_gillespie(const std::array<T, ekrc_N>& kappa, long long ET, long long ST) noexcept
		: kappa(kappa), ET(ET), ST(ST)
	// constructor
	//	kappa: the set of the three rate constants associated to the three reactions
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[ekrc_f] = {1, 0};
		nu[ekrc_b] = {-1, 0};
		nu[ekrc_cat] = {-1, 1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		if (x[eks_C] > ET || x[eks_C] + x[eks_P] > ST)
			throw std::domain_error(std::string("Current state ")
				+ std::to_string(x[eks_C]) + ", " + std::to_string(x[eks_P])
				+ " is incompatible with constants of motion.");
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

enum tqssa_substance : std::size_t {tqs_P, tqs_N};
// tQSSA substances

template <std::floating_point T = double>
class tqssa_gillespie : public gillespie<1, 1, T>
// Gillespie algorithm applied to tQSSA (total quasi-steady state approximation)
{
	using base = gillespie<1, 1, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	T kcat, kM;
	long long ET, ST;

	tqssa_gillespie(T kcat, T kM, long long ET, long long ST) noexcept
		: kcat(kcat), kM(kM), ET(ET), ST(ST)
	// constructor
	//	kcat: catalysis rate constant
	//	kM: Michaelis-Menten constant ( (kb+kcat) / kf )
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[0] = {1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		using std::sqrt;

		if (x[tqs_P] > ST)
			throw std::domain_error(std::string("Current state ")
				+ std::to_string(x[tqs_P])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case 0:
			{
				long long S_hat = ST - x[tqs_P];
				long long c = 2*ET*S_hat;
				T b = ET + S_hat + kM;
				T Delta = b*b - 2*c;
				return kcat*c / (b + sqrt(Delta));
			}
			default:
				throw std::out_of_range("Reaction channel index out of bounds");
		}
	}
};

#endif // SEK_GILLESPIE
















