//  Stochastic enzyme kinetics: Gillespie algorithm
//  Copyright (C) 2023 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SEK_GILLESPIE
#define SEK_GILLESPIE

#include <random> // uniform_real_distribution, mt19937_64
#include <valarray>
#include <concepts> // floating_point
#include <functional> // function
#include <stdexcept> // domain_error, out_of_range
#include <string> // string, to_string
#include <vector>
#include <array>
#include <cmath> // log, sqrt

#include "common.hpp"

template <std::floating_point T = double>
struct gillespie_state
{
	std::valarray<long long> x;
	T t;
};

// template argument deduction guide
gillespie_state() -> gillespie_state<>;

template <std::size_t N_s, std::size_t N_r, std::floating_point T = double>
requires (N_r > 0 && N_s > 0)
// N_s: number of substances (chemical species)
// N_r: number of reaction channels
// T: underlying floating-point type
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

	virtual ~gillespie() = default;

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
	// set t_final to 0 or negative number for infinity
	{
		using std::log;

		T a_tot = total_propensity();

		if (a_tot == 0)
			return false; // no reaction is possible

		T r1 = u_dist(gen);
		T r2 = u_dist(gen);

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

	void simulate(std::vector<gillespie_state<T>>& states, std::size_t n, T t_final = 0, bool b_initial = true)
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

template <std::floating_point T = double>
class ssek_gillespie : public gillespie<sss_N, ssrc_N, T>
// Gillespie algorithm applied to single-substrate enzyme kinetics
{
	using base = gillespie<sss_N, ssrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	std::array<T, ssrc_N> kappa;
	long long ET, ST;

	ssek_gillespie(const std::array<T, ssrc_N>& kappa, long long ET, long long ST) noexcept
		: kappa(kappa), ET(ET), ST(ST)
	// constructor
	//	kappa: the set of the three rate constants associated to the three reactions (f, b, cat)
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[ssrc_f] = {1, 0};
		nu[ssrc_b] = {-1, 0};
		nu[ssrc_cat] = {-1, 1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		if (x[sss_C] > ET || x[sss_C] + x[sss_P] > ST)
			throw std::domain_error("Current state "
				+ std::to_string(x[sss_C]) + ", " + std::to_string(x[sss_P])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case ssrc_f:
				return kappa[ssrc_f] * ((ET - x[sss_C]) * (ST - x[sss_C] - x[sss_P]));
			case ssrc_b:
				return kappa[ssrc_b] * x[sss_C];
			case ssrc_cat:
				return kappa[ssrc_cat] * x[sss_C];
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

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

	tqssa_gillespie(T kM, T kcat, long long ET, long long ST) noexcept
		: kcat(kcat), kM(kM), ET(ET), ST(ST)
	// constructor
	//	kM: Michaelis-Menten constant ( (kb+kcat) / kf )
	//	kcat: catalysis rate constant
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
			throw std::domain_error("Current state "
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
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

template <std::floating_point T = double>
class sqssa_gillespie : public gillespie<1, 1, T>
// Gillespie algorithm applied to sQSSA (standard quasi-steady state approximation)
// untested
{
	using base = gillespie<1, 1, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	T kcat, kM;
	long long ET, ST;

	sqssa_gillespie(T kM, T kcat, long long ET, long long ST) noexcept
		: kcat(kcat), kM(kM), ET(ET), ST(ST)
	// constructor
	//	kM: Michaelis-Menten constant ( (kb+kcat) / kf )
	//	kcat: catalysis rate constant
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

		if (x[sqs_P] > ST)
			throw std::domain_error("Current state "
				+ std::to_string(x[sqs_P])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case 0:
			{
				long long S = ST - x[sqs_P];
				return kcat*(ET*S) / (S + kM);
			}
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

#endif // SEK_GILLESPIE
















