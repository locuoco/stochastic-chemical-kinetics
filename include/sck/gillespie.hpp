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
		// stoichiometric vectors
		//               C, P
		nu[ssrc_f]   = { 1, 0};
		nu[ssrc_b]   = {-1, 0};
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
		// stoichiometric vectors
		//       P
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
		// stoichiometric vectors
		//       P
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

template <std::floating_point T = double>
class gks_gillespie : public gillespie<gks_N, gkrc_N, T>
// Gillespie algorithm applied to Goldbeter-Koshland switch
{
	using base = gillespie<gks_N, gkrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	std::array<T, gkrc_N> kappa;
	long long ET, DT, ST;

	gks_gillespie(const std::array<T, gkrc_N>& kappa, long long ET, long long DT, long long ST) noexcept
		: kappa(kappa), ET(ET), DT(DT), ST(ST)
	// constructor
	//	kappa: the set of the six rate constants associated to the six reactions
	//	ET: total kinase concentration constant (it is conserved)
	//	DT: total phosphatase concentration constant (it is conserved)
	//	ST: total substrate concentration constant (it is conserved)
	{
		// stoichiometric vectors
		//             SP,  C, CP
		nu[gkrc_fe] = { 0,  1,  0};
		nu[gkrc_be] = { 0, -1,  0};
		nu[gkrc_e]  = { 1, -1,  0};
		nu[gkrc_fd] = {-1,  0,  1};
		nu[gkrc_bd] = { 1,  0, -1};
		nu[gkrc_d]  = { 0,  0, -1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		if (x[gks_C] > ET || x[gks_CP] > DT || x[gks_SP] + x[gks_C] + x[gks_CP] > ST)
			throw std::domain_error("Current state "
				+ std::to_string(x[gks_SP]) + ", " + std::to_string(x[gks_C]) + ", " + std::to_string(x[gks_CP])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case gkrc_fe:
				return kappa[gkrc_fe] * ((ET - x[gks_C]) * (ST - x[gks_SP] - x[gks_C] - x[gks_CP]));
			case gkrc_be:
				return kappa[gkrc_be] * x[gks_C];
			case gkrc_e:
				return kappa[gkrc_e] * x[gks_C];
			case gkrc_fd:
				return kappa[gkrc_fd] * ((DT - x[gks_CP]) * x[gks_SP]);
			case gkrc_bd:
				return kappa[gkrc_bd] * x[gks_CP];
			case gkrc_d:
				return kappa[gkrc_d] * x[gks_CP];
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

template <std::floating_point T = double>
class gktq_gillespie : public gillespie<gts_N, gtrc_N, T>
// Gillespie algorithm applied to Goldbeter-Koshland switch tQSSA
{
	using base = gillespie<gts_N, gtrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	T kME, ke, kMD, kd;
	long long ET, DT, ST;

	gktq_gillespie(T kME, T ke, T kMD, T kd, long long ET, long long DT, long long ST) noexcept
		: kME(kME), ke(ke), kMD(kMD), kd(kd), ET(ET), DT(DT), ST(ST)
	// constructor
	//	kME, ke: phosphorylation constants (Michaelis-Menten + catalysis)
	//	kMD, kd: dephosphorylation constants (Michaelis-Menten + catalysis)
	//	ET: total kinase concentration constant (it is conserved)
	//	DT: total phosphatase concentration constant (it is conserved)
	//	ST: total substrate concentration constant (it is conserved)
	{
		// stoichiometric vectors
		//           SP_hat
		nu[gtrc_e]  = { 1};
		nu[gtrc_d]  = {-1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		if (x[gts_SP_hat] > ST)
			throw std::domain_error("Current state "
				+ std::to_string(x[gts_SP_hat])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case gtrc_e:
			{
				long long S_hat = ST - x[gts_SP_hat];
				long long c = 2*ET*S_hat;
				T b = ET + S_hat + kME;
				T Delta = b*b - 2*c;
				return ke*c / (b + sqrt(Delta));
			}
			case gtrc_d:
			{
				long long c = 2*DT*x[gts_SP_hat];
				T b = DT + x[gts_SP_hat] + kMD;
				T Delta = b*b - 2*c;
				return kd*c / (b + sqrt(Delta));
			}
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

template <std::floating_point T = double>
class gksq_gillespie : public gillespie<gss_N, gsrc_N, T>
// Gillespie algorithm applied to Goldbeter-Koshland switch sQSSA
{
	using base = gillespie<gss_N, gsrc_N, T>;

	using base::nu;

public:

	using base::x;
	using base::t;

	T kME, ke, kMD, kd;
	long long ET, DT, ST;

	gksq_gillespie(T kME, T ke, T kMD, T kd, long long ET, long long DT, long long ST) noexcept
		: kME(kME), ke(ke), kMD(kMD), kd(kd), ET(ET), DT(DT), ST(ST)
	// constructor
	//	kME, ke: phosphorylation constants (Michaelis-Menten + catalysis)
	//	kMD, kd: dephosphorylation constants (Michaelis-Menten + catalysis)
	//	ET: total kinase concentration constant (it is conserved)
	//	DT: total phosphatase concentration constant (it is conserved)
	//	ST: total substrate concentration constant (it is conserved)
	{
		// stoichiometric vectors
		//             SP
		nu[gsrc_e]  = { 1};
		nu[gsrc_d]  = {-1};
	}

	T a(std::size_t i) const final override
	// propensity functions
	//	i: reaction channel index
	{
		if (x[gss_SP] > ST)
			throw std::domain_error("Current state "
				+ std::to_string(x[gss_SP])
				+ " is incompatible with constants of motion.");
		switch (i)
		{
			case gsrc_e:
			{
				long long S = ST - x[gss_SP];
				return ke*(ET*S) / (S + kME);
			}
			case gkrc_d:
				return kd*(DT*x[gss_SP]) / (x[gss_SP] + kMD);
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

#endif // SEK_GILLESPIE
















