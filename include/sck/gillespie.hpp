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

namespace gillespie
{
	template <std::floating_point T = double>
	struct state
	{
		std::valarray<long long> x;
		T t;
	};

	// template argument deduction guide
	state() -> state<>;

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

		void simulate(std::vector<state<T>>& states, std::size_t n, T t_final = 0)
		// simulate for n steps or until t >= t_final, and save the states inside a list (initial and final state are included)
		// set t_final to 0 or negative number for infinity
		// if the total propensity gets to zero, the simulation will be terminated
		{
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
	class single_substrate : public gillespie<2, 3, T>
	// Gillespie algorithm applied to single-substrate enzyme kinetics
	{
		using base = gillespie<2, 3, T>;

		using base::nu;

	public:

		enum species : std::size_t {C, P, num_species};
		enum reaction_channels : std::size_t {f, b, cat, num_rc};

		using base::x;
		using base::t;

		std::array<T, num_rc> kappa;
		long long ET, ST;

		single_substrate(T kf, T kb, T kcat, long long ET, long long ST) noexcept
			: kappa{kf, kb, kcat}, ET(ET), ST(ST)
		// constructor
		//	kf, kb, kcat: the set of the three rate constants associated to the three reactions
		//	ET: total enzyme concentration constant (it is conserved)
		//	ST: total substrate/product concentration constant (it is conserved)
		{
			// stoichiometric vectors
			//          C, P
			nu[f]   = { 1, 0};
			nu[b]   = {-1, 0};
			nu[cat] = {-1, 1};
		}

		T a(std::size_t i) const final override
		// propensity functions
		//	i: reaction channel index
		{
			if (x[C] > ET || x[C] + x[P] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[C]) + ", " + std::to_string(x[P])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case f:
					return kappa[f] * ((ET - x[C]) * (ST - x[C] - x[P]));
				case b:
					return kappa[b] * x[C];
				case cat:
					return kappa[cat] * x[C];
				default:
					throw std::out_of_range("Reaction channel index out of bounds.");
			}
		}
	};

	template <std::floating_point T = double>
	class single_substrate_tqssa : public gillespie<1, 1, T>
	// Gillespie algorithm applied to tQSSA (total quasi-steady state approximation)
	{
		using base = gillespie<1, 1, T>;

		using base::nu;

	public:

		enum species : std::size_t {P, num_species};

		using base::x;
		using base::t;

		T kcat, kM;
		long long ET, ST;

		single_substrate_tqssa(T kM, T kcat, long long ET, long long ST) noexcept
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

			if (x[P] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[P])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case 0:
				{
					long long S_hat = ST - x[P];
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
	class single_substrate_sqssa : public gillespie<1, 1, T>
	// Gillespie algorithm applied to sQSSA (standard quasi-steady state approximation)
	// untested
	{
		using base = gillespie<1, 1, T>;

		using base::nu;

	public:

		enum species : std::size_t {P, num_species};

		using base::x;
		using base::t;

		T kcat, kM;
		long long ET, ST;

		single_substrate_sqssa(T kM, T kcat, long long ET, long long ST) noexcept
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

			if (x[P] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[P])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case 0:
				{
					long long S = ST - x[P];
					return kcat*(ET*S) / (S + kM);
				}
				default:
					throw std::out_of_range("Reaction channel index out of bounds.");
			}
		}
	};

	template <std::floating_point T = double>
	class goldbeter_koshland : public gillespie<3, 6, T>
	// Gillespie algorithm applied to Goldbeter-Koshland switch
	{
		using base = gillespie<3, 6, T>;

		using base::nu;

	public:

		enum species : std::size_t {SP, C, CP, num_species};
		enum reaction_channels : std::size_t {fe, be, e, fd, bd, d, num_rc};

		using base::x;
		using base::t;

		std::array<T, num_rc> kappa;
		long long ET, DT, ST;

		goldbeter_koshland(T kfe, T kbe, T ke, T kfd, T kbd, T kd, long long ET, long long DT, long long ST) noexcept
			: kappa{kfe, kbe, ke, kfd, kbd, kd}, ET(ET), DT(DT), ST(ST)
		// constructor
		//	kfe, kbe, ke, kfd, kbd, kd: the set of the six rate constants associated to the six reactions
		//	ET: total kinase concentration constant (it is conserved)
		//	DT: total phosphatase concentration constant (it is conserved)
		//	ST: total substrate concentration constant (it is conserved)
		{
			// stoichiometric vectors
			//        SP,  C, CP
			nu[fe] = { 0,  1,  0};
			nu[be] = { 0, -1,  0};
			nu[e]  = { 1, -1,  0};
			nu[fd] = {-1,  0,  1};
			nu[bd] = { 1,  0, -1};
			nu[d]  = { 0,  0, -1};
		}

		T a(std::size_t i) const final override
		// propensity functions
		//	i: reaction channel index
		{
			if (x[C] > ET || x[CP] > DT || x[SP] + x[C] + x[CP] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[SP]) + ", " + std::to_string(x[C]) + ", " + std::to_string(x[CP])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case fe:
					return kappa[fe] * ((ET - x[C]) * (ST - x[SP] - x[C] - x[CP]));
				case be:
					return kappa[be] * x[C];
				case e:
					return kappa[e] * x[C];
				case fd:
					return kappa[fd] * ((DT - x[CP]) * x[SP]);
				case bd:
					return kappa[bd] * x[CP];
				case d:
					return kappa[d] * x[CP];
				default:
					throw std::out_of_range("Reaction channel index out of bounds.");
			}
		}
	};

	template <std::floating_point T = double>
	class goldbeter_koshland_tqssa : public gillespie<1, 2, T>
	// Gillespie algorithm applied to Goldbeter-Koshland switch tQSSA
	{
		using base = gillespie<1, 2, T>;

		using base::nu;

	public:

		enum species : std::size_t {SP_hat, num_species};
		enum reaction_channels : std::size_t {e, d, num_rc};

		using base::x;
		using base::t;

		T kME, ke, kMD, kd;
		long long ET, DT, ST;

		goldbeter_koshland_tqssa(T kME, T ke, T kMD, T kd, long long ET, long long DT, long long ST) noexcept
			: kME(kME), ke(ke), kMD(kMD), kd(kd), ET(ET), DT(DT), ST(ST)
		// constructor
		//	kME, ke: phosphorylation constants (Michaelis-Menten + catalysis)
		//	kMD, kd: dephosphorylation constants (Michaelis-Menten + catalysis)
		//	ET: total kinase concentration constant (it is conserved)
		//	DT: total phosphatase concentration constant (it is conserved)
		//	ST: total substrate concentration constant (it is conserved)
		{
			// stoichiometric vectors
			//      SP_hat
			nu[e]  = { 1};
			nu[d]  = {-1};
		}

		T a(std::size_t i) const final override
		// propensity functions
		//	i: reaction channel index
		{
			if (x[SP_hat] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[SP_hat])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case e:
				{
					long long S_hat = ST - x[SP_hat];
					long long c = 2*ET*S_hat;
					T b = ET + S_hat + kME;
					T Delta = b*b - 2*c;
					return ke*c / (b + sqrt(Delta));
				}
				case d:
				{
					long long c = 2*DT*x[SP_hat];
					T b = DT + x[SP_hat] + kMD;
					T Delta = b*b - 2*c;
					return kd*c / (b + sqrt(Delta));
				}
				default:
					throw std::out_of_range("Reaction channel index out of bounds.");
			}
		}
	};

	template <std::floating_point T = double>
	class goldbeter_koshland_sqssa : public gillespie<1, 2, T>
	// Gillespie algorithm applied to Goldbeter-Koshland switch sQSSA
	{
		using base = gillespie<1, 2, T>;

		using base::nu;

	public:

		enum species : std::size_t {SP, num_species};
		enum reaction_channels : std::size_t {e, d, num_rc};

		using base::x;
		using base::t;

		T kME, ke, kMD, kd;
		long long ET, DT, ST;

		goldbeter_koshland_sqssa(T kME, T ke, T kMD, T kd, long long ET, long long DT, long long ST) noexcept
			: kME(kME), ke(ke), kMD(kMD), kd(kd), ET(ET), DT(DT), ST(ST)
		// constructor
		//	kME, ke: phosphorylation constants (Michaelis-Menten + catalysis)
		//	kMD, kd: dephosphorylation constants (Michaelis-Menten + catalysis)
		//	ET: total kinase concentration constant (it is conserved)
		//	DT: total phosphatase concentration constant (it is conserved)
		//	ST: total substrate concentration constant (it is conserved)
		{
			// stoichiometric vectors
			//        SP
			nu[e]  = { 1};
			nu[d]  = {-1};
		}

		T a(std::size_t i) const final override
		// propensity functions
		//	i: reaction channel index
		{
			if (x[SP] > ST)
				throw std::domain_error("Current state "
					+ std::to_string(x[SP])
					+ " is incompatible with constants of motion.");
			switch (i)
			{
				case e:
				{
					long long S = ST - x[SP];
					return ke*(ET*S) / (S + kME);
				}
				case d:
					return kd*(DT*x[SP]) / (x[SP] + kMD);
				default:
					throw std::out_of_range("Reaction channel index out of bounds.");
			}
		}
	};
} // namespace gillespie

#endif // SEK_GILLESPIE
















