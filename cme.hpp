//  Stochastic enzyme kinetics: chemical master equation (CME)
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

#ifndef SEK_CME
#define SEK_CME

#include <valarray>
#include <concepts> // floating_point
#include <stdexcept> // out_of_range
#include <array>
#include <algorithm> // max, min
#include <vector>
#include <cmath> // sqrt
#include <string>

#include "common.hpp"
#include "tensor.hpp"

template <std::floating_point T = double>
struct cme_state
{
	std::valarray<T> p;
	T t;
};

// template argument deduction guide
cme_state() -> cme_state<>;

template <std::size_t N_s, std::size_t N_r, std::floating_point T = double>
requires (N_r > 0 && N_s > 0)
// N_s: number of substances (chemical species)
// N_r: number of reaction channels
// T: underlying floating-point type
class cme
// General chemical master equation (CME) solver
{
	static constexpr std::size_t moments_max_order = 3;

	physics::vec<long long, N_s> n_max;
	mutable physics::vec<long long, N_s> x;
	mutable std::array<std::array<T, moments_max_order+1>, N_s> m;
	std::valarray<T> p, dp;
	T t;
	mutable bool calculated_stats = false;

	long long n_elems() const noexcept
	{
		long long n = 1;
		for (std::size_t i = 0; i < N_s; ++i)
			n *= n_max[i];
		return n;
	}

	bool out_of_bounds(const physics::vec<long long, N_s>& y) const noexcept
	{
		for (std::size_t i = 0; i < N_s; ++i)
			if (y[i] < 0 || y[i] >= n_max[i])
				return true;
		return false;
	}

	void calc_moments() const noexcept
	// compute moments up to prefixed order
	{
		if (calculated_stats)
			return;
		using std::size_t;
		size_t n = n_elems();
		for (size_t j = 0; j < N_s; ++j)
		{
			m[j][0] = 1;
			for (size_t i = 1; i <= moments_max_order; ++i)
				m[j][i] = 0;
		}
		x = 0;
		for (size_t pop_i = 0; pop_i < n; ++pop_i)
		{
			for (size_t j = 0; j < N_s; ++j)
			{
				long long xn = 1;
				for (size_t i = 1; i <= moments_max_order; ++i)
				{
					xn *= x[j];
					m[j][i] += p[pop_i] * xn;
				}
			}
			++x[N_s-1];
			for (size_t j = N_s; j --> 1; )
				if (x[j] == n_max[j])
				{
					x[j] = 0;
					++x[j-1];
				}
		}
		calculated_stats = true;
	}

protected:

	std::array<physics::vec<long long, N_s>, N_r> nu; // stoichiometric vector

public:

	cme(const physics::vec<long long, N_s>& n_max)
		: n_max(n_max)
	// constructor
	//	n_max: max population numbers (we must consider a finite number of possible states to make the problem computable)
	//	dt: integration time step
	{
		for (std::size_t i = 0; i < N_s; ++i)
			if (n_max[i] <= 0)
				throw std::out_of_range("The maximum population numbers must be greater than 0.");
		p.resize(n_elems(), 0);
		dp.resize(n_elems());
		p[0] = 1;
	}

	virtual ~cme() = default;

	std::size_t get_index(const physics::vec<long long, N_s>& y) const noexcept
	// get index from population numbers y
	{
		long long index = 0;
		for (std::size_t i = 0; i < N_s; ++i)
			index = index*n_max[i] + y[i];
		return index;
	}

	physics::vec<long long, N_s> get_pop(std::size_t index) const noexcept
	// get population numbers from index
	{
		physics::vec<long long, N_s> y;
		for (std::size_t i = N_s; i --> 0; )
		{
			y[i] = index % n_max[i];
			index /= n_max[i];
		}
		return y;
	}

	template <typename Integ>
	void step(Integ& integ, T dt)
	// a single step integrating the chemical master equation
	//	integ: integrator
	//	dt: integration time step
	{
		using std::size_t;

		calculated_stats = false;

		auto f = [this](const std::valarray<T>& p) -> const std::valarray<T>&
		// this lambda function calculates the time derivative of the probabilities
		// associated to each possible state, as dictated by the CME.
		{
			size_t n = n_elems();
			x = 0;
			for (size_t pop_i = 0; pop_i < n; ++pop_i)
			{
				dp[pop_i] = 0;
				for (size_t r_i = 0; r_i < N_r; ++r_i)
				{
					if (!out_of_bounds(x - nu[r_i]))
						dp[pop_i] += a(x - nu[r_i], r_i) * p[get_index(x - nu[r_i])];
					dp[pop_i] -= a(x, r_i) * p[pop_i];
				}
				++x[N_s-1];
				for (size_t j = N_s; j --> 1; )
					if (x[j] == n_max[j])
					{
						x[j] = 0;
						++x[j-1];
					}
			}
			return dp;
		};
		integ.step(p, dt, f);
		p += dp * dt;
		t += dt;
	}

	template <typename Integ>
	[[maybe_unused]] std::size_t simulate(Integ& integ, T dt, T t_final)
	// simulate until t >= t_final
	// return the number of steps
	{
		std::size_t i;
		for (i = 0; t <= t_final; ++i)
			step(integ, dt);
		return i;
	}

	template <typename Integ>
	[[maybe_unused]] std::size_t simulate(Integ& integ, std::vector<cme_state<T>>& states, T dt, T t_final, bool b_initial = true)
	// simulate until t >= t_final, and save the states inside a list (final state is always included)
	// set b_initial to false if the initial state should not be included (default is true)
	// return the number of steps
	{
		if (b_initial)
			states.push_back({p, t});
		std::size_t i;
		for (i = 0; t <= t_final; ++i)
		{
			step(integ, dt);
			states.push_back({p, t});
		}
		return i;
	}

	cme_state<T> get_state() const noexcept
	// return the current state
	{
		return {p, t};
	}

	T mean(std::size_t s_i) const
	// return the mean of the s_i-th variable (substance)
	{
		if (s_i >= N_s)
			throw std::out_of_range("Unknown substance with index " + std::to_string(s_i) + ".");
		calc_moments();
		return m[s_i][1];
	}

	T msq(std::size_t s_i) const
	// return the mean square of the s_i-th variable (substance)
	{
		if (s_i >= N_s)
			throw std::out_of_range("Unknown substance with index " + std::to_string(s_i) + ".");
		calc_moments();
		return m[s_i][2];
	}

	T sd(std::size_t s_i) const
	{
		using std::sqrt;
		if (s_i >= N_s)
			throw std::out_of_range("Unknown substance with index " + std::to_string(s_i) + ".");
		calc_moments();
		T arg = m[s_i][2] - m[s_i][1]*m[s_i][1];
		return arg > 0 ? sqrt(arg) : 0;
	}

	T nth_moment(std::size_t s_i, std::size_t n) const
	// return the n-th moment of the s_i-th variable (substance)
	{
		using std::pow;
		if (s_i >= N_s)
			throw std::out_of_range("Unknown substance with index " + std::to_string(s_i) + ".");
		if (n <= moments_max_order)
		{
			calc_moments();
			return m[s_i][n];
		}
		size_t n_ = n_elems();
		T mom = 0;
		x = 0;
		for (size_t pop_i = 0; pop_i < n_; ++pop_i)
		{
			mom += p[pop_i] * pow(x[s_i], n);
			++x[N_s-1];
			for (size_t j = N_s; j --> 1; )
				if (x[j] == n_max[j])
				{
					x[j] = 0;
					++x[j-1];
				}
		}
		return mom;
	}

	virtual T a(const physics::vec<long long, N_s>& y, std::size_t r_i) const = 0;
	// propensity functions
	//	y: population numbers
	//	r_i: reaction channel index
};

template <std::floating_point T = double>
class ekinetics_cme : public cme<eks_N, ekrc_N, T>
// Chemical master equation applied to enzyme kinetics
{
	using base = cme<eks_N, ekrc_N, T>;

	using base::nu;

public:

	std::array<T, ekrc_N> kappa;
	long long ET, ST;

	ekinetics_cme(const std::array<T, ekrc_N>& kappa, long long ET, long long ST) noexcept
		: base({ET+1, ST+1}), kappa(kappa), ET(ET), ST(ST)
	// constructor
	//	kappa: the set of the three rate constants associated to the three reactions (f, b, cat)
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[ekrc_f] = {1, 0};
		nu[ekrc_b] = {-1, 0};
		nu[ekrc_cat] = {-1, 1};
	}

	T a(const physics::vec<long long, eks_N>& y, std::size_t i) const final override
	// propensity functions
	//	y: population numbers
	//	i: reaction channel index
	{
		if (y[eks_C] + y[eks_P] > ST)
			return 0;
		switch (i)
		{
			case ekrc_f:
				return kappa[ekrc_f] * ((ET - y[eks_C]) * (ST - y[eks_C] - y[eks_P]));
			case ekrc_b:
				return kappa[ekrc_b] * y[eks_C];
			case ekrc_cat:
				return kappa[ekrc_cat] * y[eks_C];
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

template <std::floating_point T = double>
class tqssa_cme : public cme<1, 1, T>
// Chemical master equation applied to tQSSA (total quasi-steady state approximation)
{
	using base = cme<1, 1, T>;

	using base::nu;

public:

	T kcat, kM;
	long long ET, ST;

	tqssa_cme(T kM, T kcat, long long ET, long long ST) noexcept
		: base({ST+1}), kcat(kcat), kM(kM), ET(ET), ST(ST)
	// constructor
	//	kM: Michaelis-Menten constant ( (kb+kcat) / kf )
	//	kcat: catalysis rate constant
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[0] = {1};
	}

	T a(const physics::vec<long long, tqs_N>& y, std::size_t i) const final override
	// propensity functions
	//	y: population numbers
	//	i: reaction channel index
	{
		using std::sqrt;

		switch (i)
		{
			case 0:
			{
				long long S_hat = ST - y[tqs_P];
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
class sqssa_cme : public cme<1, 1, T>
// Chemical master equation applied to sQSSA (standard quasi-steady state approximation)
{
	using base = cme<1, 1, T>;

	using base::nu;

public:

	T kcat, kM;
	long long ET, ST;

	sqssa_cme(T kM, T kcat, long long ET, long long ST) noexcept
		: base({ST+1}), kcat(kcat), kM(kM), ET(ET), ST(ST)
	// constructor
	//	kM: Michaelis-Menten constant ( (kb+kcat) / kf )
	//	kcat: catalysis rate constant
	//	ET: total enzyme concentration constant (it is conserved)
	//	ST: total substrate/product concentration constant (it is conserved)
	{
		nu[0] = {1};
	}

	T a(const physics::vec<long long, sqs_N>& y, std::size_t i) const final override
	// propensity functions
	//	y: population numbers
	//	i: reaction channel index
	{
		using std::sqrt;

		switch (i)
		{
			case 0:
			{
				long long S = ST - y[sqs_P];
				return kcat*(ET*S) / (S + kM);
			}
			default:
				throw std::out_of_range("Reaction channel index out of bounds.");
		}
	}
};

#endif // SEK_CME


























