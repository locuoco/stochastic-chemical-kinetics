//  Stochastic enzyme kinetics: Runge-Kutta methods
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

#ifndef SEK_RUNGE_KUTTA
#define SEK_RUNGE_KUTTA

#include <concepts> // floating_point
#include <functional> // function
#include <valarray>
#include <array>

template <std::floating_point T = double>
struct integrator
{
	virtual ~integrator() = default;

	virtual void step(std::valarray<T>& x, T dt, std::function<const std::valarray<T>& (const std::valarray<T>&)> f) = 0;
};

namespace runge_kutta
{
	template <std::size_t Stages, std::floating_point T = double>
	class runge_kutta : public integrator<T>
	// all explicit Runge-Kutta methods inherit from this class
	// `Stages` is the number of stages of the method (number of force evaluations).
	{
		std::valarray<T> k[Stages], prev_x;

	public:

		// parameters of the Runge-Kutta method
		const std::array<T, Stages*Stages> pars;

		template <typename ... Ts>
		requires (sizeof...(Ts) == Stages*Stages)
		runge_kutta(Ts ... pars) : pars{T(pars)...} {}
		// constructor:
		//	pars... is a variadic argument which contains the parameters for the explicit Runge-Kutta method.

		virtual ~runge_kutta() = default;

		void step(std::valarray<T>& x, T dt, std::function<const std::valarray<T>& (const std::valarray<T>&)> f)
		//	x is the system state
		//	dt is the integration step
		//	f is the function that gives the derivative of x wrt time, i.e. f(x) = dx/dt
		{
			using std::size_t;

			prev_x = x;
			k[0] = f(x);
			for (size_t i = 1; i < Stages; ++i)
			{
				x = 0;
				// sum smaller contributes (in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
					x += pars[(i-1)*Stages + j+1] * k[j] * dt;
				x += prev_x;

				k[i] = f(x);
			}

			x = 0;
			for (size_t i = 0; i < Stages; ++i)
				x += pars[(Stages-1)*Stages + i] * k[i] * dt;
			x += prev_x;
		}
	};

	template <std::floating_point T = double>
	struct euler : runge_kutta<1, T>
	// Euler method (1st order, 1 stage)
	{
		euler() : runge_kutta<1, T>{1}
		{}
	};

	euler() -> euler<>;

	template <std::floating_point T = double>
	struct rk2 : runge_kutta<2, T>
	// Parametrized Runge-Kutta 2 method (2nd order, 2 stages)
	// All second-order explicit Runge-Kutta methods can be written
	// as setting the parameter of this method properly
	{
		rk2(long double a) : runge_kutta<2, T>
			{
				a, a,
				1 - 1/(2*a), 1/(2*a)
			}
		// constructor: `a` is the parameter of the parametrized Runge-Kutta 2 method
		{}
	};

	rk2(long double) -> rk2<>;

	template <std::floating_point T = double>
	struct midpoint : rk2<T>
	// Midpoint method (2nd order, 2 stages)
	{
		midpoint() : rk2<T>(.5L)
		{}
	};

	midpoint() -> midpoint<>;

	template <std::floating_point T = double>
	struct heun2 : rk2<T>
	// Heun method (2nd order, 2 stages)
	{
		heun2() : rk2<T>(1)
		{}
	};

	heun2() -> heun2<>;

	template <std::floating_point T = double>
	struct ralston2 : rk2<T>
	// Ralston method (2nd order, 2 stages)
	{
		ralston2() : rk2<T>(2/3.L)
		{}
	};

	ralston2() -> ralston2<>;

	template <std::floating_point T = double>
	struct rk4 : runge_kutta<4, T>
	// Classical Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4() : runge_kutta<4, T>
			{
				.5L, .5L, 0, 0,
				.5L, 0, .5L, 0,
				1, 0, 0, 1,
				1/6.L, 1/3.L, 1/3.L, 1/6.L
			}
		{}
	};

	rk4() -> rk4<>;

	template <std::floating_point T = double>
	struct rk4_3_8 : runge_kutta<4, T>
	// 3/8-rule Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4_3_8() : runge_kutta<4, T>
			{
				1/3.L, 1/3.L, 0, 0,
				2/3.L, -1/3.L, 1, 0,
				1, 1, -1, 1,
				1/8.L, 3/8.L, 3/8.L, 1/8.L
			}
		{}
	};

	rk4_3_8() -> rk4_3_8<>;

	template <std::floating_point T = double>
	struct ralston4 : runge_kutta<4, T>
	// Ralston method (4th order, 4 stages)
	{
		ralston4() : runge_kutta<4, T>
			{
				.4L, .4L, 0, 0,
				.45573725L, .29697761L, .15875964L, 0,
				1, .2181004L, -3.05096516L, 3.83286476L,
				.17476028L, -.55148066L, 1.2055356L, .17118478L
			}
		{}
	};

	ralston4() -> ralston4<>;

	template <std::floating_point T = double>
	struct butcher6 : runge_kutta<7, T>
	// Butcher method (6th order, 7 stages)
	{
		butcher6() : runge_kutta<7, T>
			{
				1/3.L, 1/3.L, 0, 0, 0, 0, 0,
				2/3.L, 0, 2/3.L, 0, 0, 0, 0,
				1/3.L, 1/12.L, 1/3.L, -1/12.L, 0, 0, 0,
				.5L, -1/16.L, 9/8.L, -3/16.L, -3/8.L, 0, 0,
				.5L, 0, 9/8.L, -3/8.L, -3/4.L, .5L, 0,
				1, 9/44.L, -9/11.L, 63/44.L, 18/11.L, -16/11.L, 0,
				11/120.L, 0, 27/40.L, 27/40.L, -4/15.L, -4/15.L, 11/120.L
			}
		{}
	};

	butcher6() -> butcher6<>;

	template <std::floating_point T = double>
	struct verner8 : runge_kutta<11, T>
	// Verner method (8th order, 11 stages)
	// S. Hippolyte, A. K. Richard, "A New Eighth Order Runge-Kutta Family Method", Journal of Mathematics Research, 2019
	{
		verner8() : runge_kutta<11, T>
			{
				.5L, .5L, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				.5L, .25L, .25L, 0, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, 1/7.L, (-7-3*s21)/98, (21+5*s21)/49, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, (11+s21)/84, 0, (18+4*s21)/63, (21-s21)/252, 0, 0, 0, 0, 0, 0,
				.5L, (5+s21)/48, 0, (9+s21)/36, (-231+14*s21)/360, (63-7*s21)/80, 0, 0, 0, 0, 0,
				(7-s21)/14, (10-s21)/42, 0, (-432+92*s21)/315, (633-145*s21)/90, (-504+115*s21)/70, (63-13*s21)/35, 0, 0, 0, 0,
				(7-s21)/14, 1/14.L, 0, 0, 0, (14-3*s21)/126, (13-3*s21)/63, 1/9.L, 0, 0, 0,
				.5L, 1/32.L, 0, 0, 0, (91-21*s21)/576, 11/72.L, (-385-75*s21)/1152, (63+13*s21)/128, 0, 0,
				(7+s21)/14, 1/14.L, 0, 0, 0, 1/9.L, (-733-147*s21)/2205, (515+111*s21)/504, (-51-11*s21)/56, (132+28*s21)/245, 0,
				1, 0, 0, 0, 0, (-42+7*s21)/18, (-18+28*s21)/45, (-273-53*s21)/72, (301+53*s21)/72, (28-28*s21)/45, (49-7*s21)/18,
				1/20.L, 0, 0, 0, 0, 0, 0, 49/180.L, 16/45.L, 49/180.L, 1/20.L
			}
		{}

		private:

			static constexpr long double s21 = 4.582575694955840006588047193728L; // sqrt(21)
	};

	verner8() -> verner8<>;
} // namespace runge_kutta

#endif // SEK_RUNGE_KUTTA
































