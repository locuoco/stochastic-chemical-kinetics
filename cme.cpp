//  Stochastic enzyme kinetics: Chemical master equation (CME) test
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

/*

Compilation (GCC/MinGW):
g++ cme.cpp -o cme -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include <iostream> // cout
#include <cassert>

#include "runge_kutta.hpp"
#include "cme.hpp"

int main()
{
	double kf = 10, kb = 9, kcat = 1, kM = (kb + kcat) / kf;
	long long ET = 10, ST = 9;
	double t = 1, dt = 1e-4, max_steps = 200'000;

	rk4_3_8 integ;

	// Full model
	ekinetics_cme system({kf, kb, kcat}, ET, ST);
	system.simulate(integ, max_steps, dt, t);

	auto state = system.get_state();
	for (std::size_t i = 0; i < state.p.size(); ++i)
	{
		auto pop = system.get_pop(i);
		std::cout << pop[eks_C] << ", " << pop[eks_P] << ": " << state.p[i] << '\n';
	}
	std::cout << "t = " << state.t << '\n';

	// tQSSA
	tqssa_cme sys2(kM, kcat, ET, ST);
	sys2.simulate(integ, max_steps, dt, t);

	state = sys2.get_state();
	for (std::size_t i = 0; i < state.p.size(); ++i)
	{
		auto pop = sys2.get_pop(i);
		std::cout << pop[tqs_P] << ": " << state.p[i] << '\n';
	}
	std::cout << "t = " << state.t << '\n';

	// sQSSA
	sqssa_cme sys3(kM, kcat, ET, ST);
	sys3.simulate(integ, max_steps, dt, t);

	state = sys3.get_state();
	for (std::size_t i = 0; i < state.p.size(); ++i)
	{
		auto pop = sys3.get_pop(i);
		std::cout << pop[sqs_P] << ": " << state.p[i] << '\n';
	}
	std::cout << "t = " << state.t << '\n';

	return 0;
}