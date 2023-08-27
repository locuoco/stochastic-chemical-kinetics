//  Stochastic enzyme kinetics: Gillespie algorithm test
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
g++ tests/gillespie.cpp -o gillespie -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include <iostream> // cout
#include <cassert>
#include <cmath> // fabs

#include "../gillespie.hpp"

void test_gillespie_tqssa_prod()
// Test for a specific combination of parameters that tQSSA agrees
// with the exact formulation using Gillespie algorithm.
// The test passes if the average products population at a certain
// time agree within 1% relative error
{
	using std::fabs;

	std::size_t n = 10'000, max_steps = 1'000, t = 2;
	double kf = 10, kb = 9, kcat = 1, kM = (kb + kcat) / kf;
	long long ET = 10, ST = 9;
	double P1 = 0, P2 = 0;

	ekinetics_gillespie sys1({kf, kb, kcat}, ET, ST);
	tqssa_gillespie sys2(kM, kcat, ET, ST);

	for (std::size_t i = 0; i < n; ++i)
	{
		sys1.x = 0;
		sys1.t = 0;
		sys1.simulate(max_steps, t);
		P1 += sys1.x[eks_P];
	}
	P1 /= n*ST;

	for (std::size_t i = 0; i < n; ++i)
	{
		sys2.x = 0;
		sys2.t = 0;
		sys2.simulate(max_steps, t);
		P2 += sys2.x[tqs_P];
	}
	P2 /= n*ST;

	std::cout << P1 << '\n';
	std::cout << P2 << '\n';

	assert(fabs(P1 - P2) / P1 < .01); // not more than 1% relative error
}

void test_gillespie_tqssa_completion()
// Test for a specific combination of parameters that tQSSA agrees
// with the exact formulation using Gillespie algorithm.
// The test passes if completion times agree within 2% relative error
{
	using std::fabs;

	std::size_t n = 10'000, max_steps = 1'000;
	double kf = 10, kb = 9, kcat = 1, kM = (kb + kcat) / kf;
	long long ET = 10, ST = 9;
	double t1 = 0, t2 = 0;

	ekinetics_gillespie sys1({kf, kb, kcat}, ET, ST);
	tqssa_gillespie sys2(kM, kcat, ET, ST);

	for (std::size_t i = 0; i < n; ++i)
	{
		sys1.x = 0;
		sys1.t = 0;
		sys1.simulate(max_steps);
		t1 += sys1.t;
	}
	t1 /= n;

	for (std::size_t i = 0; i < n; ++i)
	{
		sys2.x = 0;
		sys2.t = 0;
		sys2.simulate(max_steps);
		t2 += sys2.t;
	}
	t2 /= n;

	std::cout << t1 << '\n';
	std::cout << t2 << '\n';

	assert(fabs(t1 - t2) / t1 < .02); // not more than 2% relative error
}

int main()
{
	test_gillespie_tqssa_prod();
	test_gillespie_tqssa_completion();

	return 0;
}

















