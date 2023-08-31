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
g++ tests/cme.cpp -o cme -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include <iostream> // cout
#include <cassert>
#include <cmath> // fabs

#include "../include/sck/runge_kutta.hpp"
#include "../include/sck/cme.hpp"

void test_cme_tqssa()
// Test at high enzyme concentration that tQSSA agrees
// with the exact formulation using chemical master equation.
// The test passes if the average products population at a certain
// time agree within 1% relative error
{
	using std::fabs;

	double kf = 10, kb = 9, kcat = 1, kM = (kb + kcat) / kf;
	long long ET = 10, ST = 9;
	double t = 1, dt = 1e-4;

	rk4_3_8 integ; // 4-th order Runge-Kutta integrator

	// Full model
	cme::single_substrate sys_ss(kf, kb, kcat, ET, ST);
	sys_ss.simulate(integ, dt, t);

	// tQSSA
	cme::single_substrate_tqssa sys_tq(kM, kcat, ET, ST);
	sys_tq.simulate(integ, dt, t);

	// sQSSA
	cme::single_substrate_sqssa sys_sq(kM, kcat, ET, ST);
	sys_sq.simulate(integ, dt, t);

	double ss_mean_prod = sys_ss.mean(sys_ss.P);
	double tq_mean_prod = sys_tq.mean(sys_tq.P);
	double sq_mean_prod = sys_sq.mean(sys_sq.P);

	std::cout << "ek: " << ss_mean_prod << " +/- " << sys_ss.sd(sys_ss.P) << '\n';
	std::cout << "tq: " << tq_mean_prod << " +/- " << sys_tq.sd(sys_tq.P) << '\n';
	std::cout << "sq: " << sq_mean_prod << " +/- " << sys_sq.sd(sys_sq.P) << '\n';

	assert(fabs(tq_mean_prod - ss_mean_prod) / ss_mean_prod < .01);
}

int main()
{
	test_cme_tqssa();

	return 0;
}




























