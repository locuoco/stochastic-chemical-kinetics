// Stochastic enzyme kinetics: Gillespie algorithm

/*

Compilation (MinGW):
g++ gillespie.cpp -o gillespie -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include <iostream> // cout, endl

#include "gillespie.hpp"

int main()
{
	ekinetics_gillespie system({1., 1., .1}, 50, 100);

	system.x = {0, 0};
	system.simulate(1000);

	std::cout << system.x[eks_P] << std::endl;
	std::cout << "time = " << system.t << std::endl;

	return 0;
}