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

	std::vector<stoch_state<>> states;

	system.x = {0, 0};
	system.simulate(states, 10);

	for (const auto& state : states)
	{
		std::cout << state.x[eks_C] << '\n';
		std::cout << "time = " << state.t << '\n';
	}

	return 0;
}