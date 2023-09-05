//  Stochastic enzyme kinetics: Gillespie algorithm python binding
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

Compilation (MinGW):
g++ -shared -static-libgcc -static-libstdc++ -std=c++20 -Wall -Wextra -pedantic -O3 -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\Users\aless\Desktop\myLib\include -IC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\include -LC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\libs gillespie_pybind.cpp -o gillespie.pyd -lPython39

*/

#include <optional>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/sck/gillespie.hpp"

using namespace gillespie;
namespace py = pybind11;

template <typename Class>
void class_defs(py::class_<Class>& c)
{
	c.def("step", Class::step, py::arg("t_final") = 0.);
	c.def("simulate",
		[](Class& self, double t_final, std::size_t max_steps, std::size_t n_sampling, bool noreturn) -> std::optional<std::vector<state<>>>
		{
			if (noreturn)
			{
				self.simulate(t_final, max_steps);
				return {};
			}
			else
			{
				std::vector<state<>> states;
				self.simulate(states, t_final, max_steps, n_sampling);
				return states;
			}
		},
		py::arg("t_final") = 0.,
		py::arg("max_steps") = Class::default_max_steps(),
		py::arg("n_sampling") = 1,
		py::arg("noreturn") = false);
	c.def_readwrite("x", &Class::x);
	c.def_readwrite("t", &Class::t);
}

PYBIND11_MODULE(gillespie, m)
{
	py::class_<state<>>(m, "state")
		.def(py::init())
		.def_readwrite("x", &state<>::x)
		.def_readwrite("t", &state<>::t);

	py::class_<single_substrate<>> c_single_substrate(m, "single_substrate");
	c_single_substrate.def(py::init<double, double, double, long long, long long>(),
		py::arg("kf"), py::arg("kb"), py::arg("kcat"), py::arg("ET"), py::arg("ST"));
	class_defs(c_single_substrate);

	py::enum_<single_substrate<>::species>(c_single_substrate, "species")
		.value("C", single_substrate<>::species::C)
		.value("P", single_substrate<>::species::P);

	py::enum_<single_substrate<>::reaction_channels>(c_single_substrate, "reaction_channels")
		.value("f",   single_substrate<>::reaction_channels::f)
		.value("b",   single_substrate<>::reaction_channels::b)
		.value("cat", single_substrate<>::reaction_channels::cat);

	py::class_<single_substrate_tqssa<>> c_single_substrate_tqssa(m, "single_substrate_tqssa");
	c_single_substrate_tqssa.def(py::init<double, double, long long, long long>(),
		py::arg("kM"), py::arg("kcat"), py::arg("ET"), py::arg("ST"));
	class_defs(c_single_substrate_tqssa);

	py::enum_<single_substrate_tqssa<>::species>(c_single_substrate_tqssa, "species")
		.value("P", single_substrate_tqssa<>::species::P);

	py::enum_<single_substrate_tqssa<>::reaction_channels>(c_single_substrate_tqssa, "reaction_channels")
		.value("f", single_substrate_tqssa<>::reaction_channels::f);

	py::class_<single_substrate_sqssa<>> c_single_substrate_sqssa(m, "single_substrate_sqssa");
	c_single_substrate_sqssa.def(py::init<double, double, long long, long long>(),
		py::arg("kM"), py::arg("kcat"), py::arg("ET"), py::arg("ST"));
	class_defs(c_single_substrate_sqssa);

	py::enum_<single_substrate_sqssa<>::species>(c_single_substrate_sqssa, "species")
		.value("P", single_substrate_sqssa<>::species::P);

	py::enum_<single_substrate_sqssa<>::reaction_channels>(c_single_substrate_sqssa, "reaction_channels")
		.value("f", single_substrate_sqssa<>::reaction_channels::f);

	py::class_<goldbeter_koshland<>> c_goldbeter_koshland(m, "goldbeter_koshland");
	c_goldbeter_koshland.def(py::init<double, double, double, double, double, double, long long, long long, long long>(),
		py::arg("kfe"), py::arg("kbe"), py::arg("ke"), py::arg("kfd"), py::arg("kbd"), py::arg("kd"),
		py::arg("ET"), py::arg("DT"), py::arg("ST"));
	class_defs(c_goldbeter_koshland);

	py::enum_<goldbeter_koshland<>::species>(c_goldbeter_koshland, "species")
		.value("SP", goldbeter_koshland<>::species::SP)
		.value("C",  goldbeter_koshland<>::species::C)
		.value("CP", goldbeter_koshland<>::species::CP);

	py::enum_<goldbeter_koshland<>::reaction_channels>(c_goldbeter_koshland, "reaction_channels")
		.value("fe", goldbeter_koshland<>::reaction_channels::fe)
		.value("be", goldbeter_koshland<>::reaction_channels::be)
		.value("e",  goldbeter_koshland<>::reaction_channels::e)
		.value("fd", goldbeter_koshland<>::reaction_channels::fd)
		.value("bd", goldbeter_koshland<>::reaction_channels::bd)
		.value("d",  goldbeter_koshland<>::reaction_channels::d);

	py::class_<goldbeter_koshland_tqssa<>> c_goldbeter_koshland_tqssa(m, "goldbeter_koshland_tqssa");
	c_goldbeter_koshland_tqssa.def(py::init<double, double, double, double, long long, long long, long long>(),
		py::arg("kME"), py::arg("ke"), py::arg("kMD"), py::arg("kd"), py::arg("ET"), py::arg("DT"), py::arg("ST"));
	class_defs(c_goldbeter_koshland_tqssa);

	py::enum_<goldbeter_koshland_tqssa<>::species>(c_goldbeter_koshland_tqssa, "species")
		.value("SP_hat", goldbeter_koshland_tqssa<>::species::SP_hat);

	py::enum_<goldbeter_koshland_tqssa<>::reaction_channels>(c_goldbeter_koshland_tqssa, "reaction_channels")
		.value("e", goldbeter_koshland_tqssa<>::reaction_channels::e)
		.value("d", goldbeter_koshland_tqssa<>::reaction_channels::d);

	py::class_<goldbeter_koshland_sqssa<>> c_goldbeter_koshland_sqssa(m, "goldbeter_koshland_sqssa");
	c_goldbeter_koshland_sqssa.def(py::init<double, double, double, double, long long, long long, long long>(),
		py::arg("kME"), py::arg("ke"), py::arg("kMD"), py::arg("kd"), py::arg("ET"), py::arg("DT"), py::arg("ST"));
	class_defs(c_goldbeter_koshland_sqssa);

	py::enum_<goldbeter_koshland_sqssa<>::species>(c_goldbeter_koshland_sqssa, "species")
		.value("SP", goldbeter_koshland_sqssa<>::species::SP);

	py::enum_<goldbeter_koshland_sqssa<>::reaction_channels>(c_goldbeter_koshland_sqssa, "reaction_channels")
		.value("e", goldbeter_koshland_sqssa<>::reaction_channels::e)
		.value("d", goldbeter_koshland_sqssa<>::reaction_channels::d);
}











