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
g++ -shared -static-libgcc -static-libstdc++ -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\Users\aless\Desktop\myLib\include -IC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\include -LC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\libs gillespie_pybind.cpp -o gillespie.pyd -lPython39

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/sck/gillespie.hpp"

using namespace gillespie;
namespace py = pybind11;

template <typename Gillespie>
std::vector<state<>> simulate(Gillespie& self, std::size_t n, double t_final, bool b_initial)
{
	std::vector<state<>> states;
	self.simulate(states, n, t_final, b_initial);
	return states;
}

PYBIND11_MODULE(gillespie, m)
{
	py::class_<state<>>(m, "state")
		.def(py::init())
		.def_readwrite("x", &state<>::x)
		.def_readwrite("t", &state<>::t);

	py::class_<single_substrate<>>(m, "single_substrate")
		.def(py::init<double, double, double, long long, long long>(),
			py::arg("kf"), py::arg("kb"), py::arg("kcat"), py::arg("ET"), py::arg("ST"))
		.def("step", &single_substrate<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<single_substrate<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &single_substrate<>::x)
		.def_readwrite("t", &single_substrate<>::t);

	py::class_<single_substrate_tqssa<>>(m, "single_substrate_tqssa")
		.def(py::init<double, double, long long, long long>(), py::arg("kM"), py::arg("kcat"), py::arg("ET"), py::arg("ST"))
		.def("step", &single_substrate_tqssa<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<single_substrate_tqssa<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &single_substrate_tqssa<>::x)
		.def_readwrite("t", &single_substrate_tqssa<>::t);

	py::class_<single_substrate_sqssa<>>(m, "single_substrate_sqssa")
		.def(py::init<double, double, long long, long long>(), py::arg("kM"), py::arg("kcat"), py::arg("ET"), py::arg("ST"))
		.def("step", &single_substrate_sqssa<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<single_substrate_sqssa<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &single_substrate_sqssa<>::x)
		.def_readwrite("t", &single_substrate_sqssa<>::t);

	py::class_<goldbeter_koshland<>>(m, "goldbeter_koshland")
		.def(py::init<double, double, double, double, double, double, long long, long long, long long>(),
			py::arg("kfe"), py::arg("kbe"), py::arg("ke"), py::arg("kfd"), py::arg("kbd"), py::arg("kd"),
			py::arg("ET"), py::arg("DT"), py::arg("ST"))
		.def("step", &goldbeter_koshland<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<goldbeter_koshland<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &goldbeter_koshland<>::x)
		.def_readwrite("t", &goldbeter_koshland<>::t);

	py::class_<goldbeter_koshland_tqssa<>>(m, "goldbeter_koshland_tqssa")
		.def(py::init<double, double, double, double, long long, long long, long long>(),
			py::arg("kME"), py::arg("ke"), py::arg("kMD"), py::arg("kd"), py::arg("ET"), py::arg("DT"), py::arg("ST"))
		.def("step", &goldbeter_koshland_tqssa<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<goldbeter_koshland_tqssa<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &goldbeter_koshland_tqssa<>::x)
		.def_readwrite("t", &goldbeter_koshland_tqssa<>::t);

	py::class_<goldbeter_koshland_sqssa<>>(m, "goldbeter_koshland_sqssa")
		.def(py::init<double, double, double, double, long long, long long, long long>(),
			py::arg("kME"), py::arg("ke"), py::arg("kMD"), py::arg("kd"), py::arg("ET"), py::arg("DT"), py::arg("ST"))
		.def("step", &goldbeter_koshland_sqssa<>::step, py::arg("t_final") = 0.)
		.def("simulate",
			&simulate<goldbeter_koshland_sqssa<>>,
			py::arg("n_steps"),
			py::arg("t_final") = 0.,
			py::arg("b_initial") = true)
		.def_readwrite("x", &goldbeter_koshland_sqssa<>::x)
		.def_readwrite("t", &goldbeter_koshland_sqssa<>::t);
}











