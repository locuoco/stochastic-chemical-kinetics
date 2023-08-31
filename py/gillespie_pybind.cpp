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

g++ -shared -static-libgcc -static-libstdc++ -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\Users\aless\Desktop\myLib\include -IC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\include -LC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\libs gillespie_pybind.cpp -o gillespie.pyd -lPython39

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/sck/gillespie.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gillespie, m)
{
	py::class_<ssek_gillespie<double>>(m, "single_substrate")
		.def(py::init<const std::array<double, ssrc_N>&, long long, long long>())
		.def("step", &ssek_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &ssek_gillespie<double>::x)
		.def_readwrite("t", &ssek_gillespie<double>::t);

	py::class_<tqssa_gillespie<double>>(m, "single_substrate_tqssa")
		.def(py::init<double, double, long long, long long>())
		.def("step", &tqssa_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &tqssa_gillespie<double>::x)
		.def_readwrite("t", &tqssa_gillespie<double>::t);

	py::class_<sqssa_gillespie<double>>(m, "single_substrate_sqssa")
		.def(py::init<double, double, long long, long long>())
		.def("step", &sqssa_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &sqssa_gillespie<double>::x)
		.def_readwrite("t", &sqssa_gillespie<double>::t);

	py::class_<gks_gillespie<double>>(m, "goldbeter_koshland")
		.def(py::init<const std::array<double, gkrc_N>&, long long, long long, long long>())
		.def("step", &gks_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &gks_gillespie<double>::x)
		.def_readwrite("t", &gks_gillespie<double>::t);

	py::class_<gktq_gillespie<double>>(m, "goldbeter_koshland_tqssa")
		.def(py::init<double, double, double, double, long long, long long, long long>())
		.def("step", &gktq_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &gktq_gillespie<double>::x)
		.def_readwrite("t", &gktq_gillespie<double>::t);

	py::class_<gksq_gillespie<double>>(m, "goldbeter_koshland_sqssa")
		.def(py::init<double, double, double, double, long long, long long, long long>())
		.def("step", &gksq_gillespie<double>::step, py::arg("t_final") = 0.)
		.def_readwrite("x", &gksq_gillespie<double>::x)
		.def_readwrite("t", &gksq_gillespie<double>::t);
}











