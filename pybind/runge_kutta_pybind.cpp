//  Stochastic enzyme kinetics: Runge-Kutta methods python binding
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
g++ -shared -static -std=c++20 -Wall -Wextra -pedantic -O3 -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\ProgramData\mambaforge\pkgs\pybind11-2.11.1-py311h005e61a_0\Lib\site-packages\pybind11\include -IC:\ProgramData\mambaforge\pkgs\python-3.11.5-h2628c8c_0_cpython\include -LC:\ProgramData\mambaforge\pkgs\python-3.11.5-h2628c8c_0_cpython\libs runge_kutta_pybind.cpp -o runge_kutta.pyd -lPython311

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/sck/runge_kutta.hpp"

using namespace runge_kutta;
namespace py = pybind11;

PYBIND11_MODULE(runge_kutta, m)
{
	py::class_<integrator<>>(m, "integrator");

	py::class_<euler<>, integrator<>>(m, "euler")
		.def(py::init());

	py::class_<midpoint<>, integrator<>>(m, "midpoint")
		.def(py::init());

	py::class_<heun2<>, integrator<>>(m, "heun2")
		.def(py::init());

	py::class_<ralston2<>, integrator<>>(m, "ralston2")
		.def(py::init());

	py::class_<rk4<>, integrator<>>(m, "rk4")
		.def(py::init());

	py::class_<rk4_3_8<>, integrator<>>(m, "rk4_3_8")
		.def(py::init());

	py::class_<ralston4<>, integrator<>>(m, "ralston4")
		.def(py::init());

	py::class_<butcher6<>, integrator<>>(m, "butcher6")
		.def(py::init());

	py::class_<verner8<>, integrator<>>(m, "verner8")
		.def(py::init());
}











