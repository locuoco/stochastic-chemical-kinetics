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
g++ -shared -static-libgcc -static-libstdc++ -std=c++20 -Wall -Wextra -pedantic -O3 -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\Users\aless\Desktop\myLib\include -IC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\include -LC:\ProgramData\Anaconda3\pkgs\python-3.9.12-h6244533_0\libs runge_kutta_pybind.cpp -o runge_kutta.pyd -lPython39

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/sck/runge_kutta.hpp"

using namespace runge_kutta;
namespace py = pybind11;

PYBIND11_MODULE(runge_kutta, m)
{
	py::class_<euler<>>(m, "euler")
		.def(py::init());

	py::class_<midpoint<>>(m, "midpoint")
		.def(py::init());

	py::class_<heun2<>>(m, "heun2")
		.def(py::init());

	py::class_<ralston2<>>(m, "ralston2")
		.def(py::init());

	py::class_<rk4<>>(m, "rk4")
		.def(py::init());

	py::class_<rk4_3_8<>>(m, "rk4_3_8")
		.def(py::init());

	py::class_<ralston4<>>(m, "ralston4")
		.def(py::init());

	py::class_<butcher6<>>(m, "butcher6")
		.def(py::init());

	py::class_<verner8<>>(m, "verner8")
		.def(py::init());
}











