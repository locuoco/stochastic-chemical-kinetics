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
g++ -shared -static -std=c++20 -Wall -Wextra -pedantic -O3 -fmax-errors=1 -DMS_WIN64 -fPIC -IC:\ProgramData\mambaforge\pkgs\pybind11-2.11.1-py311h005e61a_0\Lib\site-packages\pybind11\include -IC:\ProgramData\mambaforge\pkgs\python-3.11.5-h2628c8c_0_cpython\include -LC:\ProgramData\mambaforge\pkgs\python-3.11.5-h2628c8c_0_cpython\libs gillespie_pybind.cpp -o gillespie.pyd -lPython311

*/

#include <optional>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../include/sck/gillespie.hpp"

using namespace gillespie;
namespace py = pybind11;

template <typename Class>
void class_defs(py::class_<Class>& c)
{
	c.def("step", Class::step, py::arg("t_final") = 0.);
	c.def("simulate",
		[](Class& self, double t_final, std::size_t max_steps, std::size_t n_sampling, bool noreturn) -> std::optional<py::tuple>
		{
			if (noreturn)
			{
				self.simulate(t_final, max_steps);
				return {};
			}
			else
			{
				list_of_states<Class::num_species>* states = new list_of_states<Class::num_species>();

				self.simulate(*states, t_final, max_steps, n_sampling);

				py::capsule free_when_done(
					states,
					[](void* f)
					{
						delete reinterpret_cast<list_of_states<Class::num_species>*>(f);
					}
				);
				return py::make_tuple<py::return_value_policy::take_ownership>(
						py::array_t<long long>(
							{(long long)states->x.size(), (long long)Class::num_species},
							{(long long)(Class::num_species*sizeof(long long)), (long long)sizeof(long long)},
							reinterpret_cast<long long*>(states->x.data()),
							free_when_done
						),
						py::array_t<double>(
							{states->t.size()},
							{sizeof(double)},
							states->t.data(),
							free_when_done
						)
					);
			}
		},
		py::arg("t_final") = 0.,
		py::arg("max_steps") = Class::default_max_steps(),
		py::arg("n_sampling") = 1,
		py::arg("noreturn") = false);
	c.def_property("x",
		[](const Class& self)
		{
			return py::array_t<long long>(
				{(long long)Class::num_species},
				{sizeof(long long)},
				self.x.data()
			);
		},
		[](Class& self, py::array_t<long long, py::array::c_style | py::array::forcecast> py_array)
		{
			std::copy(py_array.data(), py_array.data() + py_array.size(), self.x.begin());
		});
	c.def_readwrite("t", &Class::t);
}

PYBIND11_MODULE(gillespie, m)
{
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











