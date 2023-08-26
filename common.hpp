//  Stochastic enzyme kinetics: common data structures
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

#ifndef SEK_COMMON
#define SEK_COMMON

enum enzyme_kinetics_substance : std::size_t {eks_C, eks_P, eks_N};
// enzyme kinetics substances

enum enzyme_kinetics_reaction_channel : std::size_t {ekrc_f, ekrc_b, ekrc_cat, ekrc_N};
// enzyme kinetics reaction channels

enum tqssa_substance : std::size_t {tqs_P, tqs_N};
// tQSSA substances

#endif // SEK_COMMON