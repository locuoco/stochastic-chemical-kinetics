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

enum ssek_substance : std::size_t {sss_C, sss_P, sss_N};
// single-substrate enzyme kinetics substances

enum ssek_reaction_channel : std::size_t {ssrc_f, ssrc_b, ssrc_cat, ssrc_N};
// single-substrate enzyme kinetics reaction channels

enum tqssa_substance : std::size_t {tqs_P, tqs_N};
// tQSSA substances

enum sqssa_substance : std::size_t {sqs_P, sqs_N};
// sQSSA substances

#endif // SEK_COMMON