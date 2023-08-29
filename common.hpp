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

enum gks_substance : std::size_t {gks_SP, gks_C, gks_CP, gks_N};
// Goldbeter-Koshland switch subtances

enum gks_reaction_channel : std::size_t {gkrc_fe, gkrc_be, gkrc_e, gkrc_fd, gkrc_bd, gkrc_d, gkrc_N};
// Goldbeter-Koshland switch reaction channels

enum gktq_substance : std::size_t {gts_SP_hat, gts_N};
// tQSSA Goldbeter-Koshland switch subtances

enum gktq_reaction_channel : std::size_t {gtrc_e, gtrc_d, gtrc_N};
// tQSSA Goldbeter-Koshland switch reaction channels

enum gksq_substance : std::size_t {gss_SP, gss_N};
// sQSSA Goldbeter-Koshland switch subtances

enum gksq_reaction_channel : std::size_t {gsrc_e, gsrc_d, gsrc_N};
// sQSSA Goldbeter-Koshland switch reaction channels

#endif // SEK_COMMON