/* -*- c++ -*- */
/*
 * Copyright 2019 Free Software Foundation, Inc.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include "log_bcjr_base.h"

log_bcjr_base::log_bcjr_base(int I, int S, int O,
		const std::vector<int> &NS,
		const std::vector<int> &OS)
	: d_I(I), d_S(S), d_O(O), d_ordered_OS(S*I)
{
	if (NS.size() != S*I) {
		throw std::runtime_error("Invalid size for NS.");
	}
	d_NS = NS;

	if (OS.size() != S*I) {
		throw std::runtime_error("Invalid size for OS.");
	}
	d_OS = OS;

	generate_PS_PI();

	//Compute ordered_OS
	std::vector<int>::iterator ordered_OS_it = d_ordered_OS.begin();
	for(int s=0 ; s < S ; ++s) {
		for(size_t i=0 ; i<(d_PS[s]).size() ; ++i) {
			*(ordered_OS_it++) = OS[d_PS[s][i]*I + d_PI[s][i]];
		}
	}
}

void
log_bcjr_base::generate_PS_PI()
{
	d_PS.resize(d_S);
	d_PI.resize(d_S);

	for(int i=0 ; i<d_S ; ++i) {
		d_PS[i].reserve(d_I*d_S); // max possible size
		d_PI[i].reserve(d_I*d_S);

		for(int ii=0 ; ii<d_S ; ++ii) {
			for(int jj=0 ; jj<d_I ; ++jj) {
				if(d_NS[ii*d_I+jj] != i) {
					continue;
				}

				d_PS[i].push_back(ii);
				d_PI[i].push_back(jj);
			}
		}
	}
}

void
log_bcjr_base::compute_fw_metrics(const std::vector<float> &G,
		const std::vector<float> &A0, std::vector<float> &A, size_t K)
{
	A.resize(d_S*(K+1), -std::numeric_limits<float>::max());

	float norm_A = -std::numeric_limits<float>::max();
	std::vector<float>::iterator A_prev, A_curr;
	std::vector<int>::iterator PS_it, ordered_OS_it;

	//Integrate initial forward metrics
	std::copy(A0.begin(), A0.end(), A.begin());

	//Initialize iterators
	A_prev = A.begin();
	A_curr = A.begin() + d_S;
	for(std::vector<float>::const_iterator G_k = G.begin() ;
			G_k != G.end() ; G_k += d_O) {

        ordered_OS_it = d_ordered_OS.begin();
		for(int s=0 ; s < d_S ; ++s) {
			//Iterators for previous state and previous input lists
			PS_it=d_PS[s].begin();

			//Loop
			for(size_t i=0 ; i<(d_PS[s]).size() ; ++i) {
				// Equivalent to:
				// *A_curr = _max_star(*A_curr,
				// A_prev[PS[s][i]] + G_k[d_OS[PS[s][i]*I + PI[s][i]]);
				*A_curr = _max_star(*A_curr,
						A_prev[*(PS_it++)] + G_k[*(ordered_OS_it++)]);
			}

			//Update iterators
			++A_curr;
		}

		//Advance A_prev
		A_prev += d_S;

		//Metrics normalization
		norm_A = _max_star(&(*(A_prev)), d_S);
		std::transform(A_prev, A_curr, A_prev,
				std::bind2nd(std::minus<float>(), norm_A));
	}
}

void
log_bcjr_base::compute_bw_metrics(const std::vector<float> &G,
		const std::vector<float> &BK, std::vector<float> &B, size_t K)
{
	B.resize(d_S*(K+1), -std::numeric_limits<float>::max());

	float norm_B = -std::numeric_limits<float>::max();
	std::vector<float>::reverse_iterator B_next, B_curr;
	std::vector<int>::reverse_iterator NS_it, OS_it;

	//Integrate initial forward metrics
	std::copy(BK.rbegin(), BK.rend(), B.rbegin());

	//Initialize iterators
	B_curr = B.rbegin() + d_S;
	B_next = B.rbegin();
	for(std::vector<float>::const_reverse_iterator G_k = G.rbegin() ;
			G_k != G.rend() ; G_k += d_O) {

		//Iterators for next state and next output lists
		NS_it=d_NS.rbegin();
		OS_it=d_OS.rbegin();
		for(int s=0 ; s < d_S ; ++s) {
			//Loop
			for(size_t i=0 ; i < d_I ; ++i) {
				*B_curr = _max_star(*B_curr,
						B_next[(d_S-1)-*(NS_it++)] + G_k[(d_O-1)-*(OS_it++)]);
			}

			//Update iterators
			++B_curr;
		}

		//Advance B_next (go back, as it is a reverse iterator...)
		B_next += d_S;

		//Metrics normalization
		norm_B = _max_star(&(*B_curr)+1, d_S);
		std::transform(B_next, B_curr, B_next,
				std::bind2nd(std::minus<float>(), norm_B));
	}
}

void
log_bcjr_base::compute_app(const std::vector<float> &A, const std::vector<float> &B,
		const std::vector<float> &G, size_t K, std::vector<float> &out)
{
	std::vector<float>::const_iterator A_it = A.begin();
	std::vector<float>::const_iterator B_it = B.begin() + d_S;

	out.reserve(d_S*d_I*K);

	for(std::vector<float>::const_iterator G_k = G.begin() ;
			G_k != G.end() ; G_k += d_O) {

		for(int s=0 ; s < d_S ; ++s) {
			for (int i=0 ; i < d_I ; ++i) {
				out.push_back(B_it[d_NS[s*d_I+i]] + G_k[d_OS[s*d_I+i]] + *A_it);
			}

			//Update forward iterator
			++A_it;
		}

		//Update backward iterator
		B_it += d_S;
	}
}

void
log_bcjr_base::log_bcjr_algorithm(const std::vector<float> &A0,
		const std::vector<float> &BK, const std::vector<float> &in,
		std::vector<float> &out)
{
	std::vector<float> A, B;
	size_t K = in.size()/d_O;

	//Forward recursion
	compute_fw_metrics(in, A0, A, K);

	//Backward recursion
	compute_bw_metrics(in, BK, B, K);

	//Compute branch APP
	compute_app(A, B, in, K, out);
}

