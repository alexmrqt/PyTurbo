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

#include "viterbi.h"

viterbi::viterbi(int I, int S, int O,
		const std::vector<int> &NS,
		const std::vector<int> &OS)
	: d_I(I), d_S(S), d_O(O)
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
}

void
viterbi::generate_PS_PI()
{
	d_PS.resize(d_S);
	d_PI.resize(d_S);

	for(int i=0;i<d_S;i++) {
		d_PS[i].resize(d_I*d_S); // max possible size
		d_PI[i].resize(d_I*d_S);
		int j=0;
		for(int ii=0;ii<d_S;ii++) for(int jj=0;jj<d_I;jj++) {
			if(d_NS[ii*d_I+jj]!=i) continue;
			d_PS[i][j]=ii;
			d_PI[i][j]=jj;
			j++;
		}
		d_PS[i].resize(j);
		d_PI[i].resize(j);
	}
}

void
viterbi::viterbi_algorithm(int K, int S0, int SK, const float *in,
		unsigned char *out)
{
	viterbi_algorithm(d_I, d_S, d_O, d_NS, d_OS, d_PS, d_PI, K, S0, SK, in, out);
}

void
viterbi::viterbi_algorithm(int I, int S, int O, const std::vector<int> &NS,
	const std::vector<int> &OS, const std::vector< std::vector<int> > &PS,
	const std::vector< std::vector<int> > &PI, int K, int S0, int SK,
	const float *in, unsigned char *out)
{
	int tb_state, pidx;
	float can_metric = std::numeric_limits<float>::max();
	float min_metric = std::numeric_limits<float>::max();

	std::vector<int> trace(K*S);
	std::vector<float> alpha_prev(S, std::numeric_limits<float>::max());
	std::vector<float> alpha_curr(S, std::numeric_limits<float>::max());

	std::vector<float>::iterator alpha_curr_it;
	std::vector<int>::const_iterator PS_it, PI_it;
	std::vector<int>::iterator trace_it = trace.begin();

	//If initial state was specified
	if(S0 != -1) {
		alpha_prev[S0] = 0.0;
	}

	for(float* in_k=(float*)in ; in_k < (float*)in + K*O ; in_k += O) {
		//Current path metric iterator
		alpha_curr_it = alpha_curr.begin();
		for(int s=0 ; s < S ; ++s) {
			//Iterators for previous state and previous input lists
			PS_it=PS[s].begin();
			PI_it=PI[s].begin();

			//ACS for state s
			//Pre-loop
			*alpha_curr_it = alpha_prev[*PS_it] + in_k[OS[(*PS_it)*I + (*PI_it)]];
			*trace_it = 0;

			//Loop
			for(size_t i=1 ; i<(PS[s]).size() ; ++i) {
				//Update PS/PI iterators
				++PS_it;
				++PI_it;

				//ADD
				can_metric = alpha_prev[*PS_it] + in_k[OS[(*PS_it)*I + (*PI_it)]];

				//COMPARE
				if(can_metric < *alpha_curr_it) {
					//SELECT
					*alpha_curr_it = can_metric;
					*trace_it = i;  //Store previous input index for traceback
				}
			}

			//Update trace and path metric iterator
			++trace_it;
			++alpha_curr_it;
		}

		//Metrics normalization
		min_metric = *std::min_element(alpha_curr.begin(), alpha_curr.end());
		std::transform(alpha_curr.begin(), alpha_curr.end(), alpha_curr.begin(),
			std::bind2nd(std::minus<double>(), min_metric));

		//At this point, current path metrics becomes previous path metrics
		alpha_prev.swap(alpha_curr);
	}

	//If final state was specified
	if(SK != -1) {
		tb_state = SK;
	}
	else{
		//at this point, alpha_prev contains the path metrics of states after time K
		tb_state = (int)(min_element(alpha_prev.begin(), alpha_prev.end()) - alpha_prev.begin());
	}

	//Traceback
	trace_it = trace.end() - S; //place trace_it at the last time index

	for(unsigned char* out_k = out+K-1 ; out_k >= out ; --out_k) {
		//Retrieve previous input index from trace
		pidx=*(trace_it + tb_state);
		//Update trace_it for next output symbol
		trace_it -= S;

		//Output previous input
		*out_k = (unsigned char) PI[tb_state][pidx];

		//Update tb_state with the previous state on the shortest path
		tb_state = PS[tb_state][pidx];
	}
}
