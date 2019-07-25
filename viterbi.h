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

#ifndef INCLUDED_VITERBI__H
#define INCLUDED_VITERBI__H

#include <algorithm>
#include <functional>
#include <limits>
#include <vector>
#include <stdexcept>

/*! A maximum likelihood decoder.
 *
 * This block implements the Viterbi algorithm in its classical form, as
 * described e.g., in: G. D. Forney, "The viterbi algorithm," in Proceedings of
 * the IEEE, vol. 61, no. 3, pp. 268-278, March 1973. doi: 10.1109/PROC.1973.9030.
 *
 * This block implements the Viterbi algorithm in its classical form,
 * as way described e.g., in \cite Forney1973.
 *
 * It takes euclidean metrics as an input and produces decoded sequences.
 */
class viterbi
{
	private:
		//! The number of possible input sequences (e.g. 2 for binary codes).
		int d_I;
		//! The number of states in the trellis.
		int	d_S;
		//! The number of possible output sequences.
		int d_O;
		/* Gives the next state ns of a branch defined by its
		 * initial state s and its input symbol i : NS[s*I+i]=ns.
		 */
		std::vector<int> d_NS;
		/* Gives the output symbol of of a branch defined by its
		 * initial state s and its input symbol i : OS[s*I+i]=os.
		 */
		std::vector<int> d_OS;
		/* Defined such that d_PS[s] contains all the previous states having a
		 * branch with state s.
		 */
		std::vector<std::vector<int> > d_PS;
		//! Defined such that d_PI[s] contains all the inputs yielding to state s.
		std::vector<std::vector<int> > d_PI;

		//! Generates PS and PI tables.
		void generate_PS_PI();

	public:
		//! Default constructor.
		viterbi();

		/*! Constructs a viterbi object.
		 * \param I The number of input sequences (e.g. 2 for binary codes).
		 * \param S The number of states in the trellis.
		 * \param O The number of output sequences (e.g. 4 for a binary code
		 *  with a coding efficiency of 1/2).
		 * \param NS Gives the next state ns of a branch defined by its
		 *  initial state s and its input symbol i : NS[s*I+i]=ns.
		 * \param OS Gives the output symbol os of a branch defined by its
		 *  initial state s and its input symbol i : OS[s*I+i]=os.
		 */
		viterbi(int I, int S, int O,
				const std::vector<int> &NS,
				const std::vector<int> &OS);

		/*! Actual Viterbi algorithm implementation.
		 *
		 * \param K Length of a block of data.
		 * \param S0 Initial state of the encoder (set to -1 if unknown).
		 * \param SK Final state of the encoder (set to -1 if unknown).
		 * \param in Input branch metrics for the algorithm.
		 * \param out Output decoded sequence.
		 */
		void viterbi_algorithm(int K, int S0, int SK,
				const float *in, unsigned int *out);

		/*! Actual Viterbi algorithm implementation.
		 *
		 * \param I The number of input sequences (e.g. 2 for binary codes).
		 * \param S The number of states in the trellis.
		 * \param O The number of output sequences (e.g. 4 for a binary code
		 *  with a coding efficiency of 1/2).
		 * \param NS Gives the next state ns of a branch defined by its
		 *  initial state s and its input symbol i : NS[s*I+i]=ns.
		 * \param OS Gives the output symbol os of a branch defined by its
		 *  initial state s and its input symbol i : OS[s*I+i]=os.
		 * \param PS Such that PS[s] contains all the previous states having a
		 *  branch with state s.
		 * \param PI Such that PI[s] contains all the inputs yielding to state s.
		 * \param K Length of a block of data.
		 * \param S0 Initial state of the encoder (set to -1 if unknown).
		 * \param SK Final state of the encoder (set to -1 if unknown).
		 * \param in Input branch metrics for the algorithm.
		 * \param out Output decoded sequence.
		 */
		static void viterbi_algorithm(int I, int S, int O,
				const std::vector<int> &NS,
				const std::vector<int> &OS,
				const std::vector< std::vector<int> > &PS,
				const std::vector< std::vector<int> > &PI,
				int K, int S0, int SK,
				const float *in, unsigned int *out);

		//! Getter for d_I.
		int get_I() { return d_I; }
		//! Getter for d_S.
		int get_S() { return d_S; }
		//! Getter for d_O.
		int get_O() { return d_O; }
		//! Getter for d_NS.
		std::vector<int>& get_NS() { return d_NS; }
		//! Getter for d_OS.
		std::vector<int>& get_OS() { return d_OS; }
};

#endif /* INCLUDED_VITERBI_H */
