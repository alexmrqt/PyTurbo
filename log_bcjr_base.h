/* -*- c++ -*- */
/*
 * Copyright 2020 Alexandre Marquet.
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


#ifndef INCLUDED_TURBO_LOG_BCJR_base_H
#define INCLUDED_TURBO_LOG_BCJR_base_H

#include <algorithm>
#include <limits>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cfloat>

/*!
* \brief <+description+>
*
*/
class log_bcjr_base
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
		
        /* Same as d_FSM.OS(), but re-ordered in the following way:
         * d_ordered_OS[s*I+i] = d_OS()[d_PS()[s][i]*I + d_PI()[s][i]]
		 */
        std::vector<int> d_ordered_OS;

		/* Defined such that d_PS[s] contains all the previous states having a
		 * branch with state s.
		 * Such a previous state may appear multiple time if there are multiple
		 * transistions between two states.
		 */
		std::vector<std::vector<int> > d_PS;

		//! Defined such that d_PI[s] contains all the inputs yielding to state s.
		std::vector<std::vector<int> > d_PI;

		//! Generates PS, PI and T tables.
		void generate_PS_PI();

	public:
		/*! Constructs a log_bcjr_base object.
		 * \param I The number of input sequences (e.g. 2 for binary codes).
		 * \param S The number of states in the trellis.
		 * \param O The number of output sequences (e.g. 4 for a binary code
		 *  with a coding efficiency of 1/2).
		 * \param NS Gives the next state ns of a branch defined by its
		 *  initial state s and its input symbol i : NS[s*I+i]=ns.
		 * \param OS Gives the output symbol os of a branch defined by its
		 *  initial state s and its input symbol i : OS[s*I+i]=os.
		 */
		log_bcjr_base(int I, int S, int O,
				const std::vector<int> &NS,
				const std::vector<int> &OS);

		//! Computes max* of two value.
		/*!
		 * \param A First operand.
		 * \param B Second operand.
		 *
		 * \return max*(A,B).
		 */
		virtual float _max_star(float A, float B) = 0;

		//! Recursively compute max* of a vector.
		/*!
		 * To compute max*(A,B,C,...), recursive calls to max* are performed.
		 * For instance: max*(A,B,C) = max*(max*(A,B),C).
		 *
		 * \param vec Input data.
		 * \param n_ele number of elements in the vector.
		 *
		 * \return: max* of vec.
		 * If axis is given, the result is an array of dimension vec.ndim - 1.
		 */
		virtual float _max_star(const float *vec, size_t n_ele) = 0;

		//! Compute forward log metrics.
		/*!
		 * From A_k(s) the forward log metric for state s at time index k, and
		 * G_k(s,i) the log metric of the branch identified by state s and
		 * input symbol i at index k, this function computes:
		 *
		 * A_k(s) = max*_{ s' \in [0 ; d_S[, i \in \tau(s',s) } G_{k-1}(s', i) + A_{k-1}(s')
		 *
		 * where \tau(s,s') regroups every input symbols that belongs to every
		 * transitions between s and s'.
		 *
		 * Note: in practice, here, we only have the metrics of every possible
		 * output symbols: G_k(o) with o \in [0 ; d_O[. The correspondance is
		 * done through d_OS: G_k(s,i) = G_k(d_OS[s*I+i]).
		 *
		 * \param G Const reference to the log metrics vector (size: d_O*K).
		 * \param A0 Const reference to the initial forward state metrics
		 * (size: d_S).
		 * \param A Reference to the forward metrics vector (will have a size
		 * of d_S*(K+1) at the end of function execution).
		 * \param K Number of observations.
		 */
		virtual void compute_fw_metrics(const std::vector<float> &G,
				const std::vector<float> &A0, std::vector<float> &A, size_t K);

		//! Compute backward log metrics.
		/*!
		 * From B_k(s) the backward log metric for state s at time index k, and
		 * G_k(s,i) the log metric of the branch identified by state s and
		 * input symbol i at index k, this function computes:
		 *
		 * B_k(s) = max*_{ s' \in [0 ; d_S[, i \in \tau(s,s') } G_k(s, i) + B_{k+1}(s').
		 *
		 * where \tau(s,s') regroups every input symbols that belongs to every
		 * transitions between s and s'.
		 *
		 * Note: in practice, here, we only have the metrics of every possible
		 * output symbols: G_k(o) with o \in [0 ; d_O[. The correspondance is
		 * done through d_OS: G_k(s,i) = G_k(d_OS[s*I+i]).
		 *
		 * \param G Const reference to the log metrics vector (size: d_O*K).
		 * \param BK Const reference to the final backward state metrics
		 * (size: d_S).
		 * \param B Reference to the backward metrics vector (will have a size
		 * of d_S*(K+1) at the end of function execution).
		 * \param K Number of observations.
		 */
		virtual void compute_bw_metrics(const std::vector<float> &G,
				const std::vector<float> &BK, std::vector<float> &B, size_t K);

		//! Compute branch log a-posteriori probabilities.
		/*!
		 * From A_k(s) the forward log metric for state s at time index k,
		 * B_k(s) the backward log metric for state s at time index k, and
		 * G_k(s,s') the branch log metric between states s' and s at time
		 * index k, this function computes:
		 *
		 * APP_k(s,i) = B_{k+1}(NS(s,i)) + G_k(s,i) + A_k(s),
		 *
		 * where s' = NS(s,i) is the next state for transition with initial
		 * state s and input symbol i (NS[s*I+i]).
		 * Which is equivalent to log a-posteriori probabilites, up to an
		 * additive constant.
		 *
		 * \param A Const reference to the forward metrics vector (size: d_S*(K+1)).
		 * \param B Const reference to the backward metrics vector (size: d_S*(K+1)).
		 * \param G Const reference to the branch log metrics vector (size: d_O*K).
		 * \param K Number of observations.
		 * \param out Reference to a posteriori branch log probabilities (will
		 *  have a size of d_S*d_I*K at the end of function execution).
		 *
		 */
		virtual void compute_app(const std::vector<float> &A,
				const std::vector<float> &B, const std::vector<float> &G,
				size_t K, std::vector<float> &out);

		/*! Actually computes logarithm of a-posteriori probabilities for a
		 * given observation sequence.
		 *
		 * \param A0 Log of initial state probabilities of the encoder (size: d_S).
		 * \param BK Log of final state probabilities of the encoder (size: d_S).
		 * \param in Log of input branch metrics for the algorithm (size: d_O*k).
		 * \param out A quantity equivalent to log a-posteriori probabilites, up
		 *  to an additive constant (will have a size of d_S*d_I*K at the end of
		 *  function execution).
		 */
		void log_bcjr_algorithm(const std::vector<float> &A0,
				const std::vector<float> &BK,
				const std::vector<float> &in,
				std::vector<float> &out);

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

#endif /* INCLUDED_TURBO_LOG_BCJR_base_H */

