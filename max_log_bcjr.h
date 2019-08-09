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

#ifndef INCLUDED_TURBO_MAX_LOG_BCJR_H
#define INCLUDED_TURBO_MAX_LOG_BCJR_H

#include "log_bcjr_base.h"

/*!
* \brief <+description+>
*
*/
class max_log_bcjr : public log_bcjr_base
{
	public:
		//! Default constructor.
		max_log_bcjr();

		/*! Constructs a log_bcjr object.
		 * \param I The number of input sequences (e.g. 2 for binary codes).
		 * \param S The number of states in the trellis.
		 * \param O The number of output sequences (e.g. 4 for a binary code
		 *  with a coding efficiency of 1/2).
		 * \param NS Gives the next state ns of a branch defined by its
		 *  initial state s and its input symbol i : NS[s*I+i]=ns.
		 * \param OS Gives the output symbol os of a branch defined by its
		 *  initial state s and its input symbol i : OS[s*I+i]=os.
		 */
		max_log_bcjr(int I, int S, int O,
				const std::vector<int> &NS,
				const std::vector<int> &OS) : log_bcjr_base(I, S, O, NS, OS) {};

		//! Computes max of two value.
		/*!
		 * \param A First operand.
		 * \param B Second operand.
		 *
		 * \return max(A,B).
		 */
		static inline float max(float A, float B)
		{
			return std::max(A, B);
		}
		// Override log_bcjr_base method
		float _max_star(float A, float B) { return max(A, B); }

		//! Compute max of a vector.
		/*!
		 * \param vec Input data.
		 * \param n_ele number of elements in the vector.
		 *
		 * \return: max of vec.
		 */
		static inline float max(const float *vec, size_t n_ele)
		{
			return *std::max_element(vec, vec+n_ele);
		}
		// Override log_bcjr_base method
		float _max_star(const float *vec, size_t n_ele) { return max(vec, n_ele); }
};

#endif /* INCLUDED_TURBO_MAX_LOG_BCJR_H */
