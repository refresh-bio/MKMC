#ifndef _STATISTICS_ENTROPY_H
#define _STATISTICS_ENTROPY_H

#include <algorithm>
#include <numeric>
#include <cmath>

#include "statistics/helper_structures.h"

namespace refresh
{
	class statistics_entropy
	{
		normalization::details::fast_log<size_t> flog;
		double ln2;

	public:
		// *************************************************************************************
		statistics_entropy() : flog(1000)
		{
			ln2 = flog.log(2);
		}

		// *************************************************************************************
		template<typename Iter>
		double entropy(Iter first, Iter last)
		{
			return entropy_n(first, std::distance(first, last));
		}

		// *************************************************************************************
		template<typename Iter>
		double entropy_n(Iter first, size_t n)
		{
			size_t sum_C = 0;
			double sum_Ci_ln_Ci = 0;

			for (size_t i = 0; i < n; ++i, ++first)
			{
				sum_C += *first;
				sum_Ci_ln_Ci += (double) *first * flog.log(*first);
			}

			double H_k_e = -sum_Ci_ln_Ci / sum_C + flog.log(sum_C);
			double H_k = H_k_e / ln2;

			return H_k;
		}
	};
}

#endif