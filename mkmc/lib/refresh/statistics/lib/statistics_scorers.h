#ifndef _STATISTICS_SCORERS_H
#define _STATISTICS_SCORERS_H

#include <algorithm>
#include <vector>
#include <cinttypes>
#include <cmath>

#define STATS_ENABLE_STDVEC_WRAPPERS
#include "../../../../stats/include/stats.hpp"

namespace refresh
{
	template<typename CLASS_T = uint32_t, typename VALUE_T = double>
	class scorers
	{
	public:
		struct anova_t
		{
			double statistic;
			double p_value;
		};

	private:
		std::vector<VALUE_T> maxes;
		std::vector<VALUE_T> sum_sqrts;

		std::vector<size_t> anova_n_items;
		std::vector<double> anova_sums;

	public:
		template<typename X_Iter, typename C_Iter>
		VALUE_T dids(X_Iter X_first, X_Iter X_last, C_Iter C_first, size_t no_classes)
		{
			return dids_n(X_first, C_first, no_classes, std::distance(X_first, X_last));
		}

		template<typename X_Iter, typename C_Iter>
		VALUE_T dids_n(X_Iter X_first, C_Iter C_first, size_t no_classes, size_t n)
		{
			maxes.clear();
			maxes.resize(no_classes, (VALUE_T) 0);

			sum_sqrts.clear();
			sum_sqrts.resize(no_classes, (VALUE_T) 0);

			auto p_C = C_first;
			auto p_X = X_first;

			// Find maximum for each class
			for (size_t i = 0; i < n; ++i, p_C++, ++p_X)
				if ((VALUE_T) *p_X > maxes[*p_C])
					maxes[*p_C] = (VALUE_T) *p_X;

			// Calculate sums of sqrts for "large" values
			p_C = C_first;
			p_X = X_first;

			for (size_t i = 0; i < n; ++i, ++p_C, ++p_X)
				for(size_t j = 0; j < no_classes; ++j)
					if (j != *p_C)
						if ((VALUE_T) *p_X > maxes[j])
							sum_sqrts[j] += sqrt((VALUE_T)*p_X - maxes[j]);

			return *std::max_element(sum_sqrts.begin(), sum_sqrts.end());
		}


		template<typename X_Iter, typename C_Iter>
		anova_t anova(X_Iter X_first, X_Iter X_last, C_Iter C_first, size_t no_classes)
		{
			return anova_n(X_first, C_first, no_classes, std::distance(X_first, X_last));
		}

		template<typename X_Iter, typename C_Iter>
		anova_t anova_n(X_Iter X_first, C_Iter C_first, size_t no_classes, size_t n)
		{
			anova_t ret{ 0, 0 };

			if (n == 0 || no_classes == 0)
				return ret;

			anova_n_items.resize(no_classes);
			anova_sums.resize(no_classes);

			std::fill_n(anova_n_items.begin(), no_classes, 0);
			std::fill_n(anova_sums.begin(), no_classes, 0.0);

			double X_sum = 0;
			double X_sum2 = 0;
			auto p_X = X_first;
			auto p_C = C_first;

			for (size_t i = 0; i < n; ++i, ++p_X, ++p_C)
			{
				X_sum += *p_X;
				X_sum2 += stats_details::pow2(*p_X);
				
				anova_sums[*p_C] += *p_X;
				anova_n_items[*p_C]++;
			}

			double sstot = X_sum2 - stats_details::pow2(X_sum) / (double)n;
			
			double X_avg = X_sum / n;
			double ssbg = 0.0;

			for (size_t i = 0; i < no_classes; ++i)
				ssbg += anova_n_items[i] * stats_details::pow2(anova_sums[i] / anova_n_items[i] - X_avg);

			double sswg = sstot - ssbg;
			double dfbg = (double)no_classes - 1.0;
			double dfwg = (double)n - (double)no_classes;

			double msb = ssbg / dfbg;
			double msw = sswg / dfwg;
			double f = msb / msw;

			ret.statistic = f;
			ret.p_value = std::clamp<double>(1.0 - stats::pf(ret.statistic, dfbg, dfwg, false), 0, 1);

			return ret;
		}
	};
}


#endif