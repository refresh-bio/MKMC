#ifndef _STATISTICS_SCORERS_H
#define _STATISTICS_SCORERS_H

#include <algorithm>
#include <vector>
#include <cinttypes>
#include <cmath>

namespace refresh
{
	template<typename CLASS_T = uint32_t, typename VALUE_T = double>
	class scorers
	{
		std::vector<VALUE_T> maxes;
		std::vector<VALUE_T> sum_sqrts;

	public:
		template<typename Iter1, typename Iter2>
		VALUE_T dids(Iter1 class_first, Iter1 class_last, Iter2 counter_first, size_t no_classes)
		{
			return dids_n(class_first, counter_first, no_classes, std::distance(class_first, class_last));
		}

		template<typename Iter1, typename Iter2>
		VALUE_T dids_n(Iter1 class_first, Iter2 counter_first, size_t no_classes, size_t n)
		{
			maxes.clear();
			maxes.resize(no_classes, (VALUE_T) 0);

			sum_sqrts.clear();
			sum_sqrts.resize(no_classes, (VALUE_T) 0);

			auto p_class = class_first;
			auto p_counter = counter_first;

			// Find maximum for each class
			for (size_t i = 0; i < n; ++i, p_class++, ++p_counter)
				if ((VALUE_T) *p_counter > maxes[*p_class])
					maxes[*p_class] = (VALUE_T) *p_counter;

			// Calculate sums of srqts for "large" values
			p_class = class_first;
			p_counter = counter_first;

			for (size_t i = 0; i < n; ++i, ++p_class, ++p_counter)
				for(size_t j = 0; j < no_classes; ++j)
					if (j != *p_class)
						if ((VALUE_T) *p_counter > maxes[j])
							sum_sqrts[j] += sqrt((VALUE_T)*p_counter - maxes[j]);

			return *std::max_element(sum_sqrts.begin(), sum_sqrts.end());
		}

	};
}


#endif