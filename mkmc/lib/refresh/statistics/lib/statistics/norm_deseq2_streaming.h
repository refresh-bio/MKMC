#ifndef _NORM_DESEQ2_STREAMING_H
#define _NORM_DESEQ2_STREAMING_H

#include <algorithm>
#include "helper_structures.h"

namespace refresh::normalization
{
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class deseq2_streaming
	{
		normalization::details::fast_log<ENTRY_T> flog;

		size_t no_series = 0;
		std::vector<std::vector<double>> series;
		std::vector<VALUE_T> mult_factor;

		// *************************************************************************************
		double median(std::vector<double>& vec)
		{
			size_t n = vec.size();

			if (n == 0)
				return 0;

			auto p_middle = vec.begin() + n / 2;
			nth_element(vec.begin(), p_middle, vec.end());

			if (n % 2 == 0)
			{
				auto p = max_element(vec.begin(), p_middle);
				return (*p_middle + *p) / 2;
			}
			else
				return *p_middle;
		}

	public:
		// *************************************************************************************
		void restart(size_t _no_series = 0)
		{
			no_series = _no_series;
			series.clear();
			series.resize(no_series);
			mult_factor.clear();
		}

		// *************************************************************************************
		template<typename Iter>
		void add(Iter first)
		{
			double avg = 0;

			for (size_t i = 0; i < no_series; ++i)
			{
				double x = flog.log(*first++);

				series[i].emplace_back(x);
				avg += x;
			}

			avg /= no_series;

			for (size_t i = 0; i < no_series; ++i)
				series[i].back() -= avg;
		}

		// *************************************************************************************
		void merge(deseq2_streaming<ENTRY_T, VALUE_T>& to_merge)
		{
			for (size_t i = 0; i < no_series; ++i)
				series[i].insert(series[i].end(), to_merge.series[i].begin(), to_merge.series[i].end());
		}

		// *************************************************************************************
		bool serialize(std::vector<uint8_t>& data)
		{
			// Determine multiplicator factors
			mult_factor.resize(no_series);

			for (size_t i = 0; i < mult_factor.size(); ++i)
			{
				double med = median(series[i]);

				mult_factor[i] = 1.0 / med;
			}

			// Serialize
			for (auto x : mult_factor)
				serialization::serialize_little_endian(x, data);

			return true;
		}

		// *************************************************************************************
		bool deserialize(std::vector<uint8_t>& data)
		{
			series.clear();
			mult_factor.resize(no_series);

			size_t pos = 0;

			for (auto& x : mult_factor)
				serialization::load_little_endian(x, data, pos);

			return true;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void norm_entry(Iter1 first1, Iter1 last1, Iter2 first2)
		{
			for (size_t i = 0; i < no_series; ++i, ++first1, ++first2)
				*first2 = (VALUE_T)*first1 * mult_factor[i];
		}
	};
}
#endif