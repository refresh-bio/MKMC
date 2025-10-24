#ifndef _NORM_FREQUENCY_COUNT_H
#define _NORM_FREQUENCY_COUNT_H

namespace refresh::normalization
{
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class frequency_count
	{
		std::vector<ENTRY_T> counts;
		std::vector<VALUE_T> mult_factor;
		size_t no_series = 0;

	public:
		// *************************************************************************************
		void restart(size_t _no_series = 0)
		{
			no_series = _no_series;
			counts.clear();
			counts.resize(no_series);
			mult_factor.clear();
		}

		// *************************************************************************************
		template<typename Iter>
		void add(Iter first)
		{
			for (size_t i = 0; i < no_series; ++i)
				counts[i] += *first++;
		}

		// *************************************************************************************
		void merge(frequency_count<ENTRY_T, VALUE_T>& to_merge)
		{
			for (size_t i = 0; i < no_series; ++i)
				counts[i] += to_merge.counts[i];
		}

		// *************************************************************************************
		bool serialize(std::vector<uint8_t>& data)
		{
			// Determine multiplicator factors
			mult_factor.resize(no_series);

			for (size_t i = 0; i < mult_factor.size(); ++i)
				if (counts[i])
					mult_factor[i] = (VALUE_T)1.0 / (VALUE_T)counts[i];
				else
					mult_factor[i] = (VALUE_T)0.0;

			// Serialize
			for (auto x : mult_factor)
				serialization::serialize_little_endian(x, data);

			return true;
		}

		// *************************************************************************************
		bool deserialize(std::vector<uint8_t>& data)
		{
			counts.clear();
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