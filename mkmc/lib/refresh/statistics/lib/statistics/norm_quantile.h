#ifndef _NORM_QUANTILE_H
#define _NORM_QUANTILE_H

#include "helper_structures.h"

namespace refresh::normalization
{
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class quantile
	{
		std::vector<normalization::details::compact_histogram<ENTRY_T>> quantile_hist;
		std::vector<std::vector<std::pair<ENTRY_T, VALUE_T>>> quantile_mappings;		// TODO: consider use unordered_map here - more RAM but faster

		size_t no_series = 0;

	public:
		// *************************************************************************************
		void restart(size_t _no_series = 0)
		{
			no_series = _no_series;
			quantile_hist.clear();
			quantile_hist.resize(no_series);
			quantile_mappings.clear();
		}

		// *************************************************************************************
		template<typename Iter>
		void add(Iter first)
		{
			for (size_t i = 0; i < no_series; ++i)
				quantile_hist[i].add(*first++);
		}

		// *************************************************************************************
		void merge(quantile<ENTRY_T, VALUE_T>& to_merge)
		{
			for (size_t i = 0; i < no_series; ++i)
				quantile_hist[i].merge(to_merge.quantile_hist[i]);
		}

		// *************************************************************************************
		bool serialize(std::vector<uint8_t>& data)
		{
			// Determine mapping from input to normalized values
			std::vector<std::vector<std::tuple<ENTRY_T, size_t, VALUE_T>>> hist_cum;
			std::vector<std::pair<ENTRY_T, size_t>> lh;

			// Find cummulative histograms
			hist_cum.resize(no_series);

			for (size_t i = 0; i < no_series; ++i)
			{
				auto& hci = hist_cum[i];

				quantile_hist[i].get_histogram(lh);
				quantile_hist[i].clear();

				hci.resize(lh.size());
				hci[0] = std::make_tuple(lh[0].first, lh[0].second, (VALUE_T)lh[0].first * (VALUE_T)lh[0].second);

				for (size_t j = 1; j < hci.size(); ++j)
					hci[j] = std::make_tuple(lh[j].first, std::get<1>(hci[j - 1]) + lh[j].second, std::get<2>(hci[j - 1]) + (VALUE_T)lh[j].first * (VALUE_T)lh[j].second);
			}

			// Find mappings and serialize
			std::vector<std::pair<ENTRY_T, VALUE_T>> mapping;

			for (size_t i = 0; i < no_series; ++i)
			{
				auto& hci = hist_cum[i];
				lh.resize(hci.size());

				// Find raw histogram
				lh[0] = std::make_pair(std::get<0>(hci[0]), std::get<1>(hci[0]));

				for (size_t j = 1; j < lh.size(); ++j)
					lh[j] = std::make_pair(std::get<0>(hci[j]), std::get<1>(hci[j]) - std::get<1>(hci[j - 1]));

				mapping.clear();
				mapping.resize(lh.size());

				size_t cum_sum = 0;

				serialization::serialize_little_endian(lh.size(), data);

				for (size_t j = 0; j < lh.size(); ++j)
				{
					size_t left_side = cum_sum;
					size_t right_side = cum_sum + lh[j].second;

					VALUE_T v = ((VALUE_T)lh[j].first) * lh[j].second;

					for (size_t k = 0; k < no_series; ++k)
					{
						if (k == i)
							continue;

						auto& hck = hist_cum[k];
						auto p = std::lower_bound(hck.begin(), hck.end(), left_side, [](const auto& v, const size_t x) {return std::get<1>(v) < x; });
						auto q = std::lower_bound(hck.begin(), hck.end(), right_side, [](const auto& v, const size_t x) {return std::get<1>(v) < x; });

						v += (VALUE_T)std::get<2>(*q) - (VALUE_T)(std::get<1>(*q) - right_side) * (VALUE_T)std::get<0>(*q);
						v -= (VALUE_T)std::get<2>(*p) - (VALUE_T)(std::get<1>(*p) - left_side) * (VALUE_T)std::get<0>(*p);
					}

					cum_sum += lh[j].second;

					mapping[j] = std::make_pair(lh[j].first, v / (VALUE_T)no_series / (VALUE_T)(right_side - left_side));
				}

				for (const auto& x : mapping)
					serialization::serialize_little_endian(x, data);
			}

			return true;
		}

		// *************************************************************************************
		bool deserialize(std::vector<uint8_t>& data)
		{
			quantile_mappings.resize(no_series);
			size_t pos = 0;

			for (auto& qm : quantile_mappings)
			{
				size_t series_size;
				serialization::load_little_endian(series_size, data, pos);

				qm.resize(series_size);
				for (auto& qd : qm)
					serialization::load_little_endian(qd, data, pos);
			}

			return true;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void norm_entry(Iter1 first1, Iter1 last1, Iter2 first2)
		{
			for (size_t i = 0; i < no_series; ++i, ++first1, ++first2)
			{
				auto p = std::lower_bound(quantile_mappings[i].begin(), quantile_mappings[i].end(), *first1, [](auto entry, auto x) {return entry.first < x; });
				if (p == quantile_mappings[i].end() || p->first != *first1)
					*first2 = 0;					// TODO: should throw exception !!!

				*first2 = p->second;
			}
		}
	};
}

#endif