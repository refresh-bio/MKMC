#ifndef _STATISTICS_NORMALIZATION
#define _STATISTICS_NORMALIZATION

#include <vector>
#include <map>
#include <cinttypes>
#include <algorithm>
#include <tuple>
#include <bit>

#include "../../serialization/lib/serialization.h"

namespace refresh
{
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_base
	{
	public:
		enum class method_t {frequency_count = 0, deseq = 1, deseq2 = 2, quantile = 3};

	protected:
		uint64_t methods = 0;
		size_t no_series = 0;
		size_t no_entries = 0;


		// *************************************************************************************
		bool method_included(method_t m)
		{
			return (bool)(methods & (1ull << (uint32_t) m));
		}
	 

	public:
		// *************************************************************************************
		normalization_base() = default;

		// *************************************************************************************
		size_t get_no_series()
		{
			return no_series;
		}

		// *************************************************************************************
		void clear()
		{
			methods = 0;
			no_series = 0;
			no_entries = 0;
		}

		// *************************************************************************************
		void register_method(method_t m)
		{
			methods |= 1ull << (uint32_t)m;
		}

		// *************************************************************************************
		void set_no_series(uint64_t x)
		{
			no_series = x;
		}

	};

	// *************************************************************************************
	namespace normalization::details
	{
		template<typename ENTRY_T>
		class compact_histogram
		{
			std::map<ENTRY_T, size_t> data;

			void add(size_t x, size_t inc)
			{
				data[x] += inc;
			}

		public:
			void add(size_t x)
			{
				data[x]++;
			}

			void get_histogram(std::vector<std::pair<ENTRY_T, size_t>>& hist)
			{
				hist.assign(data.begin(), data.end());
			}

			void clear()
			{
				data.clear();
			}

			void merge(const compact_histogram<ENTRY_T>& ch)
			{
				for (const auto p : ch.hist)
					add(p.first, p.second);
			}
		};

		template<>
		class compact_histogram<size_t>
		{
			const size_t thr = 32 << 10;
//			const size_t thr = 8;				// for testing

			std::vector<size_t> small;
			std::map<size_t, size_t> big;

			void add(size_t x, size_t inc)
			{
				if (x >= thr)
					big[x] += inc;
				else
				{
					if (x >= small.size())
						small.resize(std::min(2 * (x + 1), thr));
					small[x] += inc;
				}
			}

		public:
			void add(size_t x)
			{
				if (x >= thr)
					big[x]++;
				else
				{
					if (x >= small.size())
						small.resize(std::min(2 * (x + 1), thr));
					small[x]++;
				}
			}

			void get_histogram(std::vector<std::pair<size_t, size_t>>& hist)
			{
				hist.clear();

				for (size_t i = 0; i < small.size(); ++i)
					if (small[i])
						hist.emplace_back(i, small[i]);

				for (const auto x : big)
					hist.emplace_back(x.first, x.second);

//				hist.shrink_to_fit();
			}

			void clear()
			{
				small.clear();
				small.shrink_to_fit();
				big.clear();
			}

			void merge(const compact_histogram<size_t>& ch)
			{
				for (size_t i = 0; i < ch.small.size(); ++i)
					if (ch.small[i] != 0)
						add(i, ch.small[i]);
				for (const auto p : ch.big)
					add(p.first, p.second);
			}
		};
	}

	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_learn : public normalization_base<ENTRY_T, VALUE_T>
	{
		using base = normalization_base<ENTRY_T, VALUE_T>;

		std::vector<ENTRY_T> counts;
		std::vector<normalization::details::compact_histogram<ENTRY_T>> quantile_hist;

		// *************************************************************************************
		void restart_method_frequency_count()
		{
			counts.clear();
			counts.resize(base::no_series);
		}

		// *************************************************************************************
		void restart_method_quantile()
		{
			quantile_hist.clear();
			quantile_hist.resize(base::no_series);
		}

		// *************************************************************************************
		template<typename Iter>
		void add_method_frequency_count(Iter first)
		{
			for (size_t i = 0; i < base::no_series; ++i)
				counts[i] += *first++;
		}

		// *************************************************************************************
		template<typename Iter>
		void add_method_quantile(Iter first)
		{
			for (size_t i = 0; i < base::no_series; ++i)
				quantile_hist[i].add(*first++);
		}

		// *************************************************************************************
		void merge_method_frequency_count(normalization_learn<ENTRY_T, VALUE_T>& to_merge)
		{
			for (size_t i = 0; i < base::no_series; ++i)
				counts[i] += to_merge.counts[i];
		}

		// *************************************************************************************
		void merge_method_quantile(normalization_learn<ENTRY_T, VALUE_T>& to_merge)
		{
			for (size_t i = 0; i < base::no_series; ++i)
				quantile_hist[i].merge(to_merge.quantile_hist[i]);
		}

		// *************************************************************************************
		bool serialize_method_frequency_count(std::vector<uint8_t>& data)
		{
			// Determine multiplicator factors
			std::vector<VALUE_T> mult_factor(base::no_series);

			for (size_t i = 0; i < mult_factor.size(); ++i)
				mult_factor[i] = (VALUE_T)1.0 / (VALUE_T) counts[i];

			// Serialize
			for (auto x : mult_factor)
				serialization::serialize_little_endian(x, data);

			return true;
		}

		// *************************************************************************************
		bool serialize_method_quantile(std::vector<uint8_t>& data)
		{
			// Determine mapping from input to normalized values
			std::vector<std::vector<std::tuple<ENTRY_T, size_t, VALUE_T>>> hist_cum;
			std::vector<std::pair<ENTRY_T, size_t>> lh;

			// Find cummulative histograms
			hist_cum.resize(base::no_series);

			for (size_t i = 0; i < base::no_series; ++i)
			{
				auto& hci = hist_cum[i];

				quantile_hist[i].get_histogram(lh);
				quantile_hist[i].clear();

				hci.resize(lh.size());
				hci[0] = std::make_tuple(lh[0].first, lh[0].second, (VALUE_T)lh[0].first * (VALUE_T)lh[0].second);

				for (size_t j = 1; j < hci.size(); ++j)
					hci[j] = std::make_tuple(lh[j].first, get<1>(hci[j-1]) + lh[j].second, get<2>(hci[j-1]) + (VALUE_T)lh[j].first * (VALUE_T)lh[j].second);
			}

			// Find mappings and serialize
			std::vector<std::pair<ENTRY_T, VALUE_T>> mapping;

			for (size_t i = 0; i < base::no_series; ++i)
			{
				auto& hci = hist_cum[i];
				lh.resize(hci.size());

				// Find raw histogram
				lh[0] = std::make_pair(get<0>(hci[0]), get<1>(hci[0]));
				
				for (size_t j = 1; j < lh.size(); ++j)
					lh[j] = std::make_pair(get<0>(hci[j]), get<1>(hci[j]) - get<1>(hci[j-1]));

				mapping.clear();
				mapping.resize(lh.size());

				size_t cum_sum = 0;

				serialization::serialize_little_endian(lh.size(), data);

				for (size_t j = 0; j < lh.size(); ++j)
				{
					size_t left_side = cum_sum;
					size_t right_side = cum_sum + lh[j].second;

					VALUE_T v = ((VALUE_T) lh[j].first) * lh[j].second;

					for (size_t k = 0; k < base::no_series; ++k)
					{
						if (k == i)
							continue;

						auto& hck = hist_cum[k];
						auto p = std::lower_bound(hck.begin(), hck.end(), left_side, [](const auto& v, const size_t x) {return get<1>(v) < x; });
						auto q = std::lower_bound(hck.begin(), hck.end(), right_side, [](const auto& v, const size_t x) {return get<1>(v) < x; });
						
						v += (VALUE_T)get<2>(*q) - (VALUE_T)(get<1>(*q) - right_side) * (VALUE_T)get<0>(*q);
						v -= (VALUE_T)get<2>(*p) - (VALUE_T)(get<1>(*p) - left_side) * (VALUE_T)get<0>(*p);
					}

					cum_sum += lh[j].second;

					mapping[j] = std::make_pair(lh[j].first, v / (VALUE_T) base::no_series / (VALUE_T) (right_side - left_side));
				}

				for (const auto x : mapping)
					serialization::serialize_little_endian(x, data);
			}

			return true;
		}


	public:
		// *************************************************************************************
		normalization_learn() = default;

		// *************************************************************************************
		void initialize()
		{
			if (base::method_included(base::method_t::frequency_count))
				restart_method_frequency_count();
			
			if (base::method_included(base::method_t::quantile))
				restart_method_quantile();
		}
		 
		// *************************************************************************************
		bool serialize(const normalization_base<ENTRY_T, VALUE_T>::method_t m, std::vector<uint8_t>& data)
		{
			if (!base::method_included(m) && base::no_entries == 0)
				return false;

			data.clear();

			if (m == base::method_t::frequency_count)
				return serialize_method_frequency_count(data);
			if (m == base::method_t::quantile)
				return serialize_method_quantile(data);
		}

		// *************************************************************************************
		template<typename Iter>
		void add_entry(Iter first, Iter last)
		{
			if (base::method_included(base::method_t::frequency_count))
				add_method_frequency_count(first);
			if (base::method_included(base::method_t::quantile))
				add_method_quantile(first);
		}

		// *************************************************************************************
		template<typename Container>
		void add_entry(const Container &cont)
		{
			//assert(cont.size() == normalization_base<ENTRY_T, VALUE_T>::method_t::no_entries);

			base::no_entries++;

			if (base::method_included(base::method_t::frequency_count))
				add_method_frequency_count(cont.begin());
			if (base::method_included(base::method_t::quantile))
				add_method_quantile(cont.begin());
		}

		template<typename Iter>
		bool merge_with(Iter first, Iter last)
		{
			if (first == last)
				return true;

			// Check if possible
			for (auto p = first; p != last; ++p)
				if (base::no_series != p->no_series || base::methods != p->methods)
					return false;

			for (auto p = first; p != last; ++p)
			{
				base::no_entries += p->no_entries;

				if(base::method_included(base::method_t::frequency_count))
					merge_method_frequency_count(*p);
				if(base::method_included(base::method_t::quantile))
					merge_method_quantile(*p);
			}
		}
	};
	
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_work : public normalization_base<ENTRY_T, VALUE_T>
	{
		using base = normalization_base<ENTRY_T, VALUE_T>;

		std::vector<VALUE_T> frequency_count_mult_factor;
		std::vector<std::vector<std::pair<ENTRY_T, VALUE_T>>> quantile_mappings;		// TODO: consider use unordered_map here - more RAM but faster

		// *************************************************************************************
		bool deserialize_method_frequency_count(std::vector<uint8_t>& data)
		{
			frequency_count_mult_factor.resize(base::no_series);

			size_t pos = 0;

			for (auto& x : frequency_count_mult_factor)
				serialization::load_little_endian(x, data, pos);

			return true;
		}

		// *************************************************************************************
		bool deserialize_method_quantile(std::vector<uint8_t>& data)
		{
			quantile_mappings.resize(base::no_series);
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
		void norm_entry_frequency_count(Iter1 first1, Iter1 last1, Iter2 first2)
		{
			for (size_t i = 0; i < base::no_series; ++i, ++first1, ++first2)
				*first2 = (VALUE_T)*first1 * frequency_count_mult_factor[i];
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void norm_entry_quantile(Iter1 first1, Iter1 last1, Iter2 first2)
		{
			for (size_t i = 0; i < base::no_series; ++i, ++first1, ++first2)
			{
				auto p = std::lower_bound(quantile_mappings[i].begin(), quantile_mappings[i].end(), *first1, [](auto entry, auto x) {return entry.first < x; });
				if (p == quantile_mappings[i].end() || p->first != *first1)
					*first2 = 0;					// TODO: should throw exception !!!

				*first2 = p->second;
			}
		}


	public:
		// *************************************************************************************
		normalization_work() = default;

		// *************************************************************************************
		void initialize()
		{
		}

		// *************************************************************************************
		bool deserialize(const normalization_base<ENTRY_T, VALUE_T>::method_t m, std::vector<uint8_t>& data)
		{
			if (!base::method_included(m))
				return false;

			if (m == base::method_t::frequency_count)
				return deserialize_method_frequency_count(data);
			if (m == base::method_t::quantile)
				return deserialize_method_quantile(data);

		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void norm_entry(normalization_base<ENTRY_T, VALUE_T>::method_t method, Iter1 first1, Iter1 last1, Iter2 first2)
		{
			if (!base::method_included(method))
				return;

			if (method == base::method_t::frequency_count)
				norm_entry_frequency_count(first1, last1, first2);
			if (method == base::method_t::quantile)
				norm_entry_quantile(first1, last1, first2);
		}

		// *************************************************************************************
		template<typename Cont1, typename Cont2>
		void norm_entry(normalization_base<ENTRY_T, VALUE_T>::method_t method, const Cont1& cont1, Cont2& cont2)
		{
			cont2.resize(cont1.size());
			norm_entry(method, cont1.begin(), cont1.end(), cont2.begin());
		}
	};


}

#endif