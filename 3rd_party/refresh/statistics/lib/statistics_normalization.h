#ifndef _STATISTICS_NORMALIZATION
#define _STATISTICS_NORMALIZATION

#include <vector>
#include <map>
#include <cinttypes>
#include <algorithm>
#include <tuple>
#include <bit>

#include "refresh/serialization/lib/serialization.h"
#include "statistics/helper_structures.h"
#include "statistics/norm_frequency_count.h"
#include "statistics/norm_quantile.h"
#include "statistics/norm_deseq2.h"
//#include "statistics/norm_deseq2_streaming.h"

namespace refresh
{
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_base
	{
	public:
//		enum class method_t { frequency_count = 0, deseq = 1, quantile = 2, deseq2 = 3, deseq2_streaming = 4 };
		enum class method_t { frequency_count = 0, deseq = 1, quantile = 2, deseq2 = 3, deseq2_fc = 4};

	protected:
		uint64_t methods = 0;
		size_t no_series = 0;
		size_t no_entries = 0;

		normalization::frequency_count<ENTRY_T, VALUE_T> meth_frequency_count;
		normalization::quantile<ENTRY_T, VALUE_T> meth_quantile;
		normalization::deseq2<ENTRY_T, VALUE_T> meth_deseq2;
//		normalization::deseq2_streaming<ENTRY_T, VALUE_T> meth_deseq2_streaming;

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
			meth_frequency_count.restart(no_series);
			meth_quantile.restart(no_series);
			meth_deseq2.restart(no_series);
//			meth_deseq2_streaming.restart(no_series);
		}
	};

	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_learn : public normalization_base<ENTRY_T, VALUE_T>
	{
		using base = normalization_base<ENTRY_T, VALUE_T>;

	public:
		// *************************************************************************************
		normalization_learn() = default;

		// *************************************************************************************
		void initialize()
		{
			if (base::method_included(base::method_t::frequency_count))
				base::meth_frequency_count.restart(base::no_series);
			if (base::method_included(base::method_t::quantile))
				base::meth_quantile.restart(base::no_series);
			if (base::method_included(base::method_t::deseq2) || base::method_included(base::method_t::deseq2_fc))
				base::meth_deseq2.restart(base::no_series);
//			if (base::method_included(base::method_t::deseq2_streaming))
//				base::meth_deseq2_streaming.restart(base::no_series);
		}
		 
		// *************************************************************************************
		bool serialize(const typename normalization_base<ENTRY_T, VALUE_T>::method_t m, std::vector<uint8_t>& data)
		{
			if (!base::method_included(m) && base::no_entries == 0)
				return false;

			data.clear();

			if (m == base::method_t::frequency_count)
				return base::meth_frequency_count.serialize(data);
			if (m == base::method_t::quantile)
				return base::meth_quantile.serialize(data);
			if (m == base::method_t::deseq2 || m == base::method_t::deseq2_fc)
				return base::meth_deseq2.serialize(data);
//			if (m == base::method_t::deseq2_streaming)
//				return base::meth_deseq2_streaming.serialize(data);

			assert(false);

			return false;
		}

		// *************************************************************************************
		template<typename Iter>
		void add_entry(Iter first, Iter last)
		{
			if (base::method_included(base::method_t::frequency_count))
				base::meth_frequency_count.add(first);
			if (base::method_included(base::method_t::quantile))
				base::meth_quantile.add(first);
			if (base::method_included(base::method_t::deseq2))
				base::meth_deseq2.add(first);
//			if (base::method_included(base::method_t::deseq2_streaming))
//				base::meth_deseq2_streaming.add(first);
		}

		// *************************************************************************************
		template<typename Container>
		void add_entry(const Container &cont)
		{
			//assert(cont.size() == normalization_base<ENTRY_T, VALUE_T>::method_t::no_entries);

			base::no_entries++;

			if (base::method_included(base::method_t::frequency_count))
				base::meth_frequency_count.add(cont.begin());
			if (base::method_included(base::method_t::quantile))
				base::meth_quantile.add(cont.begin());
			if (base::method_included(base::method_t::deseq2))
				base::meth_deseq2.add(cont.begin());
//			if (base::method_included(base::method_t::deseq2_streaming))
//				base::meth_deseq2_streaming.add(cont.begin());
		}

		// *************************************************************************************
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
					base::meth_frequency_count.merge(p->meth_frequency_count);
				if(base::method_included(base::method_t::quantile))
					base::meth_quantile.merge(p->meth_quantile);
				if (base::method_included(base::method_t::deseq2))
					base::meth_deseq2.merge(p->meth_deseq2);
//				if (base::method_included(base::method_t::deseq2_streaming))
//					base::meth_deseq2_streaming.merge(p->meth_deseq2_streaming);
			}

			return true;
		}
	};
	
	// *************************************************************************************
	//
	// *************************************************************************************
	template<typename ENTRY_T, typename VALUE_T>
	class normalization_work : public normalization_base<ENTRY_T, VALUE_T>
	{
		using base = normalization_base<ENTRY_T, VALUE_T>;

	public:
		// *************************************************************************************
		normalization_work() = default;

		// *************************************************************************************
		void initialize()
		{
		}

		// *************************************************************************************
		bool deserialize(const typename normalization_base<ENTRY_T, VALUE_T>::method_t m, std::vector<uint8_t>& data)
		{
			if (!base::method_included(m))
				return false;

			if (m == base::method_t::frequency_count)
				return base::meth_frequency_count.deserialize(data);
			if (m == base::method_t::quantile)
				return base::meth_quantile.deserialize(data);
			if (m == base::method_t::deseq2)
				return base::meth_deseq2.deserialize(data);
//			if (m == base::method_t::deseq2_streaming)
//				return base::meth_deseq2_streaming.deserialize(data);

			assert(false);
			
			return false;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void norm_entry(typename normalization_base<ENTRY_T, VALUE_T>::method_t method, Iter1 first1, Iter1 last1, Iter2 first2)
		{
			if (!base::method_included(method))
				return;

			if (method == base::method_t::frequency_count)
				base::meth_frequency_count.norm_entry(first1, last1, first2);
			if (method == base::method_t::quantile)
				base::meth_quantile.norm_entry(first1, last1, first2);
			if (method == base::method_t::deseq2)
				base::meth_deseq2.norm_entry_basic(first1, last1, first2);
			if (method == base::method_t::deseq2_fc)
				base::meth_deseq2.norm_entry_fc(first1, last1, first2);
//			if (method == base::method_t::deseq2_streaming)
//				base::meth_deseq2_streaming.norm_entry(first1, last1, first2);
		}

		// *************************************************************************************
		template<typename Cont1, typename Cont2>
		void norm_entry(typename normalization_base<ENTRY_T, VALUE_T>::method_t method, const Cont1& cont1, Cont2& cont2)
		{
			cont2.resize(cont1.size());
			norm_entry(method, cont1.begin(), cont1.end(), cont2.begin());
		}
	};
}

#endif