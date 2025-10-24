#ifndef _HELPER_STRUCTURES_H
#define _HELPER_STRUCTURES_H

#include <concepts>
#include <limits>
#include <algorithm>

namespace refresh
{
	// *************************************************************************************
	namespace normalization::details
	{
		// *************************************************************************************
		// 
		// *************************************************************************************
		template<typename ENTRY_T>
		class compact_histogram
		{
			std::map<ENTRY_T, size_t> data;

			void add(ENTRY_T x, size_t inc)
			{
				data[x] += inc;
			}

		public:
			void add(ENTRY_T x)
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
				for (const auto p : ch.data)
					add(p.first, p.second);
			}
		};

		template<typename T>
		concept U64OrSizeT = std::same_as<T, uint64_t> || std::same_as<T, size_t>;

		template<U64OrSizeT ENTRY_T>
		class compact_histogram<ENTRY_T>
		{
			const size_t thr = 32 << 10;
			//			const size_t thr = 8;				// for testing

			std::vector<ENTRY_T> small;
			std::map<ENTRY_T, size_t> big;

			void add(ENTRY_T x, size_t inc)
			{
				if (x >= thr)
					big[x] += inc;
				else
				{
					if (x >= small.size())
						small.resize(std::min(2 * (static_cast<size_t>(x) + 1), thr));
					small[x] += inc;
				}
			}

		public:
			void add(ENTRY_T x)
			{
				if (x >= thr)
					big[x]++;
				else
				{
					if (x >= small.size())
						small.resize(std::min(2 * (static_cast<size_t>(x) + 1), thr));
					small[x]++;
				}
			}

			void get_histogram(std::vector<std::pair<ENTRY_T, size_t>>& hist)
			{
				hist.clear();

				for (size_t i = 0; i < small.size(); ++i)
					if (small[i])
						hist.emplace_back(i, small[i]);

				for (const auto& x : big)
					hist.emplace_back(x.first, x.second);

				//				hist.shrink_to_fit();
			}

			void clear()
			{
				small.clear();
				small.shrink_to_fit();
				big.clear();
			}

			void merge(const compact_histogram<ENTRY_T>& ch)
			{
				for (size_t i = 0; i < ch.small.size(); ++i)
					if (ch.small[i] != 0)
						add(i, ch.small[i]);
				for (const auto& p : ch.big)
					add(p.first, p.second);
			}
		};

		// *************************************************************************************
		// 
		// *************************************************************************************
		class interval_histogram
		{
		public:
			using hist_t = std::vector<std::pair<double, size_t>>;
		
		protected:
			hist_t hist;

		public:
			interval_histogram()
			{
				clear();
			}

			void init_uniform(double min_val, double max_val, size_t no_split_points)
			{
				hist.clear();
				hist.emplace_back(std::numeric_limits<double>::lowest(), 0);
				hist.emplace_back(min_val, 0);

				for (size_t i = 1; i <= no_split_points; ++i)
					hist.emplace_back(min_val + (max_val * i - min_val * i) / (no_split_points + 1), 0);

				hist.emplace_back(max_val, 0);
				hist.emplace_back(std::numeric_limits<double>::max(), 0);
			}

			template<typename Iter>
			void init_list(Iter first, Iter last)
			{
				hist.clear();
				hist.emplace_back(std::numeric_limits<double>::lowest(), 0);

				for (auto p = first; p != last; ++p)
					hist.emplace_back((double) *p, 0);

				hist.emplace_back(std::numeric_limits<double>::max(), 0);
			}

			size_t get_no_intervals() const
			{
				return hist.size() - 1;
			}

			const hist_t& get_histogram() const
			{
				return hist;
			}

			void clear()
			{
				hist.clear();
				hist.emplace_back(std::numeric_limits<double>::lowest(), 0);
				hist.emplace_back(std::numeric_limits<double>::max(), 0);
			}

			void add(double x)
			{
				auto p = std::lower_bound(hist.begin(), hist.end(), x, 
					[](const auto& entry, double x) {return entry.first < x; });

				p->second++;
			}
		};

		// *************************************************************************************
		// 
		// *************************************************************************************
		template<typename ENTRY_T, typename = void>
		class fast_log
		{
		};

		template<typename ENTRY_T>
		class fast_log<ENTRY_T, typename std::enable_if<std::is_integral<ENTRY_T>::value>::type>
		{
			size_t thr;				

			std::vector<double> precalc_logs;

		public:
//			fast_log(size_t thr = 10)		// for debugging
			fast_log(size_t thr = 32 << 10) : thr(thr)
			{
				precalc_logs.resize(thr);
				for (size_t i = 0; i < thr; ++i)
					precalc_logs[i] = std::log((double)i);
			}

			double log(const ENTRY_T x)
			{
				if ((size_t)x < thr)
					return precalc_logs[x];
				return std::log((double)x);
			}
		};

		template<typename ENTRY_T>
		class fast_log<ENTRY_T, typename std::enable_if<std::is_floating_point<ENTRY_T>::value>::type>
		{
		public:
			fast_log(size_t thr)
			{
			}

			double log(const ENTRY_T x)
			{
				return std::log((double)x);
			}
		};
	}
}

#endif