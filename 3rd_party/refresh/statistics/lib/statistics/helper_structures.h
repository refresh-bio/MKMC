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

#if 0		// Old version
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
#else
		template<U64OrSizeT ENTRY_T>
		class compact_histogram<ENTRY_T>
		{
			static constexpr size_t init_plain_size = 32;
			static constexpr size_t init_buffer_size = 32;
			//	static constexpr size_t init_buffer_size = 1024;

			std::vector<size_t> plain;
			std::vector<std::pair<ENTRY_T, size_t>> data;
			std::vector<std::pair<ENTRY_T, size_t>> buffer;

			void compact_buffer()
			{
				if(buffer.size() <= 1)
					return;

				std::sort(buffer.begin(), buffer.end(), [](const auto& a, const auto& b) {return a.first < b.first; });
				size_t i = 0, j = 1;

				while(j < buffer.size())
				{
					if (buffer[i].first == buffer[j].first)
						buffer[i].second += buffer[j].second;
					else
						buffer[++i] = buffer[j];
					++j;
				}

				buffer.resize(i + 1);
			}

			void merge_buffer()
			{
				if (buffer.empty())
					return;

				compact_buffer();

				int64_t i = int64_t(data.size()) - 1;
				int64_t j = int64_t(buffer.size()) - 1;
				int64_t k = int64_t(data.size() + buffer.size()) - 1;

				data.resize(data.size() + buffer.size());

				while (i >= 0 && j >= 0)
				{
					if (data[i].first > buffer[j].first)
						data[k--] = data[i--];
					else
						data[k--] = buffer[j--];
				}

				while (i >= 0)
					data[k--] = data[i--];
				while (j >= 0)
					data[k--] = buffer[j--];

				buffer.clear();
				buffer.reserve(std::max<size_t>(init_buffer_size, (data.size() + plain.size()) / 4));
			}

			void update_plain()
			{
				size_t new_plain_size = plain.size();
				size_t i_data;

				for (i_data = 0; i_data < data.size(); ++i_data)
				{
					auto current_space = plain.size() * 8 + (i_data + 1) * 16;
					auto new_space = (data[i_data].first + 1) * 8;

					if (new_space <= current_space)
						new_plain_size = data[i_data].first + 1;
					else
						break;
				}

				if(new_plain_size == plain.size())
					return;

				plain.resize(new_plain_size, 0);

				for (size_t i = 0; i < i_data; ++i)
					plain[data[i].first] = data[i].second;

				data.erase(data.begin(), data.begin() + i_data);
			}

			void _add(ENTRY_T x, const size_t cnt)
			{
				if (x < plain.size())
				{
					plain[x] += cnt;

					return;
				}

				auto p = std::lower_bound(data.begin(), data.end(), x, [](const auto& entry, ENTRY_T x) {return entry.first < x; });

				if (p != data.end() && p->first == x)
					p->second += cnt;
				else
				{
					buffer.emplace_back(x, cnt);
					if (buffer.size() == buffer.capacity())
					{
						merge_buffer();
						update_plain();
					}
				}
			}

		public:
			compact_histogram()
			{
				plain.resize(init_plain_size);
				buffer.reserve(init_buffer_size);
			}

			void add(ENTRY_T x)
			{
				_add(x, 1);
			}

			void get_histogram(std::vector<std::pair<ENTRY_T, size_t>>& hist)
			{
				merge_buffer();

				hist.clear();

				for(size_t i = 0; i < plain.size(); ++i)
					if (plain[i])
						hist.emplace_back(i, plain[i]);

				hist.reserve(hist.size() + data.size());
				hist.insert(hist.end(), data.begin(), data.end());
			}

			void clear()
			{
				plain.clear();
				plain.resize(init_plain_size);
				data.clear();
				buffer.clear();
				buffer.reserve(init_buffer_size);
			}

			void merge(const compact_histogram<ENTRY_T>& ch)
			{
				for(size_t i = 0; i < ch.plain.size(); ++i)
					if (ch.plain[i] != 0)
						_add(i, ch.plain[i]);

				for (const auto& p : ch.data)
					_add(p.first, p.second);

				for(const auto &p : ch.buffer)
					_add(p.first, p.second);

				merge_buffer();
				update_plain();
			}
		};
#endif

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