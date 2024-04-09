#ifndef _STATISTICS_CORRELATION
#define _STATISTICS_CORRELATION

#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

namespace refresh
{
	// *************************************************************************************
	class correlation
	{
		std::vector<std::pair<double, size_t>> spearman_mapping;
		std::vector<double> spearman_X_ranks, spearman_Y_ranks;

		// *************************************************************************************
		void set_spearman_ranks(std::vector<double>& ranks, size_t from, size_t to)
		{
			double rank = (to + from) / 2.0;

			for (auto i = from; i <= to; ++i)
				ranks[spearman_mapping[i].second] = rank;
		}

		template<typename Iter1>
		void spearman_remap(Iter1 first, size_t n, std::vector<double>& ranks)
		{
			spearman_mapping.resize(n);
			ranks.resize(n);

			auto p = first;

			for (size_t i = 0; i < n; ++i, ++p)
			{
				spearman_mapping[i].first = (double)*p;
				spearman_mapping[i].second = i;
			}

			std::sort(spearman_mapping.begin(), spearman_mapping.end());

			size_t j = 0;

			for (size_t i = 1; i < n; ++i)
				if (spearman_mapping[i].first != spearman_mapping[j].first)
				{
					set_spearman_ranks(ranks, j, i - 1);
					j = i;
				}

			set_spearman_ranks(ranks, j, n - 1);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double kendall_tau_b_imp_naive(Iter1 X_first, Iter2 Y_first, size_t n)
		{
			size_t n_concordant = 0;
			size_t n_discordant = 0;
			size_t n_ties_X = 0;
			size_t n_ties_Y = 0;

			auto p_Xi = X_first;
			auto p_Yi = Y_first;

			for (size_t i = 0; i < n - 1; ++i, ++p_Xi, ++p_Yi)
			{
				auto p_Xj = p_Xi;		++p_Xj;
				auto p_Yj = p_Yi;		++p_Yj;

				for (size_t j = i + 1; j < n; ++j, ++p_Xj, ++p_Yj)
				{
					if (*p_Xi < *p_Xj && *p_Yi < *p_Yj)
						++n_concordant;
					else if (*p_Xi > *p_Xj && *p_Yi > *p_Yj)
						++n_concordant;
					else if (*p_Xi > *p_Xj && *p_Yi < *p_Yj)
						++n_discordant;
					else if (*p_Xi < *p_Xj && *p_Yi > *p_Yj)
						++n_discordant;
					else
					{
						n_ties_X += *p_Xi == *p_Xj;
						n_ties_Y += *p_Yi == *p_Yj;
					}
				}
			}

			size_t n_pairs = n * (n - 1) / 2;

			return ((double)(n_concordant)-(double)n_discordant) / sqrt(((double)n_pairs - (double)n_ties_X) * ((double)n_pairs - (double)n_ties_Y));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double kendall_tau_b_imp_fast(Iter1 X_first, Iter2 Y_first, size_t n)
		{
			using X_t = typename Iter1::value_type;
			using Y_t = typename Iter2::value_type;

			std::vector<std::pair<X_t, Y_t>> v1(n), v2(n);

			auto p_X = X_first;
			auto p_Y = Y_first;

			// Pairing
			for (size_t i = 0; i < n; ++i, ++p_X, ++p_Y)
			{
				v1[i].first = *p_X;
				v1[i].second = *p_Y;
			}

			// Sorting by X
			std::sort(v1.begin(), v1.end());

			size_t n_ties_X = 0;
			size_t n_ties_Y = 0;
			size_t n_ties_both = 0;

			// Fiding no. of X and both ties
			size_t j_X = 0;
			size_t j_both = 0;
			for (size_t i = 1; i < n; ++i)
			{
				if (v1[i] != v1[j_both])
				{
					n_ties_both += (i - j_both) * (i - j_both - 1) / 2;
					j_both = i;
				}

				if (v1[i].first != v1[j_X].first)
				{
					n_ties_X += (i - j_X) * (i - j_X - 1) / 2;
					j_X = i;
				}
			}

			n_ties_both += (n - j_both) * (n - j_both - 1) / 2;
			n_ties_X += (n - j_X) * (n - j_X - 1) / 2;

			size_t n_swaps = 0;

			// Merge sort
			for (size_t list_size = 1; list_size < n; list_size *= 2)
			{
				size_t p1, p2, p3;
				p1 = 0;

				while (p1 < n)
				{
					p2 = p1 + list_size;
					if (p2 >= n)
					{
						std::copy(v1.begin() + p1, v1.end(), v2.begin() + p1);
						break;
					}
					p3 = p2 + list_size;
					if (p3 > n)
						p3 = n;

					size_t i = 0;
					size_t j = 0;
					size_t k = 0;
					size_t n = p2 - p1;
					size_t m = p3 - p2;

					while (i < n && j < m)
					{
						if (v1[p1 + i].second > v1[p2 + j].second)
						{
							n_swaps += n - i;
							v2[p1 + k] = v1[p2 + j];
							++j;
						}
						else
						{
							v2[p1 + k] = v1[p1 + i];
							++i;
						}
						++k;
					}

					// Only one of the copy executions will be made, so no need to increment k
					std::copy(v1.begin() + p1 + i, v1.begin() + p2, v2.begin() + p1 + k);
					std::copy(v1.begin() + p2 + j, v1.begin() + p3, v2.begin() + p1 + k);

					p1 = p3;
				}

				swap(v1, v2);
			}

			// Finding Y ties
			size_t j_Y = 0;
			for (size_t i = 1; i < n; ++i)
			{
				if (v1[i].second != v1[j_Y].second)
				{
					n_ties_Y += (i - j_Y) * (i - j_Y - 1) / 2;
					j_Y = i;
				}
			}

			n_ties_Y += (n - j_Y) * (n - j_Y - 1) / 2;

			double n_pairs = (double)(n * (n - 1) / 2);
			double numerator = n_pairs - n_ties_X - n_ties_Y + n_ties_both - 2 * n_swaps;

			return numerator / sqrt((n_pairs - (double)n_ties_X) * (n_pairs - (double)n_ties_Y));
		}

	public:
		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double pearson(Iter1 X_first, Iter1 X_last, Iter2 Y_first)
		{
			return pearson_n(X_first, Y_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double pearson_n(Iter1 X_first, Iter2 Y_first, size_t n)
		{
			if (n == 0)
				return 0.0;

			double X_avg = 0;
			double Y_avg = 0;

			auto p_X = X_first;
			auto p_Y = Y_first;

			for (size_t i = 0; i < n; ++i, ++p_X, ++p_Y)
			{
				X_avg += *p_X;
				Y_avg += *p_Y;
			}

			X_avg /= n;
			Y_avg /= n;

			p_X = X_first;
			p_Y = Y_first;

			double A = 0;
			double B = 0;
			double C = 0;

			for (size_t i = 0; i < n; ++i, ++p_X, ++p_Y)
			{
				double dx = *p_X - X_avg;
				double dy = *p_Y - Y_avg;

				A += dx * dy;
				B += dx * dx;
				C += dy * dy;
			}

			double D = sqrt(B) * sqrt(C);

			if (D == 0)
				return 1.0;

			return A / D;
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		double spearman(Iter1 X_first, Iter1 X_last, Iter2 Y_first, bool release_memory = false)
		{
			return spearman_n(X_first, Y_first, std::distance(X_first, X_last), release_memory);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		double spearman_n(Iter1 X_first, Iter2 Y_first, size_t n, bool release_memory = false)
		{
			if (n == 0)
				return 0.0;

			spearman_remap(X_first, n, spearman_X_ranks);
			spearman_remap(Y_first, n, spearman_Y_ranks);

			auto r = pearson_n(spearman_X_ranks.begin(), spearman_Y_ranks.begin(), n);

			if (release_memory)
				release_memory_spearman();

			return r;
		}

		// *************************************************************************************
		void release_memory_spearman()
		{
			spearman_mapping.clear();
			spearman_mapping.shrink_to_fit();
			spearman_X_ranks.clear();
			spearman_X_ranks.shrink_to_fit();
			spearman_Y_ranks.clear();
			spearman_Y_ranks.shrink_to_fit();
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double kendall_tau(Iter1 X_first, Iter1 X_last, Iter2 Y_first)
		{
			return kendall_tau_n(X_first, Y_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static double kendall_tau_n(Iter1 X_first, Iter2 Y_first, size_t n)
		{
			const size_t naive_variant_thr = 10;

			if (n == 0)
				return 0;

			if (n < naive_variant_thr)
				return kendall_tau_b_imp_naive(X_first, Y_first, n);
			else
				return kendall_tau_b_imp_fast(X_first, Y_first, n);
		}

	};

}

#endif