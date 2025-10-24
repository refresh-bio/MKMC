#ifndef _STATISTICS_TEST_H
#define _STATISTICS_TEST_H

#include <algorithm>
#include <numeric>
#include <cmath>

#define STATS_ENABLE_STDVEC_WRAPPERS
#include "stats.hpp"

namespace refresh
{
	namespace stats_details
	{
		struct mean_sd_t
		{
			double avg;
			double sd;
			size_t n;
		};

		inline double pow2(double x)
		{
			return x * x;
		}

		template<typename T>
		double pow2(T x)
		{
			return pow2(static_cast<double>(x));
		}

		inline double pow3(double x)
		{
			return x * x * x;
		}

		template<typename T>
		double pow3(T x)
		{
			return pow3(static_cast<double>(x));
		}

		template<typename Iter>
		double mean(Iter first, size_t n)
		{
			if (n == 0)
				return 0;

			return std::accumulate(first, first + n, (double)0.0) / (double)n;
		}

		template<typename Iter>
		double std_dev(double x, Iter first, size_t n)
		{
			if (n < 2)
				return 0;

			double res = 0;

			for (size_t i = 0; i < n; ++i)
				res += pow2(*first++ - x);

			return sqrt(res / (double)(n - 1));
		}

		template<typename X_Iter, typename C_Iter>
		bool mean_std_dev(X_Iter X_first, C_Iter C_first, size_t n_X, std::vector<mean_sd_t>& mean_sd)
		{
			size_t n_class = mean_sd.size();

			if (n_X == 0 || n_class == 0)
				return true;

			for (size_t i = 0; i < n_class; ++i)
			{
				mean_sd[i].avg = 0.0;
				mean_sd[i].sd = 0.0;
				mean_sd[i].n = 0;
			}

			// Calculate mean
			auto p_X = X_first;
			auto p_C = C_first;

			for (size_t i = 0; i < n_X; ++i, ++p_X, ++p_C)
			{
				mean_sd[*p_C].avg += *p_X;
				mean_sd[*p_C].n++;
			}

			for (size_t i = 0; i < n_class; ++i)
				mean_sd[i].avg /= (double)mean_sd[i].n;

			// Calculate std.dev.
			p_X = X_first;
			p_C = C_first;

			for (size_t i = 0; i < n_X; ++i, ++p_X, ++p_C)
				mean_sd[*p_C].sd += pow2(*p_X - mean_sd[*p_C].avg);

			for (size_t i = 0; i < n_class; ++i)
				if (mean_sd[i].n > 1)
					mean_sd[i].sd = sqrt(mean_sd[i].sd / (double)(mean_sd[i].n - 1));

			return true;
		}
	}

	// *************************************************************************************
	// 
	// *************************************************************************************
	class statistical_test
	{
	public:
		struct t_test_t
		{
			double statistic;
			double p_value;
			double df;

			bool isfinite() const
			{
				return std::isfinite(statistic) && std::isfinite(p_value) && std::isfinite(df);
			}

			bool isnormal() const
			{
				return std::isnormal(statistic) && std::isnormal(p_value) && std::isnormal(df);
			}

			bool operator==(const t_test_t& other) const
			{
				return statistic == other.statistic && p_value == other.p_value && df == other.df;
			}

			bool isinf() const
			{
				return std::isinf(statistic) || std::isinf(p_value) || std::isinf(df);
			}

			bool is_nan() const
			{
				return std::isnan(statistic) || std::isnan(p_value) || std::isnan(df);
			}
		};

		struct mann_whitney_t
		{
			double statistic_U1;
			double statistic_U2;
			double p_value;

			bool isfinite() const
			{
				return std::isfinite(statistic_U1) && std::isfinite(statistic_U2) && std::isfinite(p_value);
			}

			bool isnormal() const
			{
				return std::isnormal(statistic_U1) && std::isnormal(statistic_U2) && std::isnormal(p_value);
			}

			bool operator==(const mann_whitney_t& other) const
			{
				return statistic_U1 == other.statistic_U1 && statistic_U2 == other.statistic_U2 && p_value == other.p_value;
			}

			bool isinf() const
			{
				return std::isinf(statistic_U1) || std::isinf(statistic_U2) || std::isinf(p_value);
			}

			bool is_nan() const
			{
				return std::isnan(statistic_U1) || std::isnan(statistic_U2) || std::isnan(p_value);
			}
		};

	private:
		std::vector<std::pair<double, int>> mem_mann_whitney;
		std::vector<stats_details::mean_sd_t> mem_mean_sd;

	public:
		statistical_test() = default;

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		t_test_t t_test(Iter1 X_first, Iter1 X_last, Iter2 C_first, bool equal_var = true)
		{
			return t_test_n(X_first, C_first, std::distance(X_first, X_last), equal_var);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		t_test_t t_test_n(Iter1 X_first, Iter2 C_first, size_t n, bool equal_var = true)
		{
			t_test_t ret;

			auto& mean_sd = mem_mean_sd;

			mean_sd.resize(2);

			stats_details::mean_std_dev(X_first, C_first, n, mean_sd);

			if (equal_var)
			{
				if (n < 3 || mean_sd[0].n < 1 || mean_sd[1].n < 1)
				{
					ret.statistic = std::numeric_limits<double>::quiet_NaN();
					ret.p_value = std::numeric_limits<double>::quiet_NaN();
					ret.df = std::numeric_limits<double>::quiet_NaN();
					return ret;
				}

				double sp = sqrt(((mean_sd[0].n - 1) * stats_details::pow2(mean_sd[0].sd) + (mean_sd[1].n - 1) * stats_details::pow2(mean_sd[1].sd)) / (double) (n - 2));

				if (sp == 0)
				{
					ret.statistic = std::numeric_limits<double>::quiet_NaN();
					ret.p_value = std::numeric_limits<double>::quiet_NaN();
					ret.df = std::numeric_limits<double>::quiet_NaN();
					return ret;
				}

				ret.statistic = (mean_sd[0].avg - mean_sd[1].avg) / (sp * sqrt(1.0 / mean_sd[0].n + 1.0 / mean_sd[1].n));
				ret.df = (double) n - 2.0;
			}
			else
			{
				if (mean_sd[0].n < 2 || mean_sd[1].n < 2)
				{
					ret.statistic = std::numeric_limits<double>::quiet_NaN();
					ret.p_value = std::numeric_limits<double>::quiet_NaN();
					ret.df = std::numeric_limits<double>::quiet_NaN();
					return ret;
				}

				double s_delta = sqrt(stats_details::pow2(mean_sd[0].sd) / (double) mean_sd[0].n + stats_details::pow2(mean_sd[1].sd) / (double)mean_sd[1].n);

				if (s_delta == 0)
				{
					ret.statistic = std::numeric_limits<double>::quiet_NaN();
					ret.p_value = std::numeric_limits<double>::quiet_NaN();
					ret.df = std::numeric_limits<double>::quiet_NaN();
					return ret;
				}

				ret.statistic = (mean_sd[0].avg - mean_sd[1].avg) / s_delta;

				double a0 = stats_details::pow2(mean_sd[0].sd) / mean_sd[0].n;
				double a1 = stats_details::pow2(mean_sd[1].sd) / mean_sd[1].n;

				double denom = (stats_details::pow2(a0) / (mean_sd[0].n - 1) + stats_details::pow2(a1) / (mean_sd[1].n - 1));
				if (denom == 0)
				{
					ret.statistic = std::numeric_limits<double>::quiet_NaN();
					ret.p_value = std::numeric_limits<double>::quiet_NaN();
					ret.df = std::numeric_limits<double>::quiet_NaN();
					return ret;
				}

				ret.df = stats_details::pow2(a0 + a1) / denom;
			}

			ret.p_value = std::clamp<double>(2.0 * (1.0 - stats::pt(fabs(ret.statistic), ret.df, false)), 0, 1);

			return ret;
		}

		// *************************************************************************************
		template<typename X_Iter, typename C_Iter>
		mann_whitney_t mann_whitney_U_test(X_Iter X_first, X_Iter X_last, C_Iter C_first)
		{
			return mann_whitney_U_test_n(X_first, C_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename X_Iter, typename C_Iter>
		mann_whitney_t mann_whitney_U_test_n(X_Iter X_first, C_Iter C_first, size_t n)
		{
//			const size_t approx_method_thr = 8;		// as in SciPy
			const size_t approx_method_thr = 0;

			mann_whitney_t ret;

			mem_mann_whitney.clear();
			mem_mann_whitney.resize(n);

			size_t nc[2] = { 0, 0 };

			for (size_t i = 0; i < n; ++i)
			{
				nc[*C_first]++;
				mem_mann_whitney[i] = std::make_pair(*X_first++, *C_first++);
			}

			std::sort(mem_mann_whitney.begin(), mem_mann_whitney.end());

			double R[2] = { 0, 0 };

			double tie_corr = 0;

			for (size_t i = 0; i < n;)
			{
				size_t j;
				for (j = i + 1; j < n; ++j)
					if (mem_mann_whitney[j].first != mem_mann_whitney[i].first)
						break;

				double adder = 1.0 + (double)(i + j - 1) / 2.0;

				for (size_t k = i; k < j; ++k)
					R[mem_mann_whitney[k].second] += adder;

				if (j > i + 1)
				{
					double t = (double) (j - i);
					tie_corr += stats_details::pow3(t) - (t);
				}

				i = j;
			}

			ret.statistic_U1 = R[0] - (double) nc[0] * (nc[0] + 1) / 2.0;
			ret.statistic_U2 = (double)nc[0] * nc[1] - ret.statistic_U1;

			if (n > approx_method_thr)
			{
				double m_u = (double) nc[0] * nc[1] / 2.0;
				double s_u = sqrt(m_u / 6.0 * ((n + 1) - tie_corr / (n * (n - 1))));

				double cont_term = 0.5;

				double z = (std::max(ret.statistic_U1, ret.statistic_U2) - m_u - cont_term) / s_u;
				ret.p_value = std::clamp<double>(2 * (1.0 - stats::pnorm(z, 0, 1)), 0, 1);
			}
			else
			{
				// !!! TODO
			}

			return ret;
		}

		// *************************************************************************************
		template<typename X_Iter, typename C_Iter>
		double SNR_test(X_Iter X_first, X_Iter X_last, C_Iter C_first)
		{
			return SNR_test_n(X_first, C_first, std::distance(X_first, X_last));
		}

		// *************************************************************************************
		template<typename X_Iter, typename C_Iter>
		double SNR_test_n(X_Iter X_first, C_Iter C_first, size_t n)
		{
			auto& mean_sd = mem_mean_sd;

			mean_sd.resize(2);

			stats_details::mean_std_dev(X_first, C_first, n, mean_sd);

			if (mean_sd[0].sd == 0 || mean_sd[1].sd == 0)
				return 0;

			return std::fabs(mean_sd[0].avg - mean_sd[1].avg) / (mean_sd[0].sd + mean_sd[1].sd);
		}
	};

	// *************************************************************************************
	// 
	// *************************************************************************************
	class p_val_correction
	{
		std::vector<std::pair<double, size_t>> pv_orig;
		std::vector<double> H_n;

		template<typename Iter>
		void load_and_sort(Iter first, size_t n)
		{
			pv_orig.resize(n);

			for (size_t i = 0; i < n; ++i)
			{
				pv_orig[i].first = (double)*first++;
				pv_orig[i].second = i;
			}

			std::sort(pv_orig.begin(), pv_orig.end());
		}

		void build_Hn(size_t n)
		{
			// H_0 = 0
			// H_1 = 1
			// H_2 = 1 + 1/2

			if (n < H_n.size())
				return;

			size_t n_ready = std::max<size_t>(1, H_n.size());
			H_n.resize(n + 1);

			for (size_t i = n_ready; i <= n; ++i)
				H_n[i] = H_n[i - 1] + 1.0 / i;
		}

		template<typename Iter1, typename Iter2>
		void benjamini_base_n(Iter1 S_first, Iter2 D_first, size_t n, double cm)
		{
			load_and_sort(S_first, n);

			double min_pv = 1.0;
			double mult = n * cm;

			for (size_t i = n; i > 0; --i)
			{
				double pv = pv_orig[i - 1].first * mult / i;
				min_pv = std::min(min_pv, pv);
				*(D_first + pv_orig[i - 1].second) = (decltype(*D_first))min_pv;
			}
		}

	public:
		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static void bonferroni(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			bonferroni_n(S_first, D_first, std::distance(S_first, S_last));
		}
			
		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		static void bonferroni_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			for (size_t i = 0; i < n; ++i)
				*D_first++ = (decltype(*D_first))std::clamp<double>(*S_first++ * (double)n, 0, 1);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void holm_bonferroni(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			holm_bonferroni_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void holm_bonferroni_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			load_and_sort(S_first, n);

			double max_pv = 0.0;
			for (size_t i = 0; i < n; ++i)
			{
				double pv = pv_orig[i].first * (double)(n - i);
				max_pv = std::min(std::max(max_pv, pv), 1.0);
				*(D_first + pv_orig[i].second) = (decltype(*D_first))max_pv;
			}
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_hochberg(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			benjamini_hochberg_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_hochberg_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			benjamini_base_n(S_first, D_first, n, 1.0);
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_yekutieli(Iter1 S_first, Iter1 S_last, Iter2 D_first)
		{
			benjamini_yekutieli_n(S_first, D_first, std::distance(S_first, S_last));
		}

		// *************************************************************************************
		template<typename Iter1, typename Iter2>
		void benjamini_yekutieli_n(Iter1 S_first, Iter2 D_first, size_t n)
		{
			build_Hn(n);

			benjamini_base_n(S_first, D_first, n, H_n[n]);
		}
	};
}

#endif
