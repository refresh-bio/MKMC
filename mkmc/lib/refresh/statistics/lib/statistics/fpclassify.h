#ifndef _FPCLASSIFY_H
#define _FPCLASSIFY_H

#include <concepts>
#include <iterator>

namespace refresh
{
	namespace classify
	{
		template<std::floating_point T>
		bool isnan(const T& x) {
			return std::isnan(x);
		}

		template<typename T>
			requires (!std::floating_point<T>)
		bool isnan(const T& x) {
			return x.isnan();
		}

		template<std::floating_point T>
		bool isfinite(const T& x) {
			return std::isfinite(x);
		}

		template<typename T>
			requires (!std::floating_point<T>)
		bool isfinite(const T& x) {
			return x.isfinite();
		}

		template<std::floating_point T>
		bool isinf(const T& x) {
			return std::isinf(x);
		}

		template<typename T>
			requires (!std::floating_point<T>)
		bool isinf(const T& x) {
			return x.isinf();
		}

		template<std::floating_point T>
		bool isnormal(const T& x) {
			return std::isnormal(x);
		}

		template<typename T>
			requires (!std::floating_point<T>)
		bool isnormal(const T& x) {
			return x.isnormal();
		}


		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool any_nan(It first, It last) {
			for (; first != last; ++first)
				if (std::isnan(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool any_nan(It first, It last) {
			for (; first != last; ++first)
				if (refresh::classify::isnan(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool any_inf(It first, It last) {
			for (; first != last; ++first)
				if (std::isinf(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool any_inf(It first, It last) {
			for (; first != last; ++first)
				if (refresh::classify::isinf(*first)) return true;
			return false;
		}
		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool any_finite(It first, It last) {
			for (; first != last; ++first)
				if (std::isfinite(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool any_finite(It first, It last) {
			for (; first != last; ++first)
				if (refresh::classify::isfinite(*first)) return true;
			return false;
		}
		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool any_normal(It first, It last) {
			for (; first != last; ++first)
				if (std::isnormal(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool any_normal(It first, It last) {
			for (; first != last; ++first)
				if (refresh::classify::isnormal(*first)) return true;
			return false;
		}

		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool all_nan(It first, It last) {
			for (; first != last; ++first)
				if (!std::isnan(*first)) return false;
			return true;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool all_nan(It first, It last) {
			for (; first != last; ++first)
				if (!refresh::classify::isnan(*first)) return false;
			return true;
		}

		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool all_inf(It first, It last) {
			for (; first != last; ++first)
				if (!std::isinf(*first)) return false;
			return true;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool all_inf(It first, It last) {
			for (; first != last; ++first)
				if (!refresh::classify::isinf(*first)) return false;
			return true;
		}
		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool all_finite(It first, It last) {
			for (; first != last; ++first)
				if (!std::isfinite(*first)) return false;
			return true;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool all_finite(It first, It last) {
			for (; first != last; ++first)
				if (!refresh::classify::isfinite(*first)) return false;
			return true;
		}
		template<std::input_iterator It>
			requires std::floating_point<std::iter_value_t<It>>
		bool all_normal(It first, It last) {
			for (; first != last; ++first)
				if (!std::isnormal(*first)) return false;
			return true;
		}

		template<std::input_iterator It>
			requires (!std::floating_point<std::iter_value_t<It>>)
		bool all_normal(It first, It last) {
			for (; first != last; ++first)
				if (!refresh::classify::isnormal(*first)) return false;
			return true;
		}


	}
}

#endif
