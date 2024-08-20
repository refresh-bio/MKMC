#ifndef _STATISTICS_UMAP_H
#define _STATISTICS_UMAP_H

#define EIGEN_NO_DEPRECATED_WARNING
#include "umappp/Umap.hpp"

#include <vector>
#include <iterator>
#include <cinttypes>

namespace refresh
{
	template<typename VALUE_T>
	class umap
	{
	public:
		struct params_t
		{
			VALUE_T local_connectivity = umappp::Umap<VALUE_T>::Defaults::local_connectivity;
			VALUE_T bandwidth = umappp::Umap<VALUE_T>::Defaults::bandwidth;
			VALUE_T mix_ratio = umappp::Umap<VALUE_T>::Defaults::mix_ratio;
			VALUE_T spread = umappp::Umap<VALUE_T>::Defaults::spread;
			VALUE_T min_dist = umappp::Umap<VALUE_T>::Defaults::min_dist;
			VALUE_T a = umappp::Umap<VALUE_T>::Defaults::a;
			VALUE_T b = umappp::Umap<VALUE_T>::Defaults::b;
			VALUE_T repulsion_strength = umappp::Umap<VALUE_T>::Defaults::repulsion_strength;
			umappp::InitMethod initialize = umappp::Umap<VALUE_T>::Defaults::initialize;
			int num_epochs = umappp::Umap<VALUE_T>::Defaults::num_epochs;
			VALUE_T learning_rate = umappp::Umap<VALUE_T>::Defaults::learning_rate;
			VALUE_T negative_sample_rate = umappp::Umap<VALUE_T>::Defaults::negative_sample_rate;
			int num_neighbors = umappp::Umap<VALUE_T>::Defaults::num_neighbors;
			uint64_t seed = umappp::Umap<VALUE_T>::Defaults::seed;
			int num_threads = umappp::Umap<VALUE_T>::Defaults::num_threads;
			int parallel_optimization = umappp::Umap<VALUE_T>::Defaults::parallel_optimization;
		};

	private:
		enum class input_mode_t {unknown, feature_oriented, object_oriented};

		using entry_t = std::vector<VALUE_T>;

		std::vector<entry_t> input_data;
		std::vector<VALUE_T> umap_in_data;
		std::vector<VALUE_T> umap_out_data;
		std::vector<entry_t> ret_data;

		input_mode_t input_mode = input_mode_t::unknown;
		size_t input_vector_size{};

		params_t params;

		template<typename VEC>
		void vec_clear(VEC& v)
		{
			v.clear();
			v.shrink_to_fit();
		}

		template<typename Iter>
		bool try_add(input_mode_t requested_mode, Iter first, Iter last)
		{
			if (input_mode == input_mode_t::unknown)
			{
				input_mode = requested_mode;
				input_vector_size = std::distance(first, last);
			}
			else if (input_mode != requested_mode)
				return false;

			if (std::distance(first, last) != input_vector_size)
				return false;

			input_data.emplace_back(first, last);

			return true;
		}

		void clear()
		{
			vec_clear(input_data);
			vec_clear(umap_in_data);
			vec_clear(umap_out_data);
			vec_clear(ret_data);

			input_mode = input_mode_t::unknown;
			input_vector_size = 0;
		}

	public:
		umap() = default;

		void reset()
		{
			clear();
		}

		template<typename Iter>
		bool add_object(Iter first, Iter last)
		{
			return try_add(input_mode_t::object_oriented, first, last);
		}

		template<typename Iter>
		bool add_feature(Iter first, Iter last)
		{
			return try_add(input_mode_t::feature_oriented, first, last);
		}

		bool run(size_t no_dimensions = 2)
		{
			if (input_mode == input_mode_t::unknown)
				return false;

			umappp::Umap um;

			um.set_local_connectivity(params.local_connectivity);
			um.set_bandwidth(params.bandwidth);
			um.set_mix_ratio(params.mix_ratio);
			um.set_spread(params.spread);
			um.set_min_dist(params.min_dist);
			um.set_a(params.a);
			um.set_b(params.b);
			um.set_repulsion_strength(params.repulsion_strength);
			um.set_initialize(params.initialize);
			um.set_num_epochs(params.num_epochs);
			um.set_learning_rate(params.learning_rate);
			um.set_negative_sample_rate(params.negative_sample_rate);
			um.set_num_neighbors(params.num_neighbors);
			um.set_seed(params.seed);
			um.set_num_threads(params.num_threads);
			um.set_parallel_optimization(params.parallel_optimization);

			// Prepare UMAP input data in column-major order
			umap_in_data.resize(input_data.size() * input_data.front().size());

			size_t k = 0;
			size_t no_objects = 0;
			size_t no_features = 0;

			switch (input_mode)
			{
			case input_mode_t::object_oriented:
				no_objects = input_data.size();
				no_features = input_data.front().size();

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						umap_in_data[k++] = input_data[i][j];

				break;
			case input_mode_t::feature_oriented:
				no_objects = input_data.front().size();
				no_features = input_data.size();

				for (size_t i = 0; i < no_objects; ++i)
					for (size_t j = 0; j < no_features; ++j)
						umap_in_data[k++] = input_data[j][i];

				break;
			default:
				assert(0);						// Never should be here
			}

			// Prepare space for UMAP results
			umap_out_data.resize(no_dimensions * no_objects);

			// Run UMAP
			um.run(no_features, no_objects, umap_in_data.data(), no_dimensions, umap_out_data.data());

			// Convert results to vector of vectors
			ret_data.clear();
			ret_data.resize(no_objects);

			k = 0;
			for (size_t i = 0; i < no_objects; ++i)
				for (size_t j = 0; j < no_dimensions; ++j)
					ret_data[i].emplace_back(umap_out_data[k++]);

			return true;
		}

		std::vector<std::vector<VALUE_T>>& result()
		{
			return ret_data;
		}

		params_t get_defaults() const
		{
			return params_t();
		}

		params_t get_params() const
		{
			return params;
		}

		void set_params(params_t& new_params)
		{
			params = new_params;
		}
	};

	template<typename VALUE_T>
	class umap_direct
	{
		using entry_t = std::vector<VALUE_T>;

		size_t no_objects;
		size_t no_features;
		using params_t = typename umap<VALUE_T>::params_t;
		params_t params;
		std::unique_ptr<VALUE_T[]> umap_in_data; //unique_ptr instead of std::vector to avoid memory initialization
		std::vector<VALUE_T> umap_out_data;
		std::vector<entry_t> ret_data;
	public:
		umap_direct(size_t no_objects, size_t no_features) :
			no_objects(no_objects),
			no_features(no_features),
			umap_in_data(std::make_unique_for_overwrite<VALUE_T[]>(no_objects* no_features))
		{

		}

		void add(size_t object_idx, size_t feature_idx, const VALUE_T value)
		{
			umap_in_data[object_idx * no_features + feature_idx] = value;
		}

		void run(size_t no_dimensions = 2)
		{
			umappp::Umap um;

			um.set_local_connectivity(params.local_connectivity);
			um.set_bandwidth(params.bandwidth);
			um.set_mix_ratio(params.mix_ratio);
			um.set_spread(params.spread);
			um.set_min_dist(params.min_dist);
			um.set_a(params.a);
			um.set_b(params.b);
			um.set_repulsion_strength(params.repulsion_strength);
			um.set_initialize(params.initialize);
			um.set_num_epochs(params.num_epochs);
			um.set_learning_rate(params.learning_rate);
			um.set_negative_sample_rate(params.negative_sample_rate);
			um.set_num_neighbors(params.num_neighbors);
			um.set_seed(params.seed);
			um.set_num_threads(params.num_threads);
			um.set_parallel_optimization(params.parallel_optimization);


			// Prepare space for UMAP results
			umap_out_data.resize(no_dimensions * no_objects);

			// Run UMAP
			um.run(no_features, no_objects, umap_in_data.get(), no_dimensions, umap_out_data.data());

			ret_data.clear();
			ret_data.shrink_to_fit();
			ret_data.resize(no_objects);

			size_t k = 0;
			for (size_t i = 0; i < no_objects; ++i)
				for (size_t j = 0; j < no_dimensions; ++j)
					ret_data[i].emplace_back(umap_out_data[k++]);
		}

		std::vector<std::vector<VALUE_T>>& result()
		{
			return ret_data;
		}

		params_t get_defaults() const
		{
			return params_t();
		}

		params_t get_params() const
		{
			return params;
		}

		void set_params(params_t& new_params)
		{
			params = new_params;
		}
	};

}

#endif
