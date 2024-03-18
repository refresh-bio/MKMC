#pragma once

#include <vector>



//when we don't know the total size
//in fact the above could be implemented using this as implementation and just wrapping it with vector of read_pos
template<typename T, typename COMPARATOR>
class BinaryHeapMergeStreams
{
	template<typename U>
	struct heap_desc_t
	{
		U elem;
		uint64_t id;
		heap_desc_t(U elem, uint64_t id) :
			elem(elem), id(id) {
		}
	};

	const COMPARATOR comparator;

	std::vector<heap_desc_t<T>> heap;

	void heap_down();

public:
	template<typename DO_WITH_ELEM_IF_EXISTS>
	BinaryHeapMergeStreams(const std::vector<size_t>& streams_to_merge, const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists, const COMPARATOR& comparator);

	template<typename DO_WITH_ELEM_IF_EXISTS, typename Callback>
	void ProcessElem(const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists, Callback&& callback);

	bool Empty() const {
		return heap.empty();
	}
};



template<typename T, typename COMPARATOR>
inline void BinaryHeapMergeStreams<T, COMPARATOR>::heap_down()
{
	uint64_t parent = 0;
	uint64_t left = 1;
	uint64_t right = 2;

	auto elem = heap[0];

	//has 2 childs
	while (right < heap.size()) {
		left = comparator(heap[right].elem, heap[left].elem) ? left : right; //left is min now

		if (!comparator(elem.elem, heap[left].elem)) {
			heap[parent] = elem;
			return;
		}

		heap[parent] = heap[left];
		parent = left;
		left = parent * 2 + 1;
		right = left + 1;
	}

	if (left < heap.size()) {
		if (comparator(heap[left].elem, elem.elem)) {
			heap[parent] = elem;
			return;
		}
		heap[parent] = heap[left];
		parent = left;
	}

	heap[parent] = elem;
}

template<typename T, typename COMPARATOR>
template<typename DO_WITH_ELEM_IF_EXISTS>
inline BinaryHeapMergeStreams<T, COMPARATOR>::BinaryHeapMergeStreams(const std::vector<size_t>& streams_to_merge, const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists, const COMPARATOR& comparator)
	: comparator(comparator)
{
	heap.reserve(streams_to_merge.size());

	for (size_t i = 0; i < streams_to_merge.size(); ++i) {
		do_with_elem_if_exists(streams_to_merge[i], [&](const T& elem) { heap.emplace_back(elem, streams_to_merge[i]); });
	}

	std::make_heap(heap.begin(), heap.end(), [&](const heap_desc_t<T>& elem1, const heap_desc_t<T>& elem2) -> bool {
		return comparator(elem1.elem, elem2.elem);
		});
}



template<typename T, typename COMPARATOR>
template<typename DO_WITH_ELEM_IF_EXISTS, typename Callback>
inline void BinaryHeapMergeStreams<T, COMPARATOR>::ProcessElem(const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists, Callback&& callback)
{
	auto id = heap[0].id;

	callback(heap[0].elem, id);

	if (do_with_elem_if_exists(id, [&](const T& elem)
		{
			heap[0].elem = elem;
			heap_down();
		}))
	{
	} //do nothing in if's body, because it is done in the callback of do_with_elem_if_exists
	else
	{
		heap[0] = heap.back();
		heap.pop_back();
		if (heap.size() > 1)
			heap_down();
	}
}
