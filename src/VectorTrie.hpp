#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <string.h>

namespace mtl {

template <typename T, int W, int N>
class PersistentVector {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    using leaf_type = PersistentVector<T, W, N - 1>;
    //using storage_type = std::array<leaf_type, W>;
    using storage_type = std::array<std::unique_ptr<leaf_type>, 0x1 << W>;

    void addNode(uint32_t id, pointer_type element)
    {
        /*/
        if (node_ptrs == nullptr) {
            node_ptrs = new char[1024 * 1024 * 1024];
            memset(node_ptrs, '\0', 1024 * 1024 * 1024);
            node_ptrss = node_ptrs;
            list.nodes.fill(nullptr);
        }
        /*/
        size_t ref_id = getReferenceId(id);
        if (!leafs_[ref_id]) {
            //            nodes++;
            leafs_[ref_id] = std::make_unique<value_type>();
            //list.nodes[ref_id] = new (node_ptrss) type();
            //node_ptrss += sizeof(type);
            //list.nodes[ref_id]->list.nodes.fill(nullptr);
        }
        //*/
        leafs_[ref_id]->addNode(id, element);
    }

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        size_t ref_id = getReferenceId(id);
        if (!leafs_[ref_id]) {
            leafs_[ref_id] = std::make_unique<leaf_type>();
        }
        return leafs_[ref_id]->emplace_back(id, args...);
    }

    pointer_type getNode(uint32_t id) const
    {
        auto& node = leafs_[getReferenceId(id)];
        return node ? node->getNode(id) : nullptr;
    }

    void clear() noexcept
    {
        for (auto& a : leafs_) {
            if (a) {
                a.reset();
            }
        }
    }

public:
    inline size_t getReferenceId(uint32_t id) const
    {
        constexpr uint64_t xoffset = W * (N - 1);
        constexpr uint64_t xbits = (0x1 << W) - 1;

        return {(id >> xoffset) & xbits};
    }

    storage_type leafs_;
};

template <typename T, int W>
class PersistentVector<T, W, 0> {
public:
    using type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    //using storage_type = std::array<pointer_type, W>;
    using storage_type = std::array<std::unique_ptr<type>, 0x1 << W>;
    //using leaf_type = or<T, W, N - 1>;

    void addNode(uint32_t id, pointer_type element)
    {
        leafs_[getReferenceId(id)] = element;
    }

    pointer_type getNode(uint32_t id) const
    {
        return leafs_[getReferenceId(id)].get();
    }

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        leafs_[getReferenceId(id)] = std::make_unique<T>(id, args...);
        return *leafs_[getReferenceId(id)].get();
    }

private:
    inline size_t getReferenceId(uint32_t id) const
    {
        constexpr size_t xbits = ((0x1 << W) - 1);

        return {(id >> W) & xbits};
    }

    storage_type leafs_;
};

} // namespace mtl
