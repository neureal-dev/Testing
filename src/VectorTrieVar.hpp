#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace mtl {

template <typename T, size_t W, size_t N>
class PersistentVectorVar {
public:
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;
    using reference_type = typename std::reference_wrapper<T>::type;
    static constexpr size_t getArraySize()
    {
        return {0x1 << W};
    }

    //class Leafs {
    //public:
    // /   std::array<, 0x1 << W> nodes;
    //}/;

    //using Node = leaf_type;

    class Nodes {
    public:
        using node_type = std::array<std::unique_ptr<Nodes>, getArraySize()>;
        using leaf_type = std::array<std::unique_ptr<value_type>, getArraySize()>;
        std::variant<node_type, leaf_type> nodes;
    };

    using leaf_type = Nodes;
    using storage_type = std::array<std::unique_ptr<Nodes>, getArraySize()>;

    template <class... Args>
    reference_type emplace_back(uint32_t id, Args&&... args)
    {
        leaf_type* leaf = leafs_[getReferenceId(id, N - 1)].get();
        if (!leaf) {
            leafs_[getReferenceId(id, N - 1)] = std::make_unique<Nodes>();
            leaf = leafs_[getReferenceId(id, N - 1)].get();
        }

        for (size_t l = N - 2; l > 0; --l) {
            size_t ref_id = getReferenceId(id, l);
            //std::cout << "for" << l << " " << ref_id << " " << id << std::endl;
            //if (!std::holds_alternative<std::unique_ptr<Nodes>>(leaf->nodes[ref_id])) {
            auto& n = std::get<Nodes::node_type>(leaf->nodes)[ref_id];
            if (!n) {
                //std::cout  << "create"<< l << " " << ref_id << " " << id << std::endl;
                n = std::make_unique<Nodes>();
                if (l == 1) {
                    n->nodes = Nodes::leaf_type{};
                }
            }
            leaf = std::get<Nodes::node_type>(leaf->nodes)[ref_id].get();
            //std::cout << l << " " << leaf << std::endl;
        }
        constexpr size_t xbits = ((getArraySize()) - 1);

        std::get<Nodes::leaf_type>(leaf->nodes)[id & xbits] = std::make_unique<T>(std::forward<Args>(args)...);
        return *std::get<Nodes::leaf_type>(leaf->nodes)[id & xbits].get();
    }

    void push_back(value_type&& value)
    {
        emplace_back(value.id_, std::move(value));
    }

    void push_back(const value_type& value)
    {
        emplace_back(value.id_, value);
    }

    [[nodiscard]] pointer_type getNode(uint32_t id) const
    {
        const leaf_type* node = leafs_[getReferenceId(id, N - 1)].get();

        for (size_t l = N - 2; node && l > 0; --l) {
            node = std::get<Nodes::node_type>(node->nodes)[getReferenceId(id, l)].get();
        }

        constexpr size_t xbits = (getArraySize() - 1);
        return node ? std::get<Nodes::leaf_type>(node->nodes)[id & xbits].get() : nullptr;
    }

    void clear() noexcept
    {
        for (auto& a : leafs_) {
            if (a) {
                a.reset();
            }
        }
    }

    void shrink_to_fit() {}

public:
    [[nodiscard]] static constexpr size_t getReferenceId(uint32_t id, size_t n)
    {
        uint64_t xoffset = W * n;
        constexpr uint64_t xbits = 
            getArraySize() - 1;

        return (id >> xoffset) & xbits;
    }

    storage_type leafs_;
};

} // namespace mtl
