#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <vector>

#include "SearchStorageFixture.hpp"
//#include "VectorTrie.hpp"

namespace {

//#include <algorithm>
//#include <iostream>
//#include <sstream>
//#include <string>


template <typename T>
class Skip_list {
public:
    using key_type = uint32_t;
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;

    Skip_list() : probability(0.5),
                  maxLevel(16)
    {
        int headKey = std::numeric_limits<int>::min();
        head = new Node(headKey, T(0), maxLevel);

        int nilKey = std::numeric_limits<int>::max();
        NIL = new Node(nilKey, T(0xFFFFFFFF), maxLevel);

        std::fill(head->forward.begin(), head->forward.end(), NIL);
    }

    ~Skip_list()
    {
        auto node = head;
        while (node->forward[0]) {
            auto tmp = node;
            node = node->forward[0];
            delete tmp;
        }
        delete node;
    }

    // non-modifying member functions

    /*
    It prints the key, value, level
    of each node of the skip list.

    Prints two nodes per line.
    */
    void print() const
    {
        Node* list = head->forward[0];
        int lineLenght = 0;

        std::cout << "{";

        while (list != NIL) {
            std::cout << "value: " << list->value
                      << ", key: " << list->key
                      << ", level: " << nodeLevel(list);

            list = list->forward[0];

            if (list != NIL) std::cout << " : ";

            if (++lineLenght % 2 == 0) std::cout << "\n";
        }
        std::cout << "}\n";
    }

    /*  
    It searches the skip list and
    returns the element corresponding
    to the searchKey; otherwise it returns
    failure, in the form of null pointer.
    */
    pointer_type find(key_type searchKey) const
    {

        pointer_type res{};
        if (auto x = lower_bound(searchKey)) {
            if (x->key == searchKey && x != NIL) {
                res = &(x->value);
            }
        }
        return res;
    }

    // modifying member functions

    /*
    It searches the skip list for elements
    with seachKey, if there is an element
    with that key its value is reassigned to the
    newValue, otherwise it creates and splices
    a new node, of random level.
    */
    void insert(key_type searchKey, const value_type& newValue)
    {
        auto preds = predecessors(searchKey);

        { //reassign value if node exists and return
            auto next = preds[0]->forward[0];
            if (next->key == searchKey && next != NIL) {
                next->value = newValue;
                return;
            }
        }

        // create new node
        const int newNodeLevel = randomLevel();
        auto newNodePtr = makeNode(searchKey, newValue, newNodeLevel);

        // connect pointers of predecessors and new node to respective successors
        for (int i = 0; i < newNodeLevel; ++i) {
            newNodePtr->forward[i] = preds[i]->forward[i];
            preds[i]->forward[i] = newNodePtr;
        }
    }

    /*
    It deletes the element containing
    searchKey, if it exists.
    */
    void erase(key_type searchKey)
    {
        auto preds = predecessors(searchKey);

        //check if the node exists
        auto node = preds[0]->forward[0];
        if (node->key != searchKey || node == NIL) {
            return;
        }

        // update pointers and delete node
        for (size_t i = 0; i < nodeLevel(node); ++i) {
            preds[i]->forward[i] = node->forward[i];
        }
        delete node;
    }
    void clear()
    {
    }

private:
    struct Node {
        key_type key;
        value_type value;

        // pointers to successor nodes
        std::vector<Node*> forward;

        Node(key_type k, const value_type& v, int level) : key(k), value(v), forward(level, nullptr)
        {
        }
    };

    // Generates node levels in the range [1, maxLevel).
    int randomLevel() const
    {
        int v = 1;
        while (((double)std::rand() / RAND_MAX) < probability && v < maxLevel) {
            v++;
        }
        return v;
    }

    //Returns number of incoming and outgoing pointers
    //    static int nodeLevel(const Node* v);
    static size_t nodeLevel(const Node* v)
    {
        return v->forward.size();
    }

    //creates a node on the heap and returns a pointer to it.
    static Node* makeNode(key_type key, const value_type& val, int level)
    {
        return new Node(key, val, level);
    }

    // Returns the first node for which node->key < searchKey is false
    Node* lower_bound(key_type searchKey) const
    {
        Node* x = head;

        for (size_t i = nodeLevel(head); i-- > 0;) {
            while (x->forward[i]->key < searchKey) {
                x = x->forward[i];
            }
        }
        return x->forward[0];
    }

    /*
    * Returns a collection of Pointers to Nodes
    * result[i] hold the last node of level i+1 for which result[i]->key < searchKey is true
    */
    std::vector<Node*> predecessors(key_type searchKey) const
    {
        std::vector<Node*> result(nodeLevel(head), nullptr);
        Node* x = head;

        for (size_t i = nodeLevel(head); i-- > 0;) {
            while (x->forward[i]->key < searchKey) {
                x = x->forward[i];
            }
            result[i] = x;
        }
        return result;
    }

    // data members
    const float probability;
    const int maxLevel;
    Node* head; // pointer to first node
    Node* NIL; // last node
};

//==============================================================================
// Class Skip_list member implementations

//###### private member functions ######

/*

//==================================================
int main() {

    // 1.Initialize an empty Skip_list object
    Skip_list s;

    // 2. insert()
    for (int i = 0; i < 50; ++i) {
        std::stringstream ss;
        ss << i;

        s.insert(i, ss.str());
    }

    // 2a. print()
    s.print();
    std::cout << std::endl;

    // 3. find()        
    auto f = s.find(10);
    if (f) std::cout << "Node found!\nvalue: " << f << '\n';
    else std::cout << "Node NOT found!\n";

    // 4. insert() - reassign
    s.insert(40, "TEST");

    // 4a. print()
    s.print();
    std::cout << std::endl;

    // 5. erase()
    s.erase(40);

    // 5a. print();
    s.print();
    std::cout << std::endl;

    std::cout << "\nDone!\n";
    getchar();
}
*/
template <typename T>
struct bmmap {
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;

    void push_back(const value_type& value)
    {
        map.insert(value.id_, value);
    }

    pointer_type getNode(uint32_t id)
    {
        return map.find(id);
    }

    void clear()
    {
        map.clear();
    }

    void shrink_to_fit()
    {
        //
    }

    Skip_list<value_type> map;
};

template <class T>
constexpr void ignore(const T&) {}

using SearchOptimization = SearchStorageFixture<bmmap<A>>;

BENCHMARK_DEFINE_F(SearchOptimization, SList)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto rec = records_.getNode(rid);
        if (rec) {
            benchmark::DoNotOptimize(sum += rec->id_);
        } else {
            std::cout << "record not found error" << std::endl;
        }
    }
}
/*
BENCHMARK_REGISTER_F(SearchOptimization, SList)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(1);
*/
} // namespace
