#include <algorithm>
#include <array>
#include <random>
#include <variant>

struct A {
    explicit A(uint32_t id)
            : id_(id), id1_{}, id2_{}, id4_{}
    {
    }

    [[nodiscard]] uint32_t GetId() const noexcept { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    std::array<uint64_t, 3> id4_;
};

struct Comp {

    inline bool operator()(const A* s, uint32_t i) const noexcept { return s->GetId() < i; }

    inline bool operator()(uint32_t i, const A* s) const noexcept { return i < s->GetId(); }

    inline bool operator()(const A& s, uint32_t i) const noexcept { return s.GetId() < i; }

    inline bool operator()(uint32_t i, const A& s) const noexcept { return i < s.GetId(); }

    inline bool operator()(const std::unique_ptr<A>& s, uint32_t i) const noexcept { return s->GetId() < i; }

    inline bool operator()(uint32_t i, const std::unique_ptr<A>& s) const noexcept { return i < s->GetId(); }
};


template <typename StorageType>
struct SearchStorageFixture : benchmark::Fixture {

    void SetUp(const ::benchmark::State& state) final
    {
        if (state.thread_index == 0) {

            searches_.resize(state.range(0), 0);
            //records_.reserve(state.range(0));

            //benchmark::DoNotOptimize(records_.data());

            benchmark::DoNotOptimize(searches_.data());

            std::random_device rd{};
            std::mt19937 rng{rd()};

            rng.seed(rd());

            std::variant<std::uniform_real_distribution<double>,
                    std::normal_distribution<double>,
                    std::exponential_distribution<double>>
                    distribution;

            if (state.range(1) == 0) {
                distribution = std::uniform_real_distribution<double>(1.5, 7.5);
            } else if (state.range(1) == 1) {
                distribution = std::normal_distribution<double>(5.0, 3.0);
            } else {
                distribution = std::exponential_distribution<double>(0.5);
            }

            for (int64_t i = 0; i < state.range(0); ++i) {
                double number = std::visit([&](auto& dist) -> double { return dist(rng); }, distribution);
                if (number > 0.0 && number < 10.0) {
                    searches_[i] = static_cast<uint32_t>((number / 10.0) * state.range(0));
                }
            }

            auto tmp = searches_;

            std::partial_sort_copy(std::begin(searches_), std::end(searches_), std::begin(tmp), std::end(tmp));

            tmp.erase(std::unique(std::begin(tmp), std::end(tmp)), std::end(tmp));

            std::transform(
                    std::begin(tmp),
                    std::end(tmp),
                    std::back_inserter(records_),
                    [](auto i) -> typename StorageType::value_type {
                        return typename StorageType::value_type(i);
                    });

            benchmark::ClobberMemory();
        }
    }

    void TearDown(const ::benchmark::State& state) final
    {
        if (state.thread_index == 0) {
            records_.clear();
            searches_.clear();
            records_.shrink_to_fit();
            searches_.shrink_to_fit();
            benchmark::ClobberMemory();
        }
    }

    StorageType records_;
    std::vector<uint32_t> searches_;
};
