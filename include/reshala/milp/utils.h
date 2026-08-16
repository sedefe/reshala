#pragma once

#include <atomic>

#include "reshala/lp/dual_simplex.h"

namespace reshala {

enum class Direction { kLeft = 0, kRight = 1 };
inline Index Dir2Index(Direction dir) { return static_cast<Index>(dir); }
inline Direction Index2Dir(Index i) { return static_cast<Direction>(i); }

struct Node {
    Node() {}
    Node(Index l, const Solution& s, const Domain& d, const DsState& st)
        : id(next_id()), level(l), sol(s), domain(d), ds_state(st) {}

    static Index next_id() {
        static std::atomic<Index> counter{0};
        return counter.fetch_add(1, std::memory_order_relaxed);
    }

    Index id;
    Index level;
    Solution sol;
    Domain domain;

    DsState ds_state;
};

inline void ReportNewPrimal(const std::string& src, Scalar y) {
    std::cout << src << " found solution " << FMT(0, 5) << y << "\n";
}

}  // namespace reshala
