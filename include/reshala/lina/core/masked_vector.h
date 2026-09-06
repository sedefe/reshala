#pragma once

#include <vector>

#include "reshala/types.h"

namespace reshala {

class BitMask {
    std::vector<uint64_t> data;
    Index size_;

   public:
    explicit BitMask(Index size = 0) : data((size + 63) / 64, 0), size_(size) {}

    void Resize(Index new_size) {
        size_ = new_size;
        data.resize((size_ + 63) / 64);
    }

    void Clear() { std::fill(data.begin(), data.end(), 0); }

    bool Get(Index pos) const { return (data[pos / 64] >> (pos % 64)) & 1; }

    void Set(Index pos) { data[pos / 64] |= (1ULL << (pos % 64)); }
    void Set(Index pos, bool val) {
        uint64_t mask = 1ULL << (pos % 64);
        data[pos / 64] = (data[pos / 64] & ~mask) | (val * mask);
    }
    void Push(bool val) {
        Index new_size = size_ + 1;
        size_t old_words = data.size();
        size_t new_words = (new_size + 63) / 64;

        if (new_words > old_words) {
            data.push_back(0);
        }
        Set(size_, val);

        size_ = new_size;
    }
};

class MaskedVector {
   public:
    MaskedVector(Index n) : mask_(n) { values_.reserve(n); }

    inline bool Empty() const { return values_.empty(); }

    inline void Clear() {
        mask_.Clear();
        values_.clear();
    }

    inline Index GetNValues() const { return values_.size(); }

    inline void SetMask(const BitMask& mask) { mask_ = mask; }

    inline bool Get(Index i) const { return mask_.Get(i); }

    inline bool Add(Index i) {
        if (!mask_.Get(i)) {
            mask_.Set(i);
            values_.push_back(i);
            return true;
        }
        return false;
    }

    inline const BitMask& GetMask() const { return mask_; }
    inline const std::vector<Index>& GetValues() const { return values_; }

   private:
    BitMask mask_;
    std::vector<Index> values_;
};

}  // namespace reshala
