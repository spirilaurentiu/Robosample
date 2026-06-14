#pragma once

// ============================================================================
//  MemoryArena — one aligned, contiguous memory slab with a bump pointer.
//
//  Purpose: give each World (live RobotState) and each Replica (compact
//  RobotState) ONE cache-line-aligned block, internally partitioned into SoA
//  arrays. No per-array malloc, no vector growth, no pointer chasing between
//  arrays at run time: every array's offset is fixed at build time.
//
//  OpenMP-future-compatible: every sub-array is 64-byte aligned, so later
//  `#pragma omp simd` / aligned AVX-512 loads are legal without re-layout.
//  We do NOT use OpenMP here.
//
//  Usage:
//    MemoryArena arena;
//    // Phase 1: reserve (no allocation yet) — record offsets.
//    auto qOff   = arena.reserve<double>(nq);
//    auto posOff = arena.reserve<double>(3 * nAtoms);
//    ...
//    arena.commit();                       // single aligned allocation
//    double* q   = arena.ptr<double>(qOff);
//    double* pos = arena.ptr<double>(posOff);
// ============================================================================

#include <cstddef>
#include <cstdlib>
#include <new>
#include <stdexcept>

class MemoryArena {
    public:
    static constexpr std::size_t kAlign = 64; // cache line / AVX-512

    MemoryArena() = default;
    ~MemoryArena() {
        std::free(base_);
    }

    MemoryArena(const MemoryArena&) = delete;
    MemoryArena& operator=(const MemoryArena&) = delete;

    MemoryArena(MemoryArena&& o) noexcept
        : base_(o.base_)
        , size_(o.size_)
        , cursor_(o.cursor_)
        , committed_(o.committed_) {
        o.base_ = nullptr;
        o.size_ = 0;
        o.cursor_ = 0;
        o.committed_ = false;
    }
    MemoryArena& operator=(MemoryArena&& o) noexcept {
        if (this != &o) {
            std::free(base_);
            base_ = o.base_;
            size_ = o.size_;
            cursor_ = o.cursor_;
            committed_ = o.committed_;
            o.base_ = nullptr;
            o.size_ = 0;
            o.cursor_ = 0;
            o.committed_ = false;
        }
        return *this;
    }

    // A handle is just a byte offset into the slab (stable across the program).
    using Handle = std::size_t;

    // Reserve space for `count` T's, 64-byte aligned. Returns a handle.
    // Call only before commit().
    template <class T>
    Handle reserve(std::size_t count) {
        if (committed_) {
            throw std::logic_error("MemoryArena::reserve after commit");
        }
        const std::size_t off = alignUp(cursor_, kAlign);
        cursor_ = off + count * sizeof(T);
        return off;
    }

    // Allocate the single aligned block. Idempotent guard.
    void commit() {
        if (committed_) {
            return;
        }
        size_ = alignUp(cursor_, kAlign);
        if (size_ == 0) {
            size_ = kAlign;
        }
        base_ = static_cast<std::byte*>(std::aligned_alloc(kAlign, size_));
        if (base_ == nullptr) {
            throw std::bad_alloc();
        }
        committed_ = true;
    }

    // Typed pointer into the committed slab. Aligned to kAlign by construction.
    template <class T>
    [[nodiscard]] T* ptr(Handle h) const noexcept {
        return reinterpret_cast<T*>(base_ + h);
    }

    [[nodiscard]] std::size_t bytes() const noexcept {
        return size_;
    }
    [[nodiscard]] bool committed() const noexcept {
        return committed_;
    }

    private:
    static std::size_t alignUp(std::size_t n, std::size_t a) noexcept {
        return (n + (a - 1)) & ~(a - 1);
    }

    std::byte* base_ = nullptr;
    std::size_t size_ = 0;
    std::size_t cursor_ = 0;
    bool committed_ = false;
};

// ----------------------------------------------------------------------------
//  Span — a non-owning typed view over a slab region. Lets engine code take
//  flat `double* __restrict` loops over arena arrays with no abstraction cost.
// ----------------------------------------------------------------------------
template <class T>
struct SlabSpan {
    T* data = nullptr;
    std::size_t n = 0;
    [[nodiscard]] T& operator[](std::size_t i) const noexcept {
        return data[i];
    }
    [[nodiscard]] T* begin() const noexcept {
        return data;
    }
    [[nodiscard]] T* end() const noexcept {
        return data + n;
    }
    [[nodiscard]] std::size_t size() const noexcept {
        return n;
    }
};