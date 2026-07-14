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

/**
 * @brief One aligned, contiguous memory slab handed out as fixed byte-offset
 *        handles: a two-phase bump allocator for struct-of-arrays layouts.
 *
 * @par Usage protocol (two phases)
 * Reserve every sub-array with reserve() (records offsets, allocates nothing),
 * then commit() once (a single aligned allocation), then obtain typed pointers
 * with ptr(). reserve() after commit() throws; commit() is idempotent.
 *
 * @par Ownership and lifetime
 * The arena @b owns its slab and frees it on destruction or move-assignment.
 * Every pointer from ptr() is @b borrowed: it stays valid only while this arena
 * is alive and un-reassigned. Move transfers the slab and nulls the source;
 * move-assigning a new arena frees the previous slab, invalidating all pointers
 * previously obtained from the destination. Non-copyable.
 *
 * @note Every sub-array begins on a @c kAlign (64-byte) boundary, so aligned
 *       SIMD loads over any reserved array are legal.
 */
class MemoryArena {
    public:
    /** @brief Alignment of the slab and of every reserved sub-array, in bytes. */
    static constexpr std::size_t kAlign = 64; // cache line / AVX-512

    MemoryArena() = default;
    ~MemoryArena() {
        std::free(base_);
    }

    MemoryArena(const MemoryArena&) = delete;
    MemoryArena& operator=(const MemoryArena&) = delete;

    /** @brief Move-construct: take @p o's slab and leave @p o empty. */
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
    /** @brief Move-assign: free this slab, take @p o's, and leave @p o empty.
     *  Invalidates every pointer previously obtained from this arena. */
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

    /** @brief A reserved region's byte offset into the slab; stable for the
     *  arena's lifetime. Resolve to a pointer with ptr(). */
    using Handle = std::size_t;

    /**
     * @brief Reserve space for @p count objects of type @c T, aligned to kAlign.
     * @param[in] count number of elements; the region spans count*sizeof(T) bytes.
     * @return handle (byte offset) to pass to ptr() after commit().
     * @pre commit() has not been called.
     * @throws std::logic_error if called after commit().
     * @note Reserves offsets only; no memory is allocated until commit().
     */
    template <class T>
    Handle reserve(std::size_t count) {
        if (committed_) {
            throw std::logic_error("MemoryArena::reserve after commit");
        }
        const std::size_t off = alignUp(cursor_, kAlign);
        cursor_ = off + count * sizeof(T);
        return off;
    }

    /**
     * @brief Allocate the single aligned slab sized to hold all reserved regions.
     * @post ptr() may be called on any prior handle. A zero-reservation arena
     *       still allocates one kAlign block.
     * @throws std::bad_alloc if the allocation fails.
     * @note Idempotent: a second call is a no-op (does not reallocate).
     */
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

    /**
     * @brief Resolve a handle to a typed, borrowed pointer into the committed slab.
     * @param[in] h a handle returned by reserve().
     * @return pointer to the region, aligned to kAlign; valid only while this
     *         arena is alive and un-reassigned.
     * @pre commit() has been called (otherwise the base is null and the result is
     *      undefined).
     */
    template <class T>
    [[nodiscard]] T* ptr(Handle h) const noexcept {
        return reinterpret_cast<T*>(base_ + h);
    }

    /** @brief Total committed slab size in bytes (0 before commit()). */
    [[nodiscard]] std::size_t bytes() const noexcept {
        return size_;
    }
    /** @brief True once commit() has run. */
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

/**
 * @brief Non-owning typed view over a slab region (pointer + element count).
 * @note Borrows; does not own or free. Lifetime is bounded by the arena the
 *       @c data pointer came from.
 */
template <class T>
struct SlabSpan {
    /** @brief Borrowed pointer to the first element. */
    T* data = nullptr;
    /** @brief Number of elements in the view. */
    std::size_t n = 0;
    /** @brief Unchecked element access. */
    [[nodiscard]] T& operator[](std::size_t i) const noexcept {
        return data[i];
    }
    /** @brief Iterator to the first element. */
    [[nodiscard]] T* begin() const noexcept {
        return data;
    }
    /** @brief Iterator one past the last element. */
    [[nodiscard]] T* end() const noexcept {
        return data + n;
    }
    /** @brief Element count. */
    [[nodiscard]] std::size_t size() const noexcept {
        return n;
    }
};