#pragma once
/// DCDWriter.hpp — minimal, high-throughput C++17 DCD trajectory writer.
///
/// CHARMM-format DCD (compatible with VMD, MDTraj, GROMACS, etc.)
/// POSIX I/O, Linux only.  Single header, no dependencies.
///
/// Coordinate contract  (matches AlignedVecD layout):
///   pos[3*j + 0] = x_j,  pos[3*j + 1] = y_j,  pos[3*j + 2] = z_j
///
/// Compile with -O3 -march=native for AVX-2/512 gather+cvt in deinterleave().

#include <fcntl.h>
#include <unistd.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>

namespace dcd {

/// Simulation box parameters.  Default: 1 Å cube, orthogonal.
struct Box {
    double sideA{1.0};
    double sideB{1.0};
    double sideC{1.0};
    double angleAlpha{90.0}; ///< angle between B and C (degrees)
    double angleBeta{90.0};  ///< angle between A and C (degrees)
    double angleGamma{90.0}; ///< angle between A and B (degrees)
};

/// Streaming DCD writer.
///
/// Open once, call append() per frame, let the destructor (or close()) patch
/// NFRAMES in the header.  Non-copyable, movable.
class Writer {
    public:
    /// @param path      Output path.
    /// @param numAtoms  Atoms per frame.
    /// @param withBox   Write CHARMM periodic-box EXTRA_BLOCK (default: true).
    Writer() noexcept = default;
    ~Writer() noexcept;

    Writer(const Writer&) = delete;
    auto operator=(const Writer&) -> Writer& = delete;

    Writer(Writer&&) noexcept;
    auto operator=(Writer&&) noexcept -> Writer&;

    /// Append one trajectory frame.
    ///
    /// @tparam Vec  Anything with .data() → const double* and size ≥ 3*numAtoms.
    ///              The interleaved layout x0 y0 z0 … xN yN zN is required.
    template <typename Vec>
    auto append(const Vec& coords /* , const Box& box */) -> void {
        assert(fileFd_ >= 0 && "dcd::Writer: not initialized");
        // fillBox(box);
        deinterleaveCoords(coords.data());
        writeAll(frameBuf_.data(), static_cast<std::size_t>(frameBytes_));
        ++frameCount_;
    }

    auto initialize(const std::string& path, int numAtoms, bool withBox = true) -> void;

    /// Flush, patch NFRAMES in the header, and close.  Idempotent.
    auto close() noexcept -> void;

    /// Number of frames written so far.
    [[nodiscard]] auto framesWritten() const noexcept -> int;

    private:
    //  DCD format constants

    /// Half-pi, used for the CHARMM cosine convention.
    static constexpr double kHalfPi{1.57079632679489661922};

    /// Pretend to be CHARMM v24 (non-zero value signals CHARMM format).
    static constexpr std::int32_t kCharmmVersion{24};

    /// Total header size for our fixed layout (CHARMM, 2 title lines, 32-bit).
    static constexpr int kHeaderBytes{276};

    /// File offset of the NFRAMES field (right after record-start + "CORD").
    static constexpr int kNframesOffset{8};

    /// File offset of the NSTEP field (= NFRAMES for istart=0, nsavc=1).
    static constexpr int kNstepOffset{20};

    //  State

    int fileFd_{-1};
    int numAtoms_{0};
    bool withBox_{false};
    int frameCount_{0};

    /// Single reusable output buffer.  Fortran record markers are pre-filled
    /// once in buildFrameTemplate(); only coordinate + box data change per frame.
    std::vector<char> frameBuf_;

    /// Byte offsets of payload regions inside frameBuf_ (past their markers).
    int xOffset_{0};
    int yOffset_{0};
    int zOffset_{0};
    int ucOffset_{0};
    int frameBytes_{0};

    //  Setup

    /// Pre-fill Fortran record markers in frameBuf_ and record data offsets.
    ///
    /// Frame layout  (all counts in bytes):
    ///   [opt]  4 | 48 (6 × double box) | 4       ← CHARMM EXTRA_BLOCK
    ///          4 | numAtoms×4 (float X) | 4
    ///          4 | numAtoms×4 (float Y) | 4
    ///          4 | numAtoms×4 (float Z) | 4
    auto buildFrameTemplate() -> void;

    //  Hot path

    /// Convert box angles to CHARMM cosine convention and memcpy into frameBuf_.
    auto fillBox(const Box& box) noexcept -> void;

    /// Deinterleave x0y0z0…xNyNzN (double) into separate X/Y/Z planes (float)
    /// directly inside frameBuf_.  The stride-3 gather + narrowing cast
    /// vectorises to AVX-2/512 with -O3 -march=native.
    auto deinterleaveCoords(const double* __restrict__ src) noexcept -> void;

    //  I/O

    /// Retry wrapper around POSIX write() — handles EINTR and partial writes.
    auto writeAll(const void* buf, std::size_t size) -> void;

    /// Build the 276-byte DCD header into a stack-allocated array.
    [[nodiscard]] auto buildHeader() const -> std::array<char, kHeaderBytes>;

    /// Patch NFRAMES + NSTEP via pwrite(), then close() the descriptor.
    auto finalize() noexcept -> void;
};

} // namespace dcd