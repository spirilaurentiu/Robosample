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

/**
 * @brief Streaming single-file CHARMM DCD trajectory writer.
 *
 * @par Lifecycle (required order)
 * Default-construct, then initialize() exactly once, then append() once per
 * frame, then close() (or let the destructor call it). The header's NFRAMES /
 * NSTEP fields are back-patched at close(), so the file is only complete after
 * close()/destruction. Non-copyable, movable.
 *
 * @par Caller contract (unit and frame boundary)
 * The writer performs @b no unit conversion and @b no periodic imaging. It
 * stores the coordinates it is handed verbatim (narrowed to 32-bit float) and
 * treats the Box it is handed as already being in the CHARMM on-disk convention
 * (side lengths, degree angles). The CHARMM DCD convention is Angstrom; the
 * caller (workflow/OutputWriter) is responsible for converting engine nm to
 * Angstrom (x10) and for whole-molecule periodic imaging before calling
 * append(). This class owns only its file descriptor and byte layout.
 */
class Writer {
    public:
    Writer() noexcept = default;
    ~Writer() noexcept;

    Writer(const Writer&) = delete;
    auto operator=(const Writer&) -> Writer& = delete;

    Writer(Writer&&) noexcept;
    auto operator=(Writer&&) noexcept -> Writer&;

    /**
     * @brief Append one trajectory frame.
     * @tparam Vec any type with .data() -> const double* and size >= 3*numAtoms.
     * @param[in] coords borrowed interleaved buffer x0 y0 z0 ... xN yN zN, read
     *        only for the duration of the call; values are stored verbatim as
     *        float (caller supplies Angstrom, imaged coordinates).
     * @param[in] box periodic box for this frame in the CHARMM on-disk convention.
     *        Ignored unless the writer was opened withBox=true.
     * @pre initialize() has been called.
     */
    template <typename Vec>
    auto append(const Vec& coords, const Box& box = Box{}) -> void {
        assert(fileFd_ >= 0 && "dcd::Writer: not initialized");
        if (withBox_) {
            fillBox(box);
        }
        deinterleaveCoords(coords.data());
        writeAll(frameBuf_.data(), static_cast<std::size_t>(frameBytes_));
        ++frameCount_;
    }

    /**
     * @brief Open the output file and write the DCD header. Call once, before any
     *        append().
     * @param[in] path output file path; opened for writing (truncating).
     * @param[in] numAtoms atoms per frame; every append() coordinate buffer must
     *        hold at least 3*numAtoms values in interleaved order.
     * @param[in] withBox write the CHARMM periodic-box EXTRA_BLOCK per frame.
     * @throws std::runtime_error if already initialized or the file cannot be opened.
     */
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