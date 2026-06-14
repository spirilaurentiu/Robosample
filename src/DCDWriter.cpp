#include "DCDWriter.hpp"

namespace dcd {

Writer::~Writer() noexcept {
    if (fileFd_ >= 0) {
        finalize();
    }
}

Writer::Writer(Writer&& other) noexcept
    : fileFd_{other.fileFd_}
    , numAtoms_{other.numAtoms_}
    , withBox_{other.withBox_}
    , frameCount_{other.frameCount_}
    , frameBuf_{std::move(other.frameBuf_)}
    , xOffset_{other.xOffset_}
    , yOffset_{other.yOffset_}
    , zOffset_{other.zOffset_}
    , ucOffset_{other.ucOffset_}
    , frameBytes_{other.frameBytes_} {
    other.fileFd_ = -1;
}

auto Writer::operator=(Writer&& other) noexcept -> Writer& {
    if (this != &other) {
        if (fileFd_ >= 0) {
            finalize();
        }
        fileFd_ = other.fileFd_;
        numAtoms_ = other.numAtoms_;
        withBox_ = other.withBox_;
        frameCount_ = other.frameCount_;
        frameBuf_ = std::move(other.frameBuf_);
        xOffset_ = other.xOffset_;
        yOffset_ = other.yOffset_;
        zOffset_ = other.zOffset_;
        ucOffset_ = other.ucOffset_;
        frameBytes_ = other.frameBytes_;
        other.fileFd_ = -1;
    }
    return *this;
}

auto Writer::initialize(const std::string& path, int numAtoms, bool withBox) -> void {
    if (fileFd_ >= 0) {
        throw std::runtime_error{"dcd::Writer: already initialized"};
    }

    this->numAtoms_ = numAtoms;
    this->withBox_ = withBox;

    fileFd_ = ::open(path.c_str(), (O_WRONLY | O_CREAT | O_TRUNC), 0666);
    if (fileFd_ < 0) {
        throw std::runtime_error{"dcd::Writer: cannot open '" + path + '\''};
    }

    buildFrameTemplate();

    const auto hdr = buildHeader();
    writeAll(hdr.data(), static_cast<std::size_t>(kHeaderBytes));
}

auto Writer::close() noexcept -> void {
    if (fileFd_ >= 0) {
        finalize();
        fileFd_ = -1;
    }
}

auto Writer::framesWritten() const noexcept -> int {
    return frameCount_;
}

auto Writer::buildFrameTemplate() -> void {
    const auto coordBytes = static_cast<std::int32_t>(numAtoms_ * 4);
    const int ucBlock = withBox_ ? (4 + 48 + 4) : 0;
    const int coordBlock = (4 + (numAtoms_ * 4) + 4);

    frameBytes_ = ucBlock + (3 * coordBlock);
    frameBuf_.assign(static_cast<std::size_t>(frameBytes_), '\0');

    auto poke32 = [this](int off, std::int32_t val) noexcept {
        std::memcpy(frameBuf_.data() + off, &val, sizeof(std::int32_t));
    };

    int off{0};
    if (withBox_) {
        poke32(off, 48);
        off += 4;
        ucOffset_ = off;
        off += 48;
        poke32(off, 48);
        off += 4;
    }
    poke32(off, coordBytes);
    off += 4;
    xOffset_ = off;
    off += (numAtoms_ * 4);
    poke32(off, coordBytes);
    off += 4;

    poke32(off, coordBytes);
    off += 4;
    yOffset_ = off;
    off += (numAtoms_ * 4);
    poke32(off, coordBytes);
    off += 4;

    poke32(off, coordBytes);
    off += 4;
    zOffset_ = off;
    off += (numAtoms_ * 4);
    poke32(off, coordBytes);
}

auto Writer::fillBox(const Box& box) noexcept -> void {
    if (!withBox_) {
        return;
    }
    const double uc[6] = {
        box.sideA,
        std::sin((kHalfPi / 90.0) * (90.0 - box.angleGamma)), // cosAB
        box.sideB,
        std::sin((kHalfPi / 90.0) * (90.0 - box.angleBeta)),  // cosAC
        std::sin((kHalfPi / 90.0) * (90.0 - box.angleAlpha)), // cosBC
        box.sideC,
    };
    std::memcpy(frameBuf_.data() + ucOffset_, uc, 48);
}

auto Writer::deinterleaveCoords(const double* __restrict__ src) noexcept -> void {
    auto* __restrict__ xPtr = reinterpret_cast<float*>(frameBuf_.data() + xOffset_);
    auto* __restrict__ yPtr = reinterpret_cast<float*>(frameBuf_.data() + yOffset_);
    auto* __restrict__ zPtr = reinterpret_cast<float*>(frameBuf_.data() + zOffset_);
    for (int idx{0}; idx < numAtoms_; ++idx) {
        xPtr[idx] = static_cast<float>(src[(3 * idx)]);
        yPtr[idx] = static_cast<float>(src[(3 * idx) + 1]);
        zPtr[idx] = static_cast<float>(src[(3 * idx) + 2]);
    }
}

auto Writer::writeAll(const void* buf, std::size_t size) -> void {
    const auto* ptr = static_cast<const char*>(buf);
    std::size_t remaining{size};
    while (remaining > 0) {
        const auto written = ::write(fileFd_, ptr, remaining);
        if (written < 0) {
            throw std::runtime_error{"dcd::Writer: write() failed"};
        }
        ptr += written;
        remaining -= static_cast<std::size_t>(written);
    }
}

auto Writer::buildHeader() const -> std::array<char, kHeaderBytes> {
    std::array<char, kHeaderBytes> buf{};
    int cur{0};

    auto pI32 = [&](std::int32_t val) noexcept {
        std::memcpy(buf.data() + cur, &val, sizeof(std::int32_t));
        cur += static_cast<int>(sizeof(std::int32_t));
    };
    auto pF32 = [&](float val) noexcept {
        std::memcpy(buf.data() + cur, &val, sizeof(float));
        cur += static_cast<int>(sizeof(float));
    };
    auto pRaw = [&](const void* src, int len) noexcept {
        std::memcpy(buf.data() + cur, src, static_cast<std::size_t>(len));
        cur += len;
    };

    //  Block 1: CORD record (84 bytes of content)
    pI32(84);
    pRaw("CORD", 4);
    pI32(0); // [8]  NFRAMES
    pI32(0); // [12] ISTART
    pI32(1); // [16] NSAVC
    pI32(0); // [20] NSTEP
    for (int idx{0}; idx < 4; ++idx) {
        pI32(0);
    } // [24–35] reserved
    pI32(0);                  // [40] NAMNF = 0
    pF32(1.0F);               // [44] DELTA (float, CHARMM style)
    pI32((withBox_) ? 1 : 0); // [48] DCD_HAS_EXTRA_BLOCK
    pI32(0);                  // [52] DCD_HAS_4DIMS = 0
    for (int idx{0}; idx < 7; ++idx) {
        pI32(0);
    } // [56–80] reserved
    pI32(kCharmmVersion); // [84] signals CHARMM format to readers
    pI32(84);             // block-1 record end

    //  Block 2: title strings (2 × 80 bytes)
    pI32(164);
    pI32(2); // NTITLE
    char titleLine[80]{};
    std::strncpy(titleLine, "Created by dcd::Writer", 79);
    pRaw(titleLine, 80);
    char tsLine[80]{};
    const time_t now{std::time(nullptr)};
    std::snprintf(tsLine, 80, "REMARKS Created %.24s", std::ctime(&now));
    pRaw(tsLine, 80);
    pI32(164);

    //  Block 3: atom count
    pI32(4);
    pI32(numAtoms_);
    pI32(4);

    return buf;
}

auto Writer::finalize() noexcept -> void {
    const auto count = static_cast<std::int32_t>(frameCount_);
    if (::pwrite(fileFd_, &count, sizeof(std::int32_t), kNframesOffset) < 0) {
    }
    if (::pwrite(fileFd_, &count, sizeof(std::int32_t), kNstepOffset) < 0) {
    }
    ::close(fileFd_);
    fileFd_ = -1;
}

} // namespace dcd