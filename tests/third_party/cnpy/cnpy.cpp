// ============================================================================
//  cnpy.cpp -- non-template half of the vendored cnpy (see cnpy.hpp banner
//  for provenance/license/API notes). Implements the ZIP/NPZ container
//  (local file headers, central directory, end-of-central-directory) and
//  the plain-.npy header parser. zlib is used only for crc32 (writing) and
//  inflate (reading, in case a fixture was ever re-saved through
//  numpy.savez_compressed by a reviewer) -- never for writing compressed
//  data ourselves (§9 / cnpy.hpp banner: every entry this writer produces is
//  STORE, method 0).
// ============================================================================
#include "cnpy.hpp"

#include <zlib.h>

#include <cctype>
#include <cstring>
#include <fstream>

namespace cnpy {
namespace detail {

std::string shapeTuple(const std::vector<size_t>& shape) {
    std::ostringstream os;
    os << '(';
    for (size_t i = 0; i < shape.size(); ++i) {
        os << shape[i];
        if (shape.size() == 1 || i + 1 < shape.size()) {
            os << ", ";
        }
    }
    os << ')';
    return os.str();
}

namespace {

void put16(std::vector<char>& buf, std::uint16_t v) {
    buf.push_back(static_cast<char>(v & 0xFF));
    buf.push_back(static_cast<char>((v >> 8) & 0xFF));
}
void put32(std::vector<char>& buf, std::uint32_t v) {
    buf.push_back(static_cast<char>(v & 0xFF));
    buf.push_back(static_cast<char>((v >> 8) & 0xFF));
    buf.push_back(static_cast<char>((v >> 16) & 0xFF));
    buf.push_back(static_cast<char>((v >> 24) & 0xFF));
}
std::uint16_t get16(const char* p) {
    return static_cast<std::uint16_t>(static_cast<unsigned char>(p[0])) |
           (static_cast<std::uint16_t>(static_cast<unsigned char>(p[1])) << 8);
}
std::uint32_t get32(const char* p) {
    return static_cast<std::uint32_t>(static_cast<unsigned char>(p[0])) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(p[1])) << 8) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(p[2])) << 16) |
           (static_cast<std::uint32_t>(static_cast<unsigned char>(p[3])) << 24);
}

constexpr std::uint32_t kLocalFileHeaderSig = 0x04034b50;
constexpr std::uint32_t kCentralDirSig = 0x02014b50;
constexpr std::uint32_t kEndOfCentralDirSig = 0x06054b50;

} // namespace
} // namespace detail

void NpzWriter::save(const std::string& path) const {
    std::vector<char> out;
    struct DirEntry {
        std::string name;
        std::uint32_t crc;
        std::uint32_t size;
        std::uint32_t offset;
    };
    std::vector<DirEntry> dir;
    dir.reserve(entries_.size());

    for (const Entry& e : entries_) {
        const std::uint32_t offset = static_cast<std::uint32_t>(out.size());
        const std::uint32_t crc = static_cast<std::uint32_t>(
            crc32(0L, reinterpret_cast<const Bytef*>(e.npyBytes.data()), static_cast<uInt>(e.npyBytes.size())));
        const std::uint32_t size = static_cast<std::uint32_t>(e.npyBytes.size());
        const std::string entryName = e.name + ".npy";

        detail::put32(out, detail::kLocalFileHeaderSig);
        detail::put16(out, 20);          // version needed
        detail::put16(out, 0);           // flags
        detail::put16(out, 0);           // compression method: 0 = stored
        detail::put16(out, 0);           // mod time
        detail::put16(out, 0x21);        // mod date (arbitrary, fixed -- byte-reproducible)
        detail::put32(out, crc);
        detail::put32(out, size);        // compressed size == uncompressed (stored)
        detail::put32(out, size);        // uncompressed size
        detail::put16(out, static_cast<std::uint16_t>(entryName.size()));
        detail::put16(out, 0);           // extra field length
        out.insert(out.end(), entryName.begin(), entryName.end());
        out.insert(out.end(), e.npyBytes.begin(), e.npyBytes.end());

        dir.push_back(DirEntry{entryName, crc, size, offset});
    }

    const std::uint32_t cdStart = static_cast<std::uint32_t>(out.size());
    for (const DirEntry& d : dir) {
        detail::put32(out, detail::kCentralDirSig);
        detail::put16(out, 20);          // version made by
        detail::put16(out, 20);          // version needed
        detail::put16(out, 0);           // flags
        detail::put16(out, 0);           // compression method
        detail::put16(out, 0);           // mod time
        detail::put16(out, 0x21);        // mod date
        detail::put32(out, d.crc);
        detail::put32(out, d.size);
        detail::put32(out, d.size);
        detail::put16(out, static_cast<std::uint16_t>(d.name.size()));
        detail::put16(out, 0); // extra field length
        detail::put16(out, 0); // comment length
        detail::put16(out, 0); // disk number start
        detail::put16(out, 0); // internal file attrs
        detail::put32(out, 0); // external file attrs
        detail::put32(out, d.offset);
        out.insert(out.end(), d.name.begin(), d.name.end());
    }
    const std::uint32_t cdSize = static_cast<std::uint32_t>(out.size()) - cdStart;

    detail::put32(out, detail::kEndOfCentralDirSig);
    detail::put16(out, 0); // disk number
    detail::put16(out, 0); // disk with central dir
    detail::put16(out, static_cast<std::uint16_t>(dir.size()));
    detail::put16(out, static_cast<std::uint16_t>(dir.size()));
    detail::put32(out, cdSize);
    detail::put32(out, cdStart);
    detail::put16(out, 0); // comment length

    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f) {
        throw std::runtime_error("cnpy::NpzWriter::save: cannot open '" + path + "' for writing");
    }
    f.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!f) {
        throw std::runtime_error("cnpy::NpzWriter::save: write failed for '" + path + "'");
    }
}

namespace {

std::vector<char> readWholeFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) {
        throw std::runtime_error("cnpy: cannot open '" + path + "'");
    }
    const std::streamsize n = f.tellg();
    f.seekg(0);
    std::vector<char> buf(static_cast<size_t>(n));
    if (n > 0) {
        f.read(buf.data(), n);
    }
    if (!f && n > 0) {
        throw std::runtime_error("cnpy: read failed for '" + path + "'");
    }
    return buf;
}

// Raw-deflate inflate (ZIP streams have no zlib/gzip wrapper), for the rare
// case a fixture was re-saved compressed (numpy.savez_compressed) -- this
// writer never produces compressed entries itself (see cnpy.hpp banner).
std::vector<char> inflateRaw(const char* src, size_t srcLen, size_t dstLen) {
    std::vector<char> dst(dstLen);
    z_stream strm{};
    if (inflateInit2(&strm, -MAX_WBITS) != Z_OK) {
        throw std::runtime_error("cnpy: inflateInit2 failed");
    }
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(src));
    strm.avail_in = static_cast<uInt>(srcLen);
    strm.next_out = reinterpret_cast<Bytef*>(dst.data());
    strm.avail_out = static_cast<uInt>(dstLen);
    const int rc = inflate(&strm, Z_FINISH);
    inflateEnd(&strm);
    if (rc != Z_STREAM_END) {
        throw std::runtime_error("cnpy: inflate failed (raw deflate stream)");
    }
    return dst;
}

} // namespace

NpyArray parseNpy(const std::vector<char>& npyBytes) {
    if (npyBytes.size() < 10 || npyBytes[0] != '\x93' || std::memcmp(npyBytes.data() + 1, "NUMPY", 5) != 0) {
        throw std::runtime_error("cnpy::parseNpy: bad NPY magic");
    }
    const std::uint8_t major = static_cast<std::uint8_t>(npyBytes[6]);
    size_t headerLenFieldSize;
    size_t headerLen;
    size_t dataStart;
    if (major == 1) {
        headerLenFieldSize = 2;
        headerLen = detail::get16(npyBytes.data() + 8);
        dataStart = 10 + headerLen;
    } else {
        // v2.0/3.0: 4-byte header length field.
        headerLenFieldSize = 4;
        headerLen = detail::get32(npyBytes.data() + 8);
        dataStart = 12 + headerLen;
    }
    (void)headerLenFieldSize;
    if (dataStart > npyBytes.size()) {
        throw std::runtime_error("cnpy::parseNpy: truncated header");
    }
    const std::string header(npyBytes.data() + (dataStart - headerLen), headerLen);

    NpyArray arr;
    // descr: '<f8' style token.
    {
        const size_t p = header.find("'descr':");
        const size_t q1 = header.find('\'', p + 8);
        const size_t q2 = header.find('\'', q1 + 1);
        if (p == std::string::npos || q1 == std::string::npos || q2 == std::string::npos) {
            throw std::runtime_error("cnpy::parseNpy: cannot find 'descr' in header");
        }
        const std::string descr = header.substr(q1 + 1, q2 - q1 - 1);
        // Word size is the trailing digit run of the descr token, e.g.
        // "<f8" -> 8, "<i4" -> 4, "|u1" -> 1.
        size_t digitsStart = descr.size();
        while (digitsStart > 0 && std::isdigit(static_cast<unsigned char>(descr[digitsStart - 1]))) {
            --digitsStart;
        }
        if (digitsStart == descr.size()) {
            throw std::runtime_error("cnpy::parseNpy: cannot parse word size from descr '" + descr + "'");
        }
        arr.word_size = static_cast<size_t>(std::stoi(descr.substr(digitsStart)));
    }
    // fortran_order
    {
        const size_t p = header.find("'fortran_order':");
        const size_t shapeP = header.find("'shape'");
        if (p == std::string::npos || shapeP == std::string::npos) {
            throw std::runtime_error("cnpy::parseNpy: cannot find 'fortran_order' in header");
        }
        const std::string between = header.substr(p, shapeP - p);
        arr.fortran_order = between.find("True") != std::string::npos;
    }
    // shape
    {
        const size_t p = header.find("'shape':");
        const size_t open = header.find('(', p);
        const size_t close = header.find(')', open);
        if (p == std::string::npos || open == std::string::npos || close == std::string::npos) {
            throw std::runtime_error("cnpy::parseNpy: cannot find 'shape' in header");
        }
        std::string dims = header.substr(open + 1, close - open - 1);
        size_t pos = 0;
        while (pos < dims.size()) {
            while (pos < dims.size() && (dims[pos] == ' ' || dims[pos] == ',')) {
                ++pos;
            }
            if (pos >= dims.size()) {
                break;
            }
            size_t end = pos;
            while (end < dims.size() && std::isdigit(static_cast<unsigned char>(dims[end]))) {
                ++end;
            }
            if (end > pos) {
                arr.shape.push_back(static_cast<size_t>(std::stoul(dims.substr(pos, end - pos))));
            }
            pos = end + 1;
        }
    }

    const size_t nVals = arr.num_vals();
    const size_t dataBytes = nVals * arr.word_size;
    if (dataStart + dataBytes > npyBytes.size()) {
        throw std::runtime_error("cnpy::parseNpy: truncated data section");
    }
    arr.bytes.assign(npyBytes.begin() + static_cast<std::ptrdiff_t>(dataStart),
                      npyBytes.begin() + static_cast<std::ptrdiff_t>(dataStart + dataBytes));
    return arr;
}

NpyArray npy_load(const std::string& path) {
    return parseNpy(readWholeFile(path));
}

npz_t npz_load(const std::string& path) {
    const std::vector<char> whole = readWholeFile(path);
    if (whole.size() < 22) {
        throw std::runtime_error("cnpy::npz_load: file too small to be a ZIP: '" + path + "'");
    }

    // Locate the End-Of-Central-Directory record by scanning backward for its
    // signature (no ZIP comment is ever written by NpzWriter::save, so it is
    // always exactly 22 bytes from EOF, but scan defensively).
    std::ptrdiff_t eocd = -1;
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(whole.size()) - 22; i >= 0; --i) {
        if (detail::get32(whole.data() + i) == detail::kEndOfCentralDirSig) {
            eocd = i;
            break;
        }
    }
    if (eocd < 0) {
        throw std::runtime_error("cnpy::npz_load: no End-Of-Central-Directory record found in '" + path + "'");
    }
    const std::uint16_t numEntries = detail::get16(whole.data() + eocd + 10);
    const std::uint32_t cdOffset = detail::get32(whole.data() + eocd + 16);

    npz_t out;
    size_t p = cdOffset;
    for (std::uint16_t i = 0; i < numEntries; ++i) {
        if (detail::get32(whole.data() + p) != detail::kCentralDirSig) {
            throw std::runtime_error("cnpy::npz_load: malformed central directory in '" + path + "'");
        }
        const std::uint16_t method = detail::get16(whole.data() + p + 10);
        const std::uint32_t compSize = detail::get32(whole.data() + p + 20);
        const std::uint32_t uncompSize = detail::get32(whole.data() + p + 24);
        const std::uint16_t nameLen = detail::get16(whole.data() + p + 28);
        const std::uint16_t extraLen = detail::get16(whole.data() + p + 30);
        const std::uint16_t commentLen = detail::get16(whole.data() + p + 32);
        const std::uint32_t localOffset = detail::get32(whole.data() + p + 42);
        std::string name(whole.data() + p + 46, nameLen);
        p += 46 + nameLen + extraLen + commentLen;

        // Local file header at localOffset: 30-byte fixed part, then name +
        // extra, then the (possibly compressed) data.
        const std::uint16_t localNameLen = detail::get16(whole.data() + localOffset + 26);
        const std::uint16_t localExtraLen = detail::get16(whole.data() + localOffset + 28);
        const size_t dataOffset = localOffset + 30 + localNameLen + localExtraLen;

        std::vector<char> npyBytes;
        if (method == 0) {
            npyBytes.assign(whole.begin() + static_cast<std::ptrdiff_t>(dataOffset),
                             whole.begin() + static_cast<std::ptrdiff_t>(dataOffset + compSize));
        } else if (method == 8) {
            npyBytes = inflateRaw(whole.data() + dataOffset, compSize, uncompSize);
        } else {
            throw std::runtime_error("cnpy::npz_load: unsupported ZIP compression method " +
                                      std::to_string(method) + " for entry '" + name + "'");
        }

        if (name.size() >= 4 && name.substr(name.size() - 4) == ".npy") {
            name.resize(name.size() - 4);
        }
        out[name] = parseNpy(npyBytes);
    }
    return out;
}

NpyArray npz_load(const std::string& path, const std::string& name) {
    npz_t all = npz_load(path);
    auto it = all.find(name);
    if (it == all.end()) {
        throw std::runtime_error("cnpy::npz_load: no array named '" + name + "' in '" + path + "'");
    }
    return it->second;
}

} // namespace cnpy
