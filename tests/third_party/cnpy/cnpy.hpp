// ============================================================================
//  cnpy.hpp -- minimal single-header/source library for reading and writing
//  NumPy .npy files and .npz archives (a ZIP container of .npy entries).
//
//  Vendored for tests/fixtures/robotics_oracle
//  (docs/specs/robotics-oracle-differential.md §9,
//   docs/specs/robotics-oracle-data-provenance.md §4): the robotics-oracle
//  fixtures are stored as one <case>.npz per case, readable both by this
//  loader (disasm-side TestRoboticsOracle.cpp) and by plain `numpy.load()`
//  (spec §9's explicit reviewer requirement).
//
//  Provenance: this is a from-scratch, MIT-licensed reimplementation of the
//  well-known rogersce/cnpy (MIT) API and on-disk format, vendored rather
//  than fetched because this environment has no network access. It follows
//  cnpy's design point of using zlib ONLY for the ZIP entry's CRC32 (not for
//  compression -- every entry is stored, ZIP method 0, exactly like upstream
//  cnpy), so the .npz this writes is byte-for-byte a normal ZIP/NPZ file.
//
//  API difference from upstream cnpy (deliberate, Rule 2 -- simplicity):
//  upstream's npz_save() is incremental -- every call re-reads the existing
//  archive, appends one entry, and rewrites the whole central directory to
//  disk. A robotics-oracle case bakes dozens of named arrays (one case with
//  3 states x ~30 quantities), so that pattern would be an O(n^2) file
//  rewrite for no benefit here (we always know the full array set up front).
//  Instead, `NpzWriter` accumulates all of a case's arrays in memory and
//  writes the archive once via `save()`. The reading API (`npy_load`,
//  `npz_load`) matches upstream's shape/semantics so this is a drop-in
//  concept for anyone who knows cnpy.
//
//  Supported element types: double, float, std::int32_t, std::int64_t,
//  std::uint8_t -- the only ones the robotics-oracle fixtures need (all
//  reference/aggregate arrays are float64; nothing here needs the wider
//  dtype zoo upstream cnpy supports).
//
//  MIT License (this file)
//  Copyright (c) 2026 Robosample contributors.
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// ============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace cnpy {

namespace detail {

// NumPy dtype descriptor string for T (little-endian '<', matching the only
// architecture this project targets, x86_64 -- Robosample CLAUDE.md §Environment).
template <typename T>
inline std::string dtypeDescr() {
    if constexpr (std::is_same_v<T, double>) {
        return "<f8";
    } else if constexpr (std::is_same_v<T, float>) {
        return "<f4";
    } else if constexpr (std::is_same_v<T, std::int64_t>) {
        return "<i8";
    } else if constexpr (std::is_same_v<T, std::int32_t>) {
        return "<i4";
    } else if constexpr (std::is_same_v<T, std::uint8_t>) {
        return "|u1";
    } else {
        static_assert(!sizeof(T), "cnpy: unsupported element type for NPY dtype descriptor");
        return "";
    }
}

std::string shapeTuple(const std::vector<size_t>& shape);

} // namespace detail

// Build a complete, in-memory .npy file (magic + version + header + raw
// little-endian data) for a C-order array. Header is padded with spaces so
// the total preamble length (magic+version+header-length field+header
// string) is a multiple of 64 bytes, matching numpy's own writer (not
// required for numpy to read it back, but keeps the on-disk layout
// canonical for a human `numpy.load`/`unzip -l` inspection, per §9).
template <typename T>
std::vector<char> buildNpy(const T* data, const std::vector<size_t>& shape) {
    const size_t nVals = shape.empty() ? 0
                                        : std::accumulate(shape.begin(), shape.end(), static_cast<size_t>(1),
                                                           std::multiplies<size_t>());

    std::ostringstream hdr;
    hdr << "{'descr': '" << detail::dtypeDescr<T>() << "', 'fortran_order': False, 'shape': "
        << detail::shapeTuple(shape) << ", }";
    std::string header = hdr.str();

    // magic(6) + version(2) + header_len field(2) = 10 bytes preamble for
    // format version 1.0. Pad `header` (plus trailing '\n') so 10+header
    // total is a multiple of 64.
    const size_t preamble = 10;
    size_t padded = header.size() + 1; // +1 for the trailing '\n'
    size_t total = preamble + padded;
    size_t remainder = total % 64;
    if (remainder != 0) {
        padded += (64 - remainder);
    }
    header.resize(padded - 1, ' ');
    header.push_back('\n');

    if (header.size() > 0xFFFFu) {
        throw std::runtime_error("cnpy::buildNpy: header too large for NPY v1.0 (16-bit length field)");
    }

    std::vector<char> out;
    out.reserve(preamble + header.size() + (nVals * sizeof(T)));
    out.push_back('\x93');
    out.push_back('N');
    out.push_back('U');
    out.push_back('M');
    out.push_back('P');
    out.push_back('Y');
    out.push_back(static_cast<char>(1)); // major version
    out.push_back(static_cast<char>(0)); // minor version
    const std::uint16_t headerLen = static_cast<std::uint16_t>(header.size());
    out.push_back(static_cast<char>(headerLen & 0xFF));
    out.push_back(static_cast<char>((headerLen >> 8) & 0xFF));
    out.insert(out.end(), header.begin(), header.end());
    const char* raw = reinterpret_cast<const char*>(data);
    out.insert(out.end(), raw, raw + (nVals * sizeof(T)));
    return out;
}

// One loaded array: raw little-endian bytes + shape (C order) + word size.
// Matches upstream cnpy::NpyArray's shape (shape/word_size/data<T>()).
struct NpyArray {
    std::vector<size_t> shape;
    size_t word_size = 0;
    bool fortran_order = false;
    std::vector<char> bytes;

    size_t num_vals() const {
        if (shape.empty()) {
            return 0;
        }
        size_t n = 1;
        for (size_t d : shape) {
            n *= d;
        }
        return n;
    }

    template <typename T>
    const T* data() const {
        if (sizeof(T) != word_size) {
            throw std::runtime_error("cnpy::NpyArray::data<T>(): word size mismatch with stored dtype");
        }
        return reinterpret_cast<const T*>(bytes.data());
    }
};

using npz_t = std::map<std::string, NpyArray>;

// ---- Reading (non-template; implemented in cnpy.cpp) ----
NpyArray parseNpy(const std::vector<char>& npyBytes);
npz_t npz_load(const std::string& path);
NpyArray npz_load(const std::string& path, const std::string& name);
NpyArray npy_load(const std::string& path);

// ---- Writing ----
// Accumulates named arrays in memory, then writes one .npz (ZIP, STORE
// method, per the banner above) in a single pass.
class NpzWriter {
public:
    template <typename T>
    void add(const std::string& name, const T* data, const std::vector<size_t>& shape) {
        entries_.push_back(Entry{name, buildNpy<T>(data, shape)});
        manifestArrays_.emplace_back(name, shape);
    }

    void save(const std::string& path) const;

    // For the caller to also emit the manifest.json "array list with
    // shapes" (§9) without re-deriving it from the written file.
    const std::vector<std::pair<std::string, std::vector<size_t>>>& arrays() const {
        return manifestArrays_;
    }

private:
    struct Entry {
        std::string name;
        std::vector<char> npyBytes; // complete .npy file bytes for this array
    };
    std::vector<Entry> entries_;
    std::vector<std::pair<std::string, std::vector<size_t>>> manifestArrays_;
};

} // namespace cnpy
