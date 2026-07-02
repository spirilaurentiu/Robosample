// ============================================================================
//  sha256.hpp -- minimal, header-only, public-domain-style SHA-256
//  implementation (FIPS 180-4), vendored for the robotics-oracle manifest's
//  `sha256_npz` provenance field (docs/specs/robotics-oracle-differential.md
//  §9, docs/specs/robotics-oracle-data-provenance.md §4). Written from the
//  FIPS 180-4 specification (public constants/algorithm, no third-party
//  source copied) because this environment has no network access to fetch
//  an existing vendored implementation.
//
//  This implementation and this file are dedicated to the public domain
//  (CC0-1.0-equivalent): no rights reserved, use freely.
//
//  API: `sha256::hex(const void* data, size_t len)` -> lowercase 64-char hex
//  digest string. That is the only entry point either side of the
//  robotics-oracle tooling needs.
// ============================================================================
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace sha256 {

namespace detail {

inline std::uint32_t rotr(std::uint32_t x, std::uint32_t n) {
    return (x >> n) | (x << (32 - n));
}

// FIPS 180-4 §4.2.2 round constants.
inline const std::array<std::uint32_t, 64>& K() {
    static const std::array<std::uint32_t, 64> k = {
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
        0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
        0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
        0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
        0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
        0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};
    return k;
}

} // namespace detail

// Streaming SHA-256; `hex()` below wraps this for the common one-shot case.
class Sha256 {
public:
    Sha256() { reset(); }

    void reset() {
        h_ = {0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
              0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U};
        buffer_.clear();
        totalLen_ = 0;
    }

    void update(const void* data, std::size_t len) {
        const auto* bytes = static_cast<const unsigned char*>(data);
        totalLen_ += len;
        buffer_.insert(buffer_.end(), bytes, bytes + len);
        while (buffer_.size() >= 64) {
            processBlock(buffer_.data());
            buffer_.erase(buffer_.begin(), buffer_.begin() + 64);
        }
    }

    std::array<std::uint8_t, 32> finalize() {
        const std::uint64_t bitLen = totalLen_ * 8;
        buffer_.push_back(0x80);
        while (buffer_.size() % 64 != 56) {
            buffer_.push_back(0x00);
        }
        for (int i = 7; i >= 0; --i) {
            buffer_.push_back(static_cast<unsigned char>((bitLen >> (i * 8)) & 0xFF));
        }
        for (std::size_t off = 0; off < buffer_.size(); off += 64) {
            processBlock(buffer_.data() + off);
        }
        std::array<std::uint8_t, 32> digest{};
        for (int i = 0; i < 8; ++i) {
            digest[static_cast<std::size_t>(i * 4) + 0] = static_cast<std::uint8_t>((h_[static_cast<size_t>(i)] >> 24) & 0xFF);
            digest[static_cast<std::size_t>(i * 4) + 1] = static_cast<std::uint8_t>((h_[static_cast<size_t>(i)] >> 16) & 0xFF);
            digest[static_cast<std::size_t>(i * 4) + 2] = static_cast<std::uint8_t>((h_[static_cast<size_t>(i)] >> 8) & 0xFF);
            digest[static_cast<std::size_t>(i * 4) + 3] = static_cast<std::uint8_t>(h_[static_cast<size_t>(i)] & 0xFF);
        }
        return digest;
    }

private:
    void processBlock(const unsigned char* block) {
        std::uint32_t w[64];
        for (int t = 0; t < 16; ++t) {
            w[t] = (static_cast<std::uint32_t>(block[t * 4]) << 24) |
                   (static_cast<std::uint32_t>(block[(t * 4) + 1]) << 16) |
                   (static_cast<std::uint32_t>(block[(t * 4) + 2]) << 8) |
                   static_cast<std::uint32_t>(block[(t * 4) + 3]);
        }
        for (int t = 16; t < 64; ++t) {
            const std::uint32_t s0 = detail::rotr(w[t - 15], 7) ^ detail::rotr(w[t - 15], 18) ^ (w[t - 15] >> 3);
            const std::uint32_t s1 = detail::rotr(w[t - 2], 17) ^ detail::rotr(w[t - 2], 19) ^ (w[t - 2] >> 10);
            w[t] = w[t - 16] + s0 + w[t - 7] + s1;
        }

        std::uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
        std::uint32_t e = h_[4], f = h_[5], g = h_[6], hh = h_[7];
        const auto& k = detail::K();
        for (int t = 0; t < 64; ++t) {
            const std::uint32_t s1 = detail::rotr(e, 6) ^ detail::rotr(e, 11) ^ detail::rotr(e, 25);
            const std::uint32_t ch = (e & f) ^ ((~e) & g);
            const std::uint32_t temp1 = hh + s1 + ch + k[static_cast<size_t>(t)] + w[t];
            const std::uint32_t s0 = detail::rotr(a, 2) ^ detail::rotr(a, 13) ^ detail::rotr(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t temp2 = s0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }
        h_[0] += a; h_[1] += b; h_[2] += c; h_[3] += d;
        h_[4] += e; h_[5] += f; h_[6] += g; h_[7] += hh;
    }

    std::array<std::uint32_t, 8> h_{};
    std::vector<unsigned char> buffer_;
    std::uint64_t totalLen_ = 0;
};

// One-shot: SHA-256 of a byte buffer, as a lowercase 64-hex-digit string.
inline std::string hex(const void* data, std::size_t len) {
    Sha256 s;
    s.update(data, len);
    const std::array<std::uint8_t, 32> digest = s.finalize();
    static const char* kHexDigits = "0123456789abcdef";
    std::string out(64, '0');
    for (std::size_t i = 0; i < 32; ++i) {
        out[(i * 2) + 0] = kHexDigits[(digest[i] >> 4) & 0xF];
        out[(i * 2) + 1] = kHexDigits[digest[i] & 0xF];
    }
    return out;
}

} // namespace sha256
