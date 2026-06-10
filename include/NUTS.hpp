#pragma once

#include "EnergySnapshot.hpp"
#include "OpenMM.h"
#include "bgeneral.hpp"

enum class StopReason : std::uint8_t {
    NotUsingNUTS,
    MaxDepth,
    UTurn,
    SubtreeUTurn,
    NoValidProposals
};

enum class NUTSDirection : std::uint8_t {
    Backward = 0,
    Forward = 1
};

enum class NUTSCoordinates : std::uint8_t {
    Cartesian = 0,
    Torsional = 1
};

struct PhasePointOpenMM {
    std::vector<OpenMMContext::Vec3> positions;
    std::vector<OpenMMContext::Vec3> momenta;
};

struct PhasePointSimbody {
    SimTK::Vector q;
    SimTK::Vector p;
};

struct NUTSNodeSimbody {
    PhasePointSimbody minus;
    PhasePointSimbody plus;
    PhasePointSimbody proposal;
    EnergySnapshot proposedEnergy;

    int numValidSlices{0};
    int proposalLeafIndex{0};
    bool stop{false};
};

struct NUTSNodeOpenMM {
    PhasePointOpenMM minus;
    PhasePointOpenMM plus;
    PhasePointOpenMM proposal;
    EnergySnapshot proposedEnergy;

    int numValidSlices{0};
    bool stop{false};
};

struct NUTSResultSimbody {
    PhasePointSimbody proposal;
    EnergySnapshot proposedEnergy;
    StopReason stopReason;
    int depth;
};

struct NUTSResult {
    PhasePointSimbody proposalSimbody;
    PhasePointOpenMM proposalOpenMM;
    EnergySnapshot proposedEnergy;
    StopReason stopReason;
    int depth;
};


namespace nuts_detail {

static constexpr SimTK::Real TWO_PI = 6.283185307179586476925;
static constexpr SimTK::Real INV_TWO_PI = 0.15915494309189534; // 1 / (2 * pi)

inline void copyVec(SimTK::Real* __restrict__ dst, const SimTK::Real* __restrict__ src, int n) {
    std::memcpy(dst, src, static_cast<std::size_t>(n) * sizeof(SimTK::Real));
}

// Compute both U-turn dot products in a single pass
//
// The U-turn criterion (Betancourt 2017) needs:
//   d1 = dot(wrap(q+ - q-), p-)
//   d2 = dot(wrap(q+ - q-), p+)
//
// Both share the same wrapped difference (q+ - q-).  Splitting into
// wrapDot() + dot() makes two passes over the data:
//   Pass 1: wrap and write scratch,  compute acc1 += d * p-
//   Pass 2: read scratch,            compute acc2 += d * p+
//
// wrapDot() fuses them:
//   Single pass: wrap d, write scratch, compute acc1 += d*p-  AND  acc2 += d*p+
//
// The compiler emits two back-to-back vfmadd231pd - both can be issued in the
// same clock cycle on AVX2 (2 FMA execution ports). The extra cost vs a single
// wrapDot() is one FMA per 4 doubles, well within the pipeline budget.
//
// Assembly of the fused hot loop (gcc -O3 -mavx2 -mfma -ffast-math):
//
//   .L38:
//     vmovupd   (%rdx,%rax), %ymm7      ; load b[i:+3]
//     vsubpd    (%rsi,%rax), %ymm7, %ymm1  ; d = b - a
//     vmulpd    %ymm4,       %ymm1, %ymm0  ; d * INV_TWO_PI
//     vroundpd  $4,          %ymm0, %ymm0  ; rint (single instr)
//     vfnmadd132pd %ymm3,    %ymm1, %ymm0  ; d - TWO_PI * rint
//     vfmadd231pd (%rcx,%rax), %ymm0, %ymm5  ; acc1 += d * p-
//     vfmadd231pd (%r8, %rax), %ymm0, %ymm6  ; acc2 += d * p+
//     jne .L38
//
// 8 data instructions per 4 doubles, producing BOTH dot products.
inline void wrapDot(const SimTK::Real* __restrict__ minus_q,
                    const SimTK::Real* __restrict__ plus_q,
                    const SimTK::Real* __restrict__ minus_p,
                    const SimTK::Real* __restrict__ plus_p,
                    int n,
                    SimTK::Real& out1,
                    SimTK::Real& out2) noexcept {
    SimTK::Real acc1 = 0.0;
    SimTK::Real acc2 = 0.0;

    for (int i = 0; i < n; ++i) {
        SimTK::Real delta = plus_q[i] - minus_q[i];
        delta -= TWO_PI * std::rint(delta * INV_TWO_PI);
        acc1 += delta * minus_p[i];
        acc2 += delta * plus_p[i]; // one extra vfmadd231pd, essentially free
    }
    out1 = acc1;
    out2 = acc2;
}

} // namespace nuts_detail

struct NUTSNodeRef {
    // Six raw pointers into the contiguous slab (set once at init, never change)
    SimTK::Real* minus_q = nullptr;
    SimTK::Real* minus_p = nullptr;
    SimTK::Real* plus_q = nullptr;
    SimTK::Real* plus_p = nullptr;
    SimTK::Real* proposal_q = nullptr;
    SimTK::Real* proposal_p = nullptr;
    int n = 0; // = nDOF, cached to avoid passing everywhere

    // Metadata (plain-old-data, cheap to copy by value)
    EnergySnapshot proposedEnergy{};
    int numValidSlices{0};
    bool stop{false};
    int proposalLeafIndex{0}; // signed step index of the held proposal

    // Copy all six vectors from another node - plain memcpy, no heap activity
    void copyVectorsFrom(const NUTSNodeRef& src) const;
    void copyProposalFrom(const NUTSNodeRef& src);
};


struct NUTSWorkspaceSimbody {
    // Slab
    SimTK::Real* slab{nullptr}; // single aligned allocation
    std::size_t slabDoubles{0};

    // Node pool - pointers into slab, valid after init()
    // nodePool[depth][side]:
    //   depth in [0, maxDepth)
    //   side  0 = first-half subtree, 1 = second-half subtree
    std::vector<std::array<NUTSNodeRef, 2>> nodePool;

    // Root trajectory node — its phase-point vectors also live in the slab
    NUTSNodeRef treeRoot;

    // Scratch buffers
    SimTK::Real* pTemp{nullptr}; // staging for matter->multiplyByM
    SimTK::Vector pTempSimTK;    // SimTK adapter (one allocation, kept alive)

    int nDOF{0};

    // Lifecycle
    ~NUTSWorkspaceSimbody() {
        deallocate();
    }

    // Disallow copy; allow move so the sampler can own this by value
    NUTSWorkspaceSimbody() = default;
    NUTSWorkspaceSimbody(const NUTSWorkspaceSimbody&) = delete;
    auto operator=(const NUTSWorkspaceSimbody&) -> NUTSWorkspaceSimbody& = delete;

    void init(int maxDepth, int nDOF);

    // U-turn check
    [[nodiscard]] static auto isUTurn(const NUTSNodeRef& node) -> bool;

    private:
    static auto allocAligned(std::size_t bytes, std::size_t align) -> SimTK::Real*;
    void deallocate();
};


struct NUTSTrajectoryLog {
    struct LeafEvent {
        int stepIndex; // signed: -N... -1 = backward, +1 ...+N = forward
        NUTSDirection direction;
        SimTK::Real H;
        bool sliceValid;
    };

    SimTK::Real H0{};
    SimTK::Real logU{};
    SimTK::Real sliceWindow{}; // = -logU - H0 = Exp(1) draw
    SimTK::Real timestep{};    // simulation timestep (for time-axis labels)
    int fwdSteps{0};
    int bwdSteps{0};
    int selectedStepIndex{0};
    int finalDepth{0};
    StopReason stopReason{StopReason::MaxDepth};
    std::vector<LeafEvent> events;

    [[nodiscard]] auto totalSteps() const -> int {
        return fwdSteps + bwdSteps;
    }

    void printOneLiner() const;
    void printTable() const;
};
