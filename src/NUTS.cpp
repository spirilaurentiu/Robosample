#include "NUTS.hpp"

void NUTSNodeRef::copyVectorsFrom(const NUTSNodeRef& src) const {
    nuts_detail::copyVec(minus_q, src.minus_q, n);
    nuts_detail::copyVec(minus_p, src.minus_p, n);
    nuts_detail::copyVec(plus_q, src.plus_q, n);
    nuts_detail::copyVec(plus_p, src.plus_p, n);
    nuts_detail::copyVec(proposal_q, src.proposal_q, n);
    nuts_detail::copyVec(proposal_p, src.proposal_p, n);
}

void NUTSNodeRef::copyProposalFrom(const NUTSNodeRef& src) {
    nuts_detail::copyVec(proposal_q, src.proposal_q, n);
    nuts_detail::copyVec(proposal_p, src.proposal_p, n);
    proposedEnergy = src.proposedEnergy;
    proposalLeafIndex = src.proposalLeafIndex;
}

void NUTSWorkspaceSimbody::init(int maxDepth, int nDOF) {
    this->nDOF = nDOF;

    // Layout (in units of nDOF doubles)
    //   pool nodes : maxDepth * 2 * 6
    //   treeRoot   : 6
    //   pTemp      : 1
    //   total      : (maxDepth * 12 + 8) * nDOF
    //
    // Alignment: 64 bytes (one cache line).  With AVX-256 each cache line
    // holds 8 doubles; aligning to 64 bytes means all dot-product inner
    // loops start on a vector-unit boundary -> compiler can emit aligned
    // vmovapd / vfmadd231pd without scalar prologue.
    constexpr std::size_t kAlign = 64;
    slabDoubles = static_cast<std::size_t>((maxDepth * 12) + 8) * nDOF;
    const std::size_t rawBytes = slabDoubles * sizeof(SimTK::Real);
    const std::size_t alignedBytes = ((rawBytes + kAlign - 1) / kAlign) * kAlign;

    deallocate();
    slab = allocAligned(alignedBytes, kAlign);
    std::memset(slab, 0, alignedBytes);

    // Carve up the slab
    SimTK::Real* ptr = slab;

    auto nextVec = [&]() -> SimTK::Real* {
        SimTK::Real* p = ptr;
        ptr += nDOF;
        return p;
    };

    nodePool.resize(maxDepth);
    for (auto& pair : nodePool) {
        for (NUTSNodeRef* ref : {pair.data(), &pair[1]}) {
            ref->minus_q = nextVec();
            ref->minus_p = nextVec();
            ref->plus_q = nextVec();
            ref->plus_p = nextVec();
            ref->proposal_q = nextVec();
            ref->proposal_p = nextVec();
            ref->n = nDOF;
        }
    }

    treeRoot.minus_q = nextVec();
    treeRoot.minus_p = nextVec();
    treeRoot.plus_q = nextVec();
    treeRoot.plus_p = nextVec();
    treeRoot.proposal_q = nextVec();
    treeRoot.proposal_p = nextVec();
    treeRoot.n = nDOF;

    pTemp = nextVec();

    assert(ptr == slab + slabDoubles); // sanity: used exactly what we planned

    // SimTK adapter - one heap allocation, reused every leaf
    pTempSimTK.resize(nDOF);
}

// U-turn check
auto NUTSWorkspaceSimbody::isUTurn(const NUTSNodeRef& node) -> bool {
    SimTK::Real dot1;
    SimTK::Real dot2;
    nuts_detail::wrapDot(node.minus_q, node.plus_q, node.minus_p, node.plus_p, node.n, dot1, dot2);
    return (dot1 < 0.0) || (dot2 < 0.0);
}

auto NUTSWorkspaceSimbody::allocAligned(std::size_t bytes, std::size_t align) -> SimTK::Real* {
    auto* ptr = static_cast<SimTK::Real*>(std::aligned_alloc(align, bytes));
    if (ptr == nullptr) {
        throw std::bad_alloc{};
    }
    return ptr;
}

void NUTSWorkspaceSimbody::deallocate() {
    if (slab != nullptr) {
        std::free(slab);
        slab = nullptr;
    }
}

void NUTSTrajectoryLog::printOneLiner() const {
    static const char* stopAbbr[] = {"NotNUTS", "MaxDepth", "UTurn", "SubUTurn", "NoValid"};

    const int valid = static_cast<int>(std::count_if(events.begin(), events.end(), [](const LeafEvent& e) {
        return e.sliceValid;
    }));
    const int total = totalSteps();

    // Time extents (in whatever units timestep is expressed)
    const SimTK::Real tMin = -static_cast<SimTK::Real>(bwdSteps) * timestep;
    const SimTK::Real tMax = static_cast<SimTK::Real>(fwdSteps) * timestep;
    const SimTK::Real tSel = static_cast<SimTK::Real>(selectedStepIndex) * timestep;

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4);
    oss << "\t - [NUTS] d=" << finalDepth << " | " << stopAbbr[static_cast<int>(stopReason)] << " | " << total
        << " steps" << " (bwd:" << -bwdSteps << " fwd:+" << fwdSteps << ")" << " | t=[" << tMin << " ps, +"
        << tMax << " ps]" << " | sel=" << (tSel >= 0 ? "+" : "") << tSel << " ps | valid=" << valid << "/"
        << total;

    std::cout << oss.str() << "\n";
}

void NUTSTrajectoryLog::printTable() const {
    static const char* stopStr[] = {"NotUsingNUTS",
                                    "MaxDepth",
                                    "U-turn",
                                    "Subtree U-turn",
                                    "NoValidProposals"};

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n  ┌─ NUTS trajectory ───────────────────────────────────────────────┐\n";
    std::cout << "  │  H0 = " << std::setw(12) << H0 << "   log u = " << std::setw(12) << logU
              << "   slice Δ = " << std::setw(8) << sliceWindow << "  │\n";
    std::cout << "  │  Depth = " << finalDepth << "   Steps = " << totalSteps()
              << "   Stop: " << stopStr[static_cast<int>(stopReason)]
              << "   Selected: " << (selectedStepIndex >= 0 ? "+" : "") << selectedStepIndex << "\n";
    std::cout << "  ├────────┬──────┬──────────────┬───────┬──────────┤\n";
    std::cout << "  │  Step  │ Dir  │      H       │ Valid │ Selected │\n";
    std::cout << "  ├────────┼──────┼──────────────┼───────┼──────────┤\n";

    std::cout << "  │  +0000 │ init │ " << std::setw(12) << H0 << " │  yes  │ "
              << (selectedStepIndex == 0 ? "   ←   " : "       ") << "│\n";

    auto sorted = events;
    std::sort(sorted.begin(), sorted.end(), [](const LeafEvent& a, const LeafEvent& b) {
        return a.stepIndex < b.stepIndex;
    });

    for (const auto& e : sorted) {
        const char* dirStr = (e.direction == NUTSDirection::Forward) ? "fwd " : "bwd ";
        char stepBuf[8];
        std::snprintf(stepBuf, sizeof(stepBuf), "%+05d", e.stepIndex);
        std::cout << "  │  " << stepBuf << " │ " << dirStr << " │ " << std::setw(12) << e.H << " │  "
                  << (e.sliceValid ? "yes" : "no ") << "  │ "
                  << (e.stepIndex == selectedStepIndex ? "   ←   " : "       ") << "│\n";
    }
    std::cout << "  └────────┴──────┴──────────────┴───────┴──────────┘\n\n";
}