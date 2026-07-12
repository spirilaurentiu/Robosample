#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <unordered_map>

#include "BatScaling.hpp"
#include "Context.hpp"
#include "OpenMMContext.hpp"
#include "ReplicaExchange.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

namespace py = pybind11;

PYBIND11_MODULE(robo_bindings, m) {
    m.doc() = "Robosample C++ bindings (robo_bindings)";

    py::enum_<NonbondedMethod>(m, "NonbondedMethod")
        .value("NoCutoff", NonbondedMethod::NoCutoff)
        .value("CutoffNonPeriodic", NonbondedMethod::CutoffNonPeriodic)
        .value("CutoffPeriodic", NonbondedMethod::CutoffPeriodic)
        .value("Ewald", NonbondedMethod::Ewald)
        .value("PME", NonbondedMethod::PME);

    py::enum_<AcceptRejectMode>(m, "AcceptRejectMode")
        .value("AlwaysAccept", AcceptRejectMode::AlwaysAccept)
        .value("MetropolisHastings", AcceptRejectMode::MetropolisHastings);

    // Proposal type of a world's sampler (a docking world auto-selects RigidKick).
    py::enum_<MoveType>(m, "MoveType")
        .value("MdHmc", MoveType::MdHmc)
        .value("RigidKick", MoveType::RigidKick)
        .value("NcmcSwitch", MoveType::NcmcSwitch);

    // Velocity-distortion option for the HMC momentum draw (simtk "NMA scaling")
    // or the position (BAT bond/angle) scaling drive (spec B4/D7). Pass
    // distort_option=DistortOption.NMA / .ScaleBendStretch to add_sampler;
    // None (default) = off. ScaleBendStretch requires mdSteps=0 (D7/INV-8 hard
    // SHALL, enforced by World::add_sampler / applyBatScalingDrive).
    py::enum_<DistortOption>(m, "DistortOption")
        .value("NMA", DistortOption::NMA, "Velocity distortion at momentum-draw time.")
        .value("ScaleBendStretch",
               DistortOption::ScaleBendStretch,
               "Deterministic BAT bond/angle position-scaling drive (B4/D7); requires mdSteps=0.");

    py::enum_<JointType>(m, "JointType")
        .value("Rigid", JointType::Rigid, "No mobility across the joint (Weld): 0 dof.")
        .value("Torsion", JointType::Torsion, "Rotation about the bond axis (the canonical dihedral): 1 dof.")
        .value("Slider", JointType::Slider, "Translation along the bond axis: 1 dof.")
        .value("Cylinder", JointType::Cylinder, "Rotation + translation about/along the bond axis: 2 dof.")
        .value("BendStretch",
               JointType::BendStretch,
               "Rotation perpendicular to the bond + translation along it: 2 dof.")
        .value("Cartesian", JointType::Cartesian, "3 translations: 3 dof.")
        .value("Ball", JointType::Ball, "Rotation, q = 4 (quaternion) by default: 3 dof.")
        .value("SphericalCoords", JointType::SphericalCoords, "BAT (azimuth, zenith, radius): 3 dof.")
        .value("FreeLine",
               JointType::FreeLine,
               "2 rotations (no spin about own line) + 3 translations, q = 7: 5 dof.")
        .value("Free", JointType::Free, "q = 7 (quaternion + translation): 6 dof.");

    py::class_<Selection>(m, "Selection").def(py::init<>());

    // A flexibility world. add_sampler returns *this so Python can chain. Note:
    // boostMDSteps was removed (it was never read). New optional kwargs:
    //   sphere_radius : docking binding-sphere radius [nm] (RigidKick worlds).
    //   use_fixman    : include the Fixman potential + logSineSqr in the
    //                   acceptance H of torsional (internal-coordinate) worlds.
    py::class_<World>(m, "World")
        .def("configure_ncmc",
             py::overload_cast<int, int, int, double>(&World::configureNcmc),
             py::arg("atom_begin"),
             py::arg("atom_end"),
             py::arg("ncmc_steps"),
             py::arg("hold_fraction") = 0.0,
             "Make this a per-molecule NCMC world: soften [atom_begin,atom_end) x "
             "rest nonbonded during a lambda:1->0->1 switch. Call AFTER add_sampler. "
             "Convenience overload for a single contiguous Region A; see "
             "configure_ncmc_region for an arbitrary atom-index set.")
        .def("configure_ncmc_region",
             py::overload_cast<std::vector<int>, int, double>(&World::configureNcmc),
             py::arg("atom_indices"),
             py::arg("ncmc_steps"),
             py::arg("hold_fraction") = 0.0,
             "Make this a per-molecule NCMC world with Region A given as an "
             "ARBITRARY atom-index set (global/OpenMM order; need not be "
             "contiguous -- docs/specs/ncmc-explicit-solvent/"
             "30-region-and-protocol-policy.md Sec.2). Sorted+deduplicated "
             "internally. Call AFTER add_sampler.")
        .def("set_ncmc_construction_ii",
             &World::setNcmcConstructionII,
             py::arg("on"),
             "Select Construction II (Metropolized-dynamics NCMC, docs/specs/"
             "ncmc-explicit-solvent/10-acceptance-construction.md): each fixed-"
             "lambda propagate substep is Metropolis accept/reject'd against the "
             "full H_lambda (inner GHMC), and the OUTER move accepts on the "
             "protocol work W alone (min(1,exp(-beta*W))) instead of the endpoint "
             "Hend-Hstart. Removes bath shadow work from acceptance -- the fix for "
             "near-zero NCMC acceptance in explicit solvent / large baths. Default "
             "off == Construction I (endpoint-DeltaH, unchanged). Call AFTER "
             "configure_ncmc/configure_ncmc_region.")
        .def("add_sampler",
             &World::add_sampler,
             py::arg("timeStep"),
             py::arg("mdSteps"),
             py::arg("acceptRejectMode"),
             py::arg("use_nuts"),
             py::arg("sphere_factor") = 1.0,     // auto R = R_receptor + factor*R_ligand
             py::arg("use_fixman") = py::none(), // None => auto: on for non-Cartesian
             py::arg("always_kick") = false,     // perturb every round (else: only when COM left sphere)
             py::arg("clash_threshold") = 10.0,  // reject if |peNew|>factor*|pePre| (relative, dimensionless)
             py::arg("max_initial_kick_tries") = 0, // >0: retry kick before round 0 until clash-free
             py::arg("distort_option") =
                 py::none(), // None => no velocity distortion; DistortOption.NMA => NMA Route B
             py::arg("nma_bias_scale") = 1.0, // alpha: NMA bias magnitude in thermal sigmas
             py::return_value_policy::reference,
             "Configure this world's sampler; returns the world for chaining. "
             "sphere_factor scales the auto-sized per-ligand binding sphere "
             "(R = R_receptor + sphere_factor*R_ligand). The docking kick relocates a "
             "ligand only when its COM leaves the sphere (always_kick=True perturbs every "
             "round). A proposed pose is rejected -- in ALL modes, including AlwaysAccept "
             "-- if its potential energy is non-finite or |PE| exceeds clash_threshold, so "
             "overlapping geometry never passes. use_fixman=None auto-enables "
             "Fixman+logSineSqr on non-Cartesian worlds. "
             "max_initial_kick_tries>0 enables a pre-round-0 retry loop that keeps "
             "drawing random placements until a clash-free starting pose is found "
             "(dPE <= maxStartPE), or raises RuntimeError after the budget is exhausted. "
             "distort_option=DistortOption.NMA draws the HMC momentum from a symmetric "
             "Gaussian mixture biased by +/- nma_bias_scale*uhat (uhat = unit NMA "
             "direction); detailed balance is preserved by a matching ln-cosh kinetic "
             "term. None (default) leaves the draw a plain Gaussian. "
             "nma_bias_scale (alpha, default 1.0) is the directed push in thermal-sigma "
             "units along uhat: the bias injects ~1/2 RT alpha^2 of directed energy, so "
             "alpha trades proposal boldness against acceptance. Guidance: alpha in "
             "[0.3, 1.0] is gentle (acceptance close to plain HMC); 1.0-3.0 is bolder; "
             ">5 collapses acceptance under a real Metropolis test. alpha=0 reproduces "
             "plain HMC. NOTE: until real soft-mode factors are supplied, uhat is the "
             "(physically meaningless) unit direction of the all-ones uScaleFactors, so "
             "the bias is safe (alpha-controlled) but not yet a useful soft-mode push.")
        // NEW: kinetic-metric mass scaling (fictitious mass; SAMPLING only).
        // OFF by default (scale 1.0 == physical). Call after the world is created.
        .def("set_mass_scale_by_joint",
             &World::setMassScaleByJoint,
             py::arg("joint_type"),
             py::arg("scale"),
             "Inflate the spatial inertia used ONLY in the proposal (momentum draw, "
             "KE, Fixman ln det M) for every body of the given JointType, raising the "
             "stable dt ~sqrt(scale) with zero configurational bias. scale=1.0 is "
             "physical (off). Typical: set_mass_scale_by_joint(JointType.Free, 16.0) "
             "on a solvent world to tame water libration.")
        .def("set_body_mass_scale",
             &World::setBodyMassScale,
             py::arg("body"),
             py::arg("scale"),
             "Per-body form of set_mass_scale_by_joint (scale=1.0 is physical/off).")
        .def("set_reversibility_check",
             &World::setReversibilityCheck,
             py::arg("interval"),
             "Low-level setter for the periodic reversibility probe. Prefer the "
             "reversibility_check_every=N argument to context.add_*_world(), which "
             "forwards here. interval<=0 disables it (default); interval>0 runs "
             "RobotEngine::checkReversibility every `interval` rounds (starting at round 0, "
             "so it doubles as a startup check), integrating mdSteps forward+back at this "
             "world's timestep from the freshly seeded state and logging the relative "
             "round-trip residual (~1e-12 ideal; a large or non-finite value means the "
             "timestep is too large for the current geometry). Non-destructive; it is a "
             "smoke test for the current configuration only -- the always-on guard is the "
             "per-step corrector throw in the integrator. See THEORY 5.7.")
        .def("set_ncmc_teleport",
             &World::setNcmcTeleport,
             py::arg("on"),
             "Enable the NCMC lambda=0 trough teleport: a rigid, KE-preserving long-range "
             "reposition of the decoupled region at the ghost trough (explicit-solvent analogue "
             "of the docking RigidKick). Default off. Active only for an acyclic Free-root region "
             "with a docking site.")
        .def("set_nma_soft_mode_from_hessian",
             &World::setNMASoftModeFromHessian,
             py::arg("atom_pos_ground"),
             py::arg("h") = 1e-5,
             py::arg("zero_tol") = 1e-6,
             "Build the mass-weighted internal-coordinate Hessian at the given minimized "
             "Ground-frame coords (nm, flat x,y,z, global atom order) and load the softest "
             "non-trivial mode into the world's NMA uScaleFactors. Call AFTER add_sampler "
             "with distort_option=DistortOption.NMA. Returns omega^2 of the chosen mode.")
        // ---- BAT-scaling drive (docs/specs/replica-exchange-nonequilibrium-
        // work.md B4/D5/D6/D7 -- Stage 2a) ----------------------------------
        // Python boundary uses FLAT [x0,y0,z0,x1,...] doubles (nm, global atom
        // order), matching set_nma_soft_mode_from_hessian's convention -- no
        // new bound value type.
        .def(
            "preview_bat_scaling",
            [](const World& w,
               const std::vector<double>& atomPosFlat,
               double s,
               const std::unordered_map<int, double>& anchorR,
               const std::unordered_map<int, double>& anchorTheta) {
                const auto n = atomPosFlat.size() / 3;
                std::vector<robo::Vec3> pos(n);
                for (std::size_t a = 0; a < n; ++a) {
                    pos[a] = robo::Vec3(atomPosFlat[(3 * a) + 0], atomPosFlat[(3 * a) + 1], atomPosFlat[(3 * a) + 2]);
                }
                const World::BatScalingResult result = w.previewBatScaling(pos, s, anchorR, anchorTheta);
                std::vector<double> outFlat(result.atomPos.size() * 3);
                for (std::size_t a = 0; a < result.atomPos.size(); ++a) {
                    outFlat[(3 * a) + 0] = result.atomPos[a][0];
                    outFlat[(3 * a) + 1] = result.atomPos[a][1];
                    outFlat[(3 * a) + 2] = result.atomPos[a][2];
                }
                return py::make_tuple(outFlat, result.nScaled, result.lnJac);
            },
            py::arg("atom_pos_ground"),
            py::arg("s"),
            py::arg("anchor_r") = std::unordered_map<int, double>{},
            py::arg("anchor_theta") = std::unordered_map<int, double>{},
            "Side-effect-free preview of the deterministic BAT-scaling drive (B4/D7): "
            "scale this world's D5-selected bond/angle DOFs by s about the given, "
            "atom-index-keyed anchors (INV-9's shared/frozen anchor snapshot -- "
            "Context.bat_anchor_snapshot()), WITHOUT mutating this world. Returns "
            "(scaled_atom_pos_ground_flat, n_scaled, ln_jac) -- the same D6 "
            "lnJac = (J(x')-J(x0)) + n_scaled*ln(s) apply_bat_scaling_drive commits.")
        .def("apply_bat_scaling_drive",
             &World::applyBatScalingDrive,
             py::arg("s"),
             py::arg("anchor_r") = std::unordered_map<int, double>{},
             py::arg("anchor_theta") = std::unordered_map<int, double>{},
             "Apply the BAT-scaling drive to THIS world's CURRENT geometry (mutates "
             "state; the world's q/frames are re-fit so the driven endpoint x^tau=x' is "
             "immediately readable). THROWS unless mdSteps==0 (D7/INV-8) and "
             "distort_option==DistortOption.ScaleBendStretch.")
        .def("get_distort_jacobian_det_log",
             &World::getDistortJacobianDetLog,
             "D6 lnJac of the last apply_bat_scaling_drive() call on this world (0.0 "
             "if never driven).")
        .def("get_last_n_scaled",
             &World::getLastNScaled,
             "D5 N_scaled of the last apply_bat_scaling_drive() call on this world.")
        .def_property_readonly("index", &World::index)
        .def_property_readonly("is_cartesian", &World::isCartesian)
        .def_property_readonly("is_docking", &World::isDocking)
        .def_property_readonly("n_dof",
             &World::nDof,
             "Number of velocity/momentum coordinates actually drawn for this world "
             "(ensemble-validation foundations spec, Sec. 2). Internal/torsional world: "
             "nu - n_C (summed joint mobilities minus removed loop-closure velocity "
             "constraints, n_C == 0 on an acyclic molecule). Cartesian world: OpenMM's own "
             "count 3*N_real - constraints - 3*[center-of-mass removal active], NOT "
             "model.nu (which is not the physical dof for a Cartesian world).")
        .def_property_readonly("num_loop_constraints",
             &World::numLoopConstraints,
             "Number of ring-closing loop DistanceConstraints this world's molecule graph "
             "produced (0 on every acyclic molecule). Structural guard for the cyclic "
             "Fixman path (docs/specs/fixman-idealized-chains-validation.md T2.0).")
        .def("current_constraint_log_det",
             &World::currentConstraintLogDet,
             "Loop-closure Fixman term ln det(G M^-1 G^T) at the CURRENT geometry (0 on an "
             "acyclic world). No sampling round required; realizes position + articulated-"
             "body inertias at the state's current q first.")
        .def("set_root_mobility",
             &World::setRootMobility,
             py::arg("molecule_index"),
             py::arg("mobility"),
             "Override ONE molecule's root attachment to Ground for THIS world only "
             "(rebuilds this world's model). Root mobility is a per-world property; "
             "use this (not a Context-level setter) to change an individual molecule. "
             "Call before add_sampler.")
        .def("set_root_mobilities",
             &World::setRootMobilities,
             py::arg("root_mobilities"),
             "Replace this world's WHOLE root-mobility vector (one entry per molecule) "
             "and rebuild ONCE -- the efficient path when many molecules change at once "
             "(e.g. a solvation shell over thousands of waters). Call before add_sampler.")
        .def("set_cartesian_solvent",
             &World::setCartesianSolvent,
             py::arg("atom_indices"),
             "Solvent-relaxing NCMC: mark these atoms (global/OpenMM order) to be "
             "advanced in FLAT Cartesian space by velocity-Verlet driven by OpenMM "
             "forces INSIDE the proposal, so the contact environment relaxes during "
             "the move instead of being a welded wall. Bodies stay welded (no Fixman/"
             "Jacobian contribution); only the per-atom (x,v) move. Empty == the "
             "welded engine, bit-for-bit. Call AFTER add_sampler/configure_ncmc.")
        .def("enable_reaction_reporter",
             &World::enableReactionReporter,
             py::arg("report_free_bodies") = true,
             py::arg("include_openmm") = true,
             py::arg("include_reaction") = false,
             "Mark this world the per-rigid-body force reporter "
             "(docs/specs/reaction-force-monitoring.md).\n"
             "\n"
             "WHAT IS RECORDED\n"
             "  For each interesting rigid body, ONE (force, torque) pair, about the body "
             "origin Bo, in the Ground frame -- the SUM of whichever term(s) below are "
             "enabled. Both terms are spatial forces about the SAME point (Bo, in Ground), so "
             "summing them is a valid spatial-force addition, not an apples-to-oranges "
             "combination. The output is ALWAYS the same 10-column CSV row shape (one force + "
             "torque per body); which term(s) went into the sum is a per-reporter-world "
             "CHOICE, not a separate set of columns.\n"
             "\n"
             "  include_openmm=True (default) adds the OpenMM NET APPLIED spatial force "
             "bodyForceG:\n"
             "      force  = sum over the body's atoms of the OpenMM per-atom force f_a\n"
             "      torque = sum over the body's atoms of (r_a - Bo) x f_a\n"
             "  i.e. the net EXTERNAL force the OpenMM force field exerts on the whole rigid "
             "body (van der Waals + electrostatics + the bonded terms crossing the body "
             "boundary). It is NOT the mobilizer/joint REACTION (the internal constraint "
             "force transmitted through the inboard joint). Being a sum of per-atom forces "
             "at fixed positions, bodyForceG is a pure function of the configuration q and "
             "carries NO velocity/momentum dependence.\n"
             "\n"
             "  include_reaction=True (default False) adds the STATIC (u=0) mobilizer "
             "REACTION: the constraint spatial force transmitted through the body's inboard "
             "joint, evaluated with all generalized speeds set to zero on the SAME accepted "
             "q. Unlike bodyForceG (this body's own applied load only), the reaction is "
             "CUMULATIVE over the body's entire outboard subtree and depends on the "
             "kinematic rooting, not just this body's own atoms -- it is a different, "
             "complementary quantity, not a refinement of bodyForceG. When False (default), "
             "no forward-dynamics step is taken and u is never touched.\n"
             "\n"
             "  DEFAULT (include_openmm=True, include_reaction=False) reproduces the "
             "OpenMM-only output this reporter shipped with, byte-for-byte. Setting "
             "include_openmm=False, include_reaction=True records the reaction ALONE (no "
             "applied-force term in the sum); setting BOTH True records their SUM (\"net\", "
             "in the sense of net field load + net constraint load). Raises if BOTH are "
             "False (nothing to report).\n"
             "\n"
             "WHICH BODIES\n"
             "  Auto-derived from THIS world's final model: every internally flexed body "
             "plus its parent (Ground itself excluded). report_free_bodies=False (default "
             "True) also drops FREE-FLOATING RIGID bodies -- a body that is BOTH "
             "Ground-rooted AND childless, i.e. a lone rigid molecule on a Free root (e.g. a "
             "Free-rooted lipid with no internal DOFs). A flexed molecule's OWN root (e.g. a "
             "receptor's root TM body) is Ground-rooted too but HAS children, so it is kept "
             "either way.\n"
             "\n"
             "WHEN IT IS SAMPLED\n"
             "  At the END of a round: AFTER the world's HMC move has drawn its random "
             "Boltzmann momenta, integrated the trajectory, and resolved accept/reject. The "
             "q read is therefore the ACCEPTED conformation -- the SAME configuration written "
             "to that round's DCD frame, so force rows join the trajectory 1:1 by frame "
             "index. It is taken at the trajectory-write cadence (write_freq, production "
             "only), not every round.\n"
             "\n"
             "  Beginning-vs-end / before-vs-after the momentum draw: the snapshot is taken "
             "at the END of the round (post-integration, post-accept), hence AFTER the "
             "momentum draw. But EVERY enabled term is velocity-free by construction -- "
             "bodyForceG is a sum of per-atom forces at fixed positions, and the reaction "
             "term is forced to u=0 before being evaluated -- so the recorded value is "
             "INDEPENDENT of the drawn momenta and of where in the ballistic-vs-diffusive "
             "trajectory the round happened to end -- evaluating it before or after the draw, "
             "at the same accepted q, yields the same result. (This invariance is why the "
             "reaction term is forced to u=0 rather than read at the LIVE trajectory-endpoint "
             "u; an instantaneous u != 0 reaction WOULD depend on the momentum draw and the "
             "trajectory point.)\n"
             "\n"
             "  The snapshot is READ-ONLY w.r.t. the sampler (it does not consume or perturb "
             "q/u), so the next round's move is bit-for-bit unaffected.\n"
             "\n"
             "USAGE\n"
             "  Prefer the want_spatial_force_history=True argument to context.add_*_world(); "
             "call this directly ONLY when the world is rebuilt after construction (e.g. "
             "add_ncmc_world / set_root_mobilities), so the derivation sees the FINAL body "
             "indexing. Raises on a Cartesian world (internal-coordinate body indexing is not "
             "meaningful there).")
        .def_property_readonly("is_reaction_reporter",
             &World::isReactionReporter,
             "Whether enable_reaction_reporter()/setReactionReporter() has been called on "
             "this world.");

    py::class_<OpenMMContext::ForceGroupEnergy>(m, "ForceGroupEnergy")
        .def_readonly("group", &OpenMMContext::ForceGroupEnergy::group)
        .def_readonly("name", &OpenMMContext::ForceGroupEnergy::name)
        .def_readonly("energy", &OpenMMContext::ForceGroupEnergy::energy);

    // -- SystemTopology (SoA) ------------------------------------------------
    py::class_<SystemTopology>(m, "SystemTopology")
        .def(py::init<>())
        .def_readwrite("num_molecules", &SystemTopology::numMolecules)
        .def_readwrite("atoms_begin", &SystemTopology::atomsBegin)
        .def_readwrite("atoms_end", &SystemTopology::atomsEnd)
        .def_readwrite("bonds_begin", &SystemTopology::bondsBegin)
        .def_readwrite("bonds_end", &SystemTopology::bondsEnd)
        .def_readwrite("angles_begin", &SystemTopology::anglesBegin)
        .def_readwrite("angles_end", &SystemTopology::anglesEnd)
        .def_readwrite("periodic_torsions_begin", &SystemTopology::periodicTorsionsBegin)
        .def_readwrite("periodic_torsions_end", &SystemTopology::periodicTorsionsEnd)
        .def_readwrite("harmonic_torsions_begin", &SystemTopology::harmonicTorsionsBegin)
        .def_readwrite("harmonic_torsions_end", &SystemTopology::harmonicTorsionsEnd)
        .def_readwrite("z_matrix_begin", &SystemTopology::zMatrixBegin)
        .def_readwrite("z_matrix_end", &SystemTopology::zMatrixEnd)
        .def_readwrite("urey_bradley_begin", &SystemTopology::ureyBradleyBegin)
        .def_readwrite("urey_bradley_end", &SystemTopology::ureyBradleyEnd)
        .def_readwrite("scaling14_begin", &SystemTopology::scaling14Begin)
        .def_readwrite("scaling14_end", &SystemTopology::scaling14End)
        .def_readwrite("exclusion_begin", &SystemTopology::exclusionBegin)
        .def_readwrite("exclusion_end", &SystemTopology::exclusionEnd)
        .def_readwrite("atoms_root_index", &SystemTopology::atomsRootIndex)
        .def_readwrite("root_mobilities", &SystemTopology::rootMobilities)
        .def_readwrite("num_atoms", &SystemTopology::numAtoms)
        .def_readwrite("atoms_unique_name", &SystemTopology::atomsUniqueName)
        .def_readwrite("atoms_nonbonded_index", &SystemTopology::atomsNonbondedIndex)
        .def_readwrite("atoms_prmtop_index", &SystemTopology::atomsPrmtopIndex)
        .def_readwrite("atoms_element_name", &SystemTopology::atomsElementName)
        .def_readwrite("atoms_element_symbol", &SystemTopology::atomsElementSymbol)
        .def_readwrite("atoms_mass", &SystemTopology::atomsMass)
        .def_readwrite("atoms_charge", &SystemTopology::atomsCharge)
        .def_readwrite("atoms_sigma", &SystemTopology::atomsSigma)
        .def_readwrite("atoms_epsilon", &SystemTopology::atomsEpsilon)
        .def_readwrite("atoms_radius", &SystemTopology::atomsRadius)
        .def_readwrite("atoms_screen", &SystemTopology::atomsScreen)
        .def_readwrite("atoms_x", &SystemTopology::atomsX)
        .def_readwrite("atoms_y", &SystemTopology::atomsY)
        .def_readwrite("atoms_z", &SystemTopology::atomsZ)
        .def_readwrite("atoms_atomic_number", &SystemTopology::atomsAtomicNumber)
        .def_readwrite("atoms_num_bonds_involved", &SystemTopology::atomsNumBondsInvolved)
        .def_readwrite("num_bonds", &SystemTopology::numBonds)
        .def_readwrite("bonds_i", &SystemTopology::bondsI)
        .def_readwrite("bonds_j", &SystemTopology::bondsJ)
        .def_readwrite("bonds_molecule_index", &SystemTopology::bondsMoleculeIndex)
        .def_readwrite("bonds_ring_closing", &SystemTopology::bondsRingClosing)
        .def_readwrite("bonds_stiffness", &SystemTopology::bondsStiffness)
        .def_readwrite("bonds_equilibrium", &SystemTopology::bondsEquilibrium)
        .def_readwrite("num_angles", &SystemTopology::numAngles)
        .def_readwrite("angles_i", &SystemTopology::anglesI)
        .def_readwrite("angles_j", &SystemTopology::anglesJ)
        .def_readwrite("angles_k", &SystemTopology::anglesK)
        .def_readwrite("angles_equilibrium", &SystemTopology::anglesEquilibrium)
        .def_readwrite("angles_stiffness", &SystemTopology::anglesStiffness)
        .def_readwrite("num_periodic_torsions", &SystemTopology::numPeriodicTorsions)
        .def_readwrite("periodic_torsions_improper", &SystemTopology::periodicTorsionsImproper)
        .def_readwrite("periodic_torsions_i", &SystemTopology::periodicTorsionsI)
        .def_readwrite("periodic_torsions_j", &SystemTopology::periodicTorsionsJ)
        .def_readwrite("periodic_torsions_k", &SystemTopology::periodicTorsionsK)
        .def_readwrite("periodic_torsions_l", &SystemTopology::periodicTorsionsL)
        .def_readwrite("periodic_torsions_n", &SystemTopology::periodicTorsionsN)
        .def_readwrite("periodic_torsions_phase", &SystemTopology::periodicTorsionsPhase)
        .def_readwrite("periodic_torsions_stiffness", &SystemTopology::periodicTorsionsStiffness)
        .def_readwrite("num_harmonic_torsions", &SystemTopology::numHarmonicTorsions)
        .def_readwrite("harmonic_torsions_i", &SystemTopology::harmonicTorsionsI)
        .def_readwrite("harmonic_torsions_j", &SystemTopology::harmonicTorsionsJ)
        .def_readwrite("harmonic_torsions_k", &SystemTopology::harmonicTorsionsK)
        .def_readwrite("harmonic_torsions_l", &SystemTopology::harmonicTorsionsL)
        .def_readwrite("harmonic_torsions_stiffness", &SystemTopology::harmonicTorsionsStiffness)
        .def_readwrite("harmonic_torsions_phase", &SystemTopology::harmonicTorsionsPhase)
        .def_readwrite("num_z_matrix_rows", &SystemTopology::numZMatrixRows)
        .def_readwrite("z_matrix_i", &SystemTopology::zMatrixI)
        .def_readwrite("z_matrix_j", &SystemTopology::zMatrixJ)
        .def_readwrite("z_matrix_k", &SystemTopology::zMatrixK)
        .def_readwrite("z_matrix_l", &SystemTopology::zMatrixL)
        .def_readwrite("num_urey_bradley", &SystemTopology::numUreyBradley)
        .def_readwrite("urey_bradley_i", &SystemTopology::ureyBradleyI)
        .def_readwrite("urey_bradley_k", &SystemTopology::ureyBradleyK)
        .def_readwrite("urey_bradley_stiffness", &SystemTopology::ureyBradleyStiffness)
        .def_readwrite("urey_bradley_equilibrium", &SystemTopology::ureyBradleyEquilibrium)
        .def_readwrite("num_scaling14", &SystemTopology::numScaling14)
        .def_readwrite("scaling14_i", &SystemTopology::scaling14I)
        .def_readwrite("scaling14_l", &SystemTopology::scaling14L)
        .def_readwrite("scaling14_charge_product", &SystemTopology::scaling14ChargeProduct)
        .def_readwrite("scaling14_epsilon", &SystemTopology::scaling14Epsilon)
        .def_readwrite("scaling14_sigma", &SystemTopology::scaling14Sigma)
        .def_readwrite("num_exclusions", &SystemTopology::numExclusions)
        .def_readwrite("exclusion_i", &SystemTopology::exclusionI)
        .def_readwrite("exclusion_j", &SystemTopology::exclusionJ)
        .def_readwrite("cmap_grid_size", &SystemTopology::cmapGridSize)
        .def_readwrite("cmap_grid_energy", &SystemTopology::cmapGridEnergy)
        .def_readwrite("cmap_torsion_map_index", &SystemTopology::cmapTorsionMapIndex)
        .def_readwrite("cmap_torsion_a1", &SystemTopology::cmapTorsionA1)
        .def_readwrite("cmap_torsion_a2", &SystemTopology::cmapTorsionA2)
        .def_readwrite("cmap_torsion_a3", &SystemTopology::cmapTorsionA3)
        .def_readwrite("cmap_torsion_a4", &SystemTopology::cmapTorsionA4)
        .def_readwrite("cmap_torsion_b1", &SystemTopology::cmapTorsionB1)
        .def_readwrite("cmap_torsion_b2", &SystemTopology::cmapTorsionB2)
        .def_readwrite("cmap_torsion_b3", &SystemTopology::cmapTorsionB3)
        .def_readwrite("cmap_torsion_b4", &SystemTopology::cmapTorsionB4)
        .def_readwrite("has_nb_fix", &SystemTopology::hasNBfix)
        .def_readwrite("num_nb_types", &SystemTopology::numNBTypes)
        .def_readwrite("a_coef", &SystemTopology::aCoef)
        .def_readwrite("b_coef", &SystemTopology::bCoef)
        .def_readwrite("use_gbsa_obc2", &SystemTopology::useGBSAOBC2)
        .def_readwrite("gbsa_solvent_dielectric", &SystemTopology::gbsaSolventDielectric)
        .def_readwrite("gbsa_solute_dielectric", &SystemTopology::gbsaSoluteDielectric)
        .def_readwrite("nonbonded_method", &SystemTopology::nonbondedMethod)
        .def_readwrite("nonbonded_cutoff", &SystemTopology::nonbondedCutoff)
        .def_readwrite("box_vectors", &SystemTopology::boxVectors)
        .def_readwrite("ewald_error_tolerance", &SystemTopology::ewaldErrorTolerance)
        .def_readwrite("num_virtual_sites", &SystemTopology::numVirtualSites)
        .def_readwrite("vs_site", &SystemTopology::vsSite)
        .def_readwrite("vs_atom1", &SystemTopology::vsAtom1)
        .def_readwrite("vs_atom2", &SystemTopology::vsAtom2)
        .def_readwrite("vs_atom3", &SystemTopology::vsAtom3)
        .def_readwrite("vs_weight1", &SystemTopology::vsWeight1)
        .def_readwrite("vs_weight2", &SystemTopology::vsWeight2)
        .def_readwrite("vs_weight3", &SystemTopology::vsWeight3)
        .def_readwrite("thermostat_temperature", &SystemTopology::thermostatTemperature)
        .def_readwrite("collision_frequency", &SystemTopology::collisionFrequency)
        .def_readwrite("seed", &SystemTopology::seed);

    // -- BAT-scaling shared anchor (spec INV-9, D2 revision 2) ---------------
    // A frozen snapshot of Context's running-mean anchor, keyed by scaled-body
    // atom index (World.preview_bat_scaling / apply_bat_scaling_drive's
    // anchor_r/anchor_theta arguments). Context-owned, NOT per-ThermodynamicState
    // -- see Context.accumulate_bat_anchor_stats / .bat_anchor_snapshot.
    py::class_<robo::BatAnchorStats::Snapshot>(m, "BatAnchorSnapshot")
        .def_readonly("mean_r", &robo::BatAnchorStats::Snapshot::meanR)
        .def_readonly("mean_theta", &robo::BatAnchorStats::Snapshot::meanTheta);

    // -- Replica exchange (docs/specs/replica-exchange-nonequilibrium-work.md,
    //    Interface I1/I3) ------------------------------------------------------
    // Exchange acceptance rule (B9). Stage 1: REMC (label-swap parallel
    // tempering). Stage 2b: RENE/REBASONTOP (driven BAT-scaling exchange) are
    // fully wired through run_rex_label_swap. RENEMC's acceptance formula is
    // wired in attempt_rex_swap, but its OWN driven round-loop (the velocity/
    // NMA drive segment) is a Stage 2c TODO -- run_rex_label_swap THROWS for
    // RunType.RENEMC (attempt_rex_swap does not; call it directly to exercise
    // RENEMC's acceptance algebra in isolation). NONE of the Stage 2b/2c C++
    // has been compiled or run (coordinator directive, 2026-07-12) -- treat
    // it as reviewed-on-paper only until a build confirms it.
    py::enum_<RUN_TYPE>(m, "RunType")
        .value("DEFAULT", RUN_TYPE::Default, "No exchange; independent replicas.")
        .value("REMC",
               RUN_TYPE::REMC,
               "Replica Exchange MC (parallel tempering): accept on -Δβ·ΔU.")
        .value("RENEMC",
               RUN_TYPE::RENEMC,
               "Replica Exchange Non-Equilibrium MC (volume-preserving velocity drive): "
               "-Δβ·ΔU on driven endpoints, no Jacobian (INV-10). Acceptance formula wired "
               "(attempt_rex_swap); the velocity/NMA driven round-loop is Stage 2c -- "
               "run_rex_label_swap throws if selected.")
        .value("RENE",
               RUN_TYPE::RENE,
               "Replica Exchange Non-Equilibrium (BAT-scaling drive): accept on nonequilibrium "
               "work -(W_X+W_Y), includes lnJac; driven worlds run mdSteps=0 (INV-8). Stage 2b: "
               "fully wired (WORK_* accumulation, F4 atomic commit, INV-7/INV-10 guards).")
        .value("REBASONTOP",
               RUN_TYPE::REBASONTOP,
               "RENE work-swaps plus interleaved REMC neighbour swaps layered on top (D4). "
               "Stage 2b: fully wired (set_interleave_remc_every/set_rebasontop_subrounds "
               "configure the interleave).");

    // Exchange topology (B7): Neighboring pairs adjacent thermodynamic states
    // with alternating parity (the Stage 1 default); All draws random pairs.
    py::enum_<ReplicaMixingScheme>(m, "ReplicaMixingScheme")
        .value("All", ReplicaMixingScheme::All)
        .value("Neighboring", ReplicaMixingScheme::Neighboring);

    // -- Context -------------------------------------------------------------
    py::class_<Context>(m, "Context")
        .def(py::init<std::string, std::uint32_t>(), py::arg("base_name"), py::arg("seed"))
        .def_readwrite("system_topology", &Context::systemTopology)
        .def("initialize",
             &Context::initialize,
             py::arg("temperatures") = std::vector<double>{},
             "Build OpenMM, set the replica temperature ladder, seed coordinates.")
        .def("add_cartesian_world",
             &Context::addCartesianWorld,
             py::arg("want_spatial_force_history") = false,
             py::return_value_policy::reference_internal,
             "Add a Cartesian (OpenMM-MD) world; returns it for .add_sampler(...). "
             "want_spatial_force_history=True always raises (docs/specs/"
             "reaction-force-monitoring.md Sec.3): a Cartesian world's internal-coordinate "
             "body indexing is not meaningful.")
        .def("add_robotic_world",
             &Context::addRoboticWorld,
             py::arg("selection"),
             py::arg("want_spatial_force_history") = false,
             py::return_value_policy::reference_internal,
             "Add an internal-coordinate (torsional) world for the given selection. "
             "want_spatial_force_history=True flags this world the net-applied-per-body-force "
             "reporter (docs/specs/reaction-force-monitoring.md); its selection's flexed "
             "bodies (plus parents) become the interesting-body set whose bodyForceG is "
             "streamed to <base>.<replica>.reactions.csv.")
        .def("add_torsional_world",
             &Context::addRoboticWorld,
             py::arg("selection"),
             py::arg("want_spatial_force_history") = false,
             py::return_value_policy::reference_internal,
             "Alias of add_robotic_world.")
        .def("add_docking_world",
             &Context::addDockingWorld,
             py::arg("ligand_molecule_indices"),
             py::return_value_policy::reference_internal,
             "Add a rigid-body docking world: the listed ligand molecules get Free roots, "
             "every other molecule is welded/rigid. Chain .add_sampler(sphere_radius=...).")
        .def("build_flexibilities",
             &Context::buildFlexibilities,
             py::arg("bonds"),
             py::arg("mobility"),
             py::arg("flag"),
             "Build a per-bond mobility selection (bonds=None => all eligible).")
        .def("run_rex",
             &Context::runREX,
             py::arg("equil_rounds"),
             py::arg("prod_rounds"),
             py::arg("write_freq"),
             py::arg("verbose"),
             "COORDINATE-swap replica exchange (legacy): Gibbs sweep over worlds + adjacent "
             "swaps, temperature-only. Retained as the INVARIANT-EQUIV oracle for "
             "run_rex_label_swap (docs/specs/replica-exchange-nonequilibrium-work.md).")
        .def("run_rex_label_swap",
             &Context::RunREX,
             py::arg("run_type"),
             py::arg("equil_rounds"),
             py::arg("prod_rounds"),
             py::arg("write_freq"),
             py::arg("verbose"),
             "LABEL-swap replica exchange (Replica/ThermodynamicState object model, docs/specs/"
             "replica-exchange-nonequilibrium-work.md B6-B7). RunType.DEFAULT/REMC/RENE/"
             "REBASONTOP are fully wired; RunType.RENEMC raises (its driven round-loop is "
             "Stage 2c -- its acceptance formula is exercised via attempt_rex_swap directly). "
             "Output CSV/DCD files are indexed by thermodynamic-state, matching run_rex's "
             "convention. NOT compiled/run since Stage 2b landed (coordinator directive) -- "
             "treat as reviewed-on-paper.")
        .def("attempt_rex_swap",
             &Context::attemptREXSwap,
             py::arg("thermo_c"),
             py::arg("thermo_h"),
             "Attempt one label swap between thermodynamic states thermo_c/thermo_h (B6: "
             "ETerm_equal for REMC, ETerm_nonequil for RENEMC, WTerm for RENE/REBASONTOP). "
             "Requires run_rex_label_swap (or setup_replica_exchange-equivalent state) to have "
             "been built first.")
        .def("check_inv7_and_inv10_guards",
             &Context::checkInv7AndInv10Guards,
             py::arg("run_type"),
             "INV-7/V9 (Fixman-on in every non-Cartesian sampler) and INV-10 (drive/run-type "
             "pairing) preconditions for a driven run type. No-op for DEFAULT/REMC. Reads only "
             "the configured worlds (no run_rex_label_swap needed first) -- callable directly "
             "to test the guard in isolation.")
        .def("set_replica_mixing_scheme",
             &Context::setReplicaMixingScheme,
             py::arg("scheme"),
             "Select ReplicaMixingScheme.Neighboring (default) or .All for run_rex_label_swap.")
        .def("set_swap_every",
             &Context::setSwapEvery,
             py::arg("n"),
             "Attempt an exchange mix only every n-th round of run_rex_label_swap (default 1).")
        .def("set_n_swap_attempts",
             &Context::setNSwapAttempts,
             py::arg("n"),
             "Number of random state pairs drawn per mix under ReplicaMixingScheme.All.")
        .def("set_swap_fixman",
             &Context::setSwapFixman,
             py::arg("enabled"),
             "OFF-by-default diagnostic flag (D3): Fixman never enters the swap acceptance "
             "(INV-7) regardless of this setting; stored for port-target interface parity (I3).")
        .def("set_interleave_remc_every",
             &Context::setInterleaveRemcEvery,
             py::arg("n"),
             "REBASONTOP (D4): run the interleaved REMC sub-round every n driven rounds (default 10).")
        .def("set_rebasontop_subrounds",
             &Context::setRebasontopSubrounds,
             py::arg("n"),
             "REBASONTOP (D4): number of alternating-parity REMC sub-rounds per interleave "
             "(default 6, matching the original's own count).")
        .def("attempted_swaps_matrix",
             &Context::attemptedSwapsMatrix,
             "T x T symmetric matrix of attempted swaps per thermodynamic-state pair.")
        .def("accepted_swaps_matrix",
             &Context::acceptedSwapsMatrix,
             "T x T symmetric matrix of accepted swaps per thermodynamic-state pair.")
        .def("accumulate_bat_anchor_stats",
             &Context::accumulateBatAnchorStats,
             py::arg("world"),
             "Update the shared/global BAT-scaling anchor (INV-9) from `world`'s CURRENT "
             "committed geometry. Call only after an EQUILIBRIUM move (never on a driven "
             "world's output -- that would feed the anchor from nonequilibrium samples).")
        .def("bat_anchor_snapshot",
             &Context::batAnchorSnapshot,
             "Frozen BatAnchorSnapshot of the current running means (INV-9): take ONE per "
             "round and pass its mean_r/mean_theta dicts to every drive that round.")
        .def("reset_bat_anchor_stats", &Context::resetBatAnchorStats, "Clear the running-mean anchor.")
        .def("set_mts",
             &Context::setMTS,
             py::arg("enabled"),
             py::arg("inner_substeps") = 4,
             "Enable r-RESPA multiple-timestep OpenMM MD (Cartesian world): slow forces once per "
             "outer step, fast bonded forces inner_substeps times. Call before initialize().")
        .def("initialize_openmm",
             &Context::initializeOpenMM,
             "Build the single OpenMM System/Context from system_topology.")
        .def("calc_openmm_potential_energy",
             &Context::calcOpenMMPotentialEnergy,
             "Set the reference coordinates and return the OpenMM potential energy [kJ/mol].")
        .def(
            "set_separate_force_groups",
            [](Context& /*ctx*/, bool enabled) {
                OpenMMContext::get().setSeparateForceGroups(enabled);
            },
            py::arg("enabled"),
            "Enable/disable separate OpenMM force groups for each Force.")
        .def(
            "set_enforce_periodic_box",
            [](Context& /*ctx*/, bool enabled) {
                OpenMMContext::get().setEnforcePeriodicBox(enabled);
            },
            py::arg("enabled"),
            "Whether OpenMM wraps coordinates into the primary box when state is "
            "pulled back. MUST stay False (the default) under explicit solvent so "
            "the robot engine receives whole molecules; energies/forces are "
            "unaffected (minimum image is always applied internally).")
        .def("calc_openmm_potential_energy_by_group",
             &Context::computePotentialEnergyByGroup,
             "Compute potential energy by OpenMM force group, returning (group, name, energy) tuples.");
}