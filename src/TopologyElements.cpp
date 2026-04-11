#include "TopologyElements.hpp"

static const SimTK::Angle TetrahedralAngle = 109.47 * SimTK::Deg2Rad; // radians
constexpr double DefaultBondLength = 0.19;                            // nm
constexpr double Cos120 = -0.5;                                       // cos(120 degrees)
constexpr double Sin120 = 0.866025;                                   // sin(120 degrees)

void RoboAtom::createSingleAtom() {
    // SimTK_ASSERT_ALWAYS(numBondsInvolved <= 4,
    // 	"Atom::createSingleAtom(): Atoms with more than 4 bonds are not supported.");

    // switch (numBondsInvolved)
    // {
    // case 1:
    // 	compoundSingleAtom = new SimTK::UnivalentAtom(uniqueAtomName, element);
    // 	break;
    // case 2:
    // 	compoundSingleAtom = new SimTK::BivalentAtom(uniqueAtomName, element);
    // 	break;
    // case 3:
    // 	compoundSingleAtom = new SimTK::TrivalentAtom(uniqueAtomName, element);
    // 	break;
    // case 4:
    // 	compoundSingleAtom = new SimTK::QuadrivalentAtom(uniqueAtomName, element);
    // 	break;

    // default:
    // 	break;
    // }

    // Create a new SingleAtom compound representing this atom.
    // SingleAtom is a building block in Molmodel that can hold multiple
    // "bond centers" (points in 3D where other atoms can attach).
    compoundSingleAtom = new SimTK::Compound::SingleAtom(identity.uniqueAtomName, elementInfo.element);

    if (connectivity.numBondsInvolved > 0) {
        if (connectivity.numBondsInvolved == 1) {
            // For an atom with a single bond (e.g. hydrogen), define only one bond center.
            compoundSingleAtom->addFirstBondCenter("bond1", identity.uniqueAtomName);
        } else {
            // Tetrahedral bond angle (109.47 degrees) in radians.
            // This is the ideal angle between bonds in sp3 hybridized atoms (like carbon).
            SimTK::Angle TetrahedralAngle = 109.47 * SimTK::Deg2Rad;

            // --- Step 1: Define the first two bonds ---
            // addFirstTwoBondCenters places two bond centers in space relative to the atom.
            // Here:
            //   bond1 lies on the +X axis: (1, 0, 0)
            //   bond2 lies rotated ~120 degrees in the XY plane: (-0.5, 0.866025, 0.0)
            //
            // Together, these form a "V" in the XY plane, like the two bonds in water.
            compoundSingleAtom->addFirstTwoBondCenters(
                "bond1",
                "bond2",
                identity.uniqueAtomName,
                SimTK::UnitVec3(1, 0, 0),            // Along +X
                SimTK::UnitVec3(Cos120, Sin120, 0.0) // 120 degrees rotated in XY
            );

            // --- Step 2: Add third bond (if needed) ---
            if (connectivity.numBondsInvolved > 2) {
                // addLeftHandedBondCenter places the third bond center ABOVE the XY plane,
                // at tetrahedral angles to both bond1 and bond2.
                //
                // "Left-handed" means that if you curl your left hand from bond1 to bond2,
                // your thumb points in the direction of this bond. This sets a chirality.
                compoundSingleAtom->addLeftHandedBondCenter("bond3",
                                                            identity.uniqueAtomName,
                                                            TetrahedralAngle,
                                                            TetrahedralAngle);
            }

            // --- Step 3: Add fourth bond (if needed) ---
            if (connectivity.numBondsInvolved > 3) {
                // addRightHandedBondCenter places the fourth bond BELOW the XY plane,
                // again at tetrahedral angles to bond1 and bond2.
                //
                // "Right-handed" is the mirror orientation: if you curl your right hand
                // from bond1 to bond2, your thumb points in this direction.
                //
                // Together, bond1–bond4 form a tetrahedron centered on the atom:
                // - bond1 and bond2 in XY plane
                // - bond3 pointing up (left-handed)
                // - bond4 pointing down (right-handed)
                compoundSingleAtom->addRightHandedBondCenter("bond4",
                                                             identity.uniqueAtomName,
                                                             TetrahedralAngle,
                                                             TetrahedralAngle);
            }
        }

        // The "inboard" bond is the one used to connect this atom into a larger structure.
        // Here, bond1 is always chosen as the inboard bond center.
        compoundSingleAtom->setInboardBondCenter("bond1");

        // Default length of this inboard bond (in nanometers, ~1.9 A).
        compoundSingleAtom->setDefaultInboardBondLength(DefaultBondLength);
    }

    // Give this SingleAtom a unique compound name within the molecule.
    compoundSingleAtom->setCompoundName(identity.uniqueAtomName);
}
