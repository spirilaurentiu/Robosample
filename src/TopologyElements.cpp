#include "TopologyElements.hpp"

constexpr SimTK::Angle TetrahedralAngle = 109.47 * SimTK::Deg2Rad; // radians
constexpr double DefaultBondLength = 0.19; // nm
constexpr double Cos120 = -0.5; // cos(120 degrees)
constexpr double Sin120 = 0.866025; // sin(120 degrees)

void Atom::createSingleAtom() {
	// Create a new SingleAtom compound representing this atom.
	// SingleAtom is a building block in Molmodel that can hold multiple
	// "bond centers" (points in 3D where other atoms can attach).
	compoundSingleAtom = new SimTK::Compound::SingleAtom(atomSpec.atomName, element);

	const int currAtomNBonds = getNumBondsInvolved();
	if (currAtomNBonds > 0) {
		if (currAtomNBonds == 1) {
			// For an atom with a single bond (e.g. hydrogen), define only one bond center.
			compoundSingleAtom->addFirstBondCenter("bond1", atomSpec.atomName);
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
				"bond1", "bond2", atomSpec.atomName,
				SimTK::UnitVec3(1, 0, 0),             // Along +X
				SimTK::UnitVec3(Cos120, Sin120, 0.0)  // 120 degrees rotated in XY
			);

			// --- Step 2: Add third bond (if needed) ---
			if (currAtomNBonds > 2) {
				// addLeftHandedBondCenter places the third bond center ABOVE the XY plane,
				// at tetrahedral angles to both bond1 and bond2.
				//
				// "Left-handed" means that if you curl your left hand from bond1 to bond2,
				// your thumb points in the direction of this bond. This sets a chirality.
				compoundSingleAtom->addLeftHandedBondCenter( "bond3", atomSpec.atomName, TetrahedralAngle, TetrahedralAngle
				);
			}

			// --- Step 3: Add fourth bond (if needed) ---
			if (currAtomNBonds > 3) {
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
				compoundSingleAtom->addRightHandedBondCenter("bond4", atomSpec.atomName, TetrahedralAngle, TetrahedralAngle
				);
			}
		}

		// The "inboard" bond is the one used to connect this atom into a larger structure.
		// Here, bond1 is always chosen as the inboard bond center.
		compoundSingleAtom->setInboardBondCenter("bond1");

		// Default length of this inboard bond (in nanometers, ~1.9 A).
		compoundSingleAtom->setDefaultInboardBondLength(DefaultBondLength);
	}

	// Give this SingleAtom a unique compound name within the molecule.
	compoundSingleAtom->setCompoundName(atomSpec.atomName);
}

std::string Atom::getElementName() const {
	switch (atomSpec.atomicNumber)
	{
		case 1: return "hydrogen";
		case 2: return "helium";
		case 3: return "lithium";
		case 4: return "beryllium";
		case 5: return "boron";
		case 6: return "carbon";
		case 7: return "nitrogen";
		case 8: return "oxygen";
		case 9: return "fluorine";
		case 10: return "neon";
		case 11: return "sodium";
		case 12: return "magnesium";
		case 13: return "aluminum";
		case 14: return "silicon";
		case 15: return "phosphorus";
		case 16: return "sulfur";
		case 17: return "chlorine";
		case 18: return "argon";
		case 19: return "potassium";
		case 20: return "calcium";
		case 21: return "scandium";
		case 22: return "titanium";
		case 23: return "vanadium";
		case 24: return "chromium";
		case 25: return "manganese";
		case 26: return "iron";
		case 27: return "cobalt";
		case 28: return "nickel";
		case 29: return "copper";
		case 30: return "zinc";
		case 31: return "gallium";
		case 32: return "germanium";
		case 33: return "arsenic";
		case 34: return "selenium";
		case 35: return "bromine";
		case 36: return "krypton";
		case 37: return "rubidium";
		case 38: return "strontium";
		case 39: return "yttrium";
		case 40: return "zirconium";
		case 41: return "niobium";
		case 42: return "molybdenum";
		case 43: return "technetium";
		case 44: return "ruthenium";
		case 45: return "rhodium";
		case 46: return "palladium";
		case 47: return "silver";
		case 48: return "cadmium";
		case 49: return "indium";
		case 50: return "tin";
		case 51: return "antimony";
		case 52: return "tellurium";
		case 53: return "iodine";
		case 54: return "xenon";
		case 55: return "cesium";
		case 56: return "barium";
		case 57: return "lanthanum";
		case 58: return "cerium";
		case 59: return "praseodymium";
		case 60: return "neodymium";
		case 61: return "promethium";
		case 62: return "samarium";
		case 63: return "europium";
		case 64: return "gadolinium";
		case 65: return "terbium";
		case 66: return "dysprosium";
		case 67: return "holmium";
		case 68: return "erbium";
		case 69: return "thulium";
		case 70: return "ytterbium";
		case 71: return "lutetium";
		case 72: return "hafnium";
		case 73: return "tantalum";
		case 74: return "tungsten";
		case 75: return "rhenium";
		case 76: return "osmium";
		case 77: return "iridium";
		case 78: return "platinum";
		case 79: return "gold";
		case 80: return "mercury";
		case 81: return "thallium";
		case 82: return "lead";
		case 83: return "bismuth";
		case 84: return "polonium";
		case 85: return "astatine";
		case 86: return "radon";
		case 87: return "francium";
		case 88: return "radium";
		case 89: return "actinium";
		case 90: return "thorium";
		case 91: return "protactinium";
		case 92: return "uranium";
		case 93: return "neptunium";
		case 94: return "plutonium";
		case 95: return "americium";
		case 96: return "curium";
		case 97: return "berkelium";
		case 98: return "californium";
		case 99: return "einsteinium";
		case 100: return "fermium";
		case 101: return "mendelevium";
		case 102: return "nobelium";
		case 103: return "lawrencium";
		case 104: return "rutherfordium";
		case 105: return "dubnium";
		case 106: return "seaborgium";
		case 107: return "bohrium";
		case 108: return "hassium";
		case 109: return "meitnerium";
		case 110: return "darmstadtium";
		case 111: return "roentgenium";
		case 112: return "ununbium";
		case 113: return "ununtrium";
		case 114: return "ununquadium";
		case 115: return "ununpentium";
		case 116: return "ununhexium";
		default: return "unknown";
	}
}
std::string Atom::getElementSymbol() const {
	switch (atomSpec.atomicNumber)
	{
		case 1: return "H";
		case 2: return "He";
		case 3: return "Li";
		case 4: return "Be";
		case 5: return "B";
		case 6: return "C";
		case 7: return "N";
		case 8: return "O";
		case 9: return "F";
		case 10: return "Ne";
		case 11: return "Na";
		case 12: return "Mg";
		case 13: return "Al";
		case 14: return "Si";
		case 15: return "P";
		case 16: return "S";
		case 17: return "Cl";
		case 18: return "Ar";
		case 19: return "K";
		case 20: return "Ca";
		case 21: return "Sc";
		case 22: return "Ti";
		case 23: return "V";
		case 24: return "Cr";
		case 25: return "Mn";
		case 26: return "Fe";
		case 27: return "Co";
		case 28: return "Ni";
		case 29: return "Cu";
		case 30: return "Zn";
		case 31: return "Ga";
		case 32: return "Ge";
		case 33: return "As";
		case 34: return "Se";
		case 35: return "Br";
		case 36: return "Kr";
		case 37: return "Rb";
		case 38: return "Sr";
		case 39: return "Y";
		case 40: return "Zr";
		case 41: return "Nb";
		case 42: return "Mo";
		case 43: return "Tc";
		case 44: return "Ru";
		case 45: return "Rh";
		case 46: return "Pd";
		case 47: return "Ag";
		case 48: return "Cd";
		case 49: return "In";
		case 50: return "Sn";
		case 51: return "Sb";
		case 52: return "Te";
		case 53: return "I";
		case 54: return "Xe";
		case 55: return "Cs";
		case 56: return "Ba";
		case 57: return "La";
		case 58: return "Ce";
		case 59: return "Pr";
		case 60: return "Nd";
		case 61: return "Pm";
		case 62: return "Sm";
		case 63: return "Eu";
		case 64: return "Gd";
		case 65: return "Tb";
		case 66: return "Dy";
		case 67: return "Ho";
		case 68: return "Er";
		case 69: return "Tm";
		case 70: return "Yb";
		case 71: return "Lu";
		case 72: return "Hf";
		case 73: return "Ta";
		case 74: return "W";
		case 75: return "Re";
		case 76: return "Os";
		case 77: return "Ir";
		case 78: return "Pt";
		case 79: return "Au";
		case 80: return "Hg";
		case 81: return "Tl";
		case 82: return "Pb";
		case 83: return "Bi";
		case 84: return "Po";
		case 85: return "At";
		case 86: return "Rn";
		case 87: return "Fr";
		case 88: return "Ra";
		case 89: return "Ac";
		case 90: return "Th";
		case 91: return "Pa";
		case 92: return "U";
		case 93: return "Np";
		case 94: return "Pu";
		case 95: return "Am";
		case 96: return "Cm";
		case 97: return "Bk";
		case 98: return "Cf";
		case 99: return "Es";
		case 100: return "Fm";
		case 101: return "Md";
		case 102: return "No";
		case 103: return "Lr";
		case 104: return "Rf";
		case 105: return "Db";
		case 106: return "Sg";
		case 107: return "Bh";
		case 108: return "Hs";
		case 109: return "Mt";
		case 110: return "Ds";
		case 111: return "Rg";
		case 112: return "Uub";
		case 113: return "Uut";
		case 114: return "Uuq";
		case 115: return "Uup";
		case 116: return "Uuh";
		default: return "Unknown";
	}
}
