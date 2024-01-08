#pragma once

#include "Molmodel.h"

using namespace SimTK;

class ELEMENT_CACHE {
public:
	bool addElement(int atomicNumber, SimTK::mdunits::Mass atomicMass);
	const SimTK::Element& getElement(int atomicNumber, SimTK::mdunits::Mass atomicMass) const;
    const SimTK::Element::Symbol& getSymbolByAtomicNumber(int atomicNumber) const;
	void print() const;

private:
	inline std::size_t key(int i, int j) const;

	std::unordered_map<int, SimTK::Element> cache;
	SimTK::Element::Symbol null = "";

	// std::array is an aggregate type and requires 1 more layer of braces
	std::array<std::pair<SimTK::Element::Name, SimTK::Element::Symbol>, 119> elements = { { // Laurentiu
		{ "Hydrogen", "H" },
		{ "Helium", "He" },
		{ "Lithium", "Li" },
		{ "Beryllium", "Be" },
		{ "Boron", "B" },
		{ "Carbon", "C" },
		{ "Nitrogen", "N" },
		{ "Oxygen", "O" },
		{ "Fluorine", "F" },
		{ "Neon", "Ne" },
		{ "Sodium", "Na" },
		{ "Magnesium", "Mg" },
		{ "Aluminum", "Al" },
		{ "Silicon", "Si" },
		{ "Phosphorus", "P" },
		{ "Sulfur", "S" },
		{ "Chlorine", "Cl" },
		{ "Argon", "Ar" },
		{ "Potassium", "K" },
		{ "Calcium", "Ca" },
		{ "Scandium", "Sc" },
		{ "Titanium", "Ti" },
		{ "Vanadium", "V" },
		{ "Chromium", "Cr" },
		{ "Manganese", "Mn" },
		{ "Iron", "Fe" },
		{ "Cobalt", "Co" },
		{ "Nickel", "Ni" },
		{ "Copper", "Cu" },
		{ "Zinc", "Zn" },
		{ "Gallium", "Ga" },
		{ "Germanium", "Ge" },
		{ "Arsenic", "As" },
		{ "Selenium", "Se" },
		{ "Bromine", "Br" },
		{ "Krypton", "Kr" },
		{ "Rubidium", "Rb" },
		{ "Strontium", "Sr" },
		{ "Yttrium", "Y" },
		{ "Zirconium", "Zr" },
		{ "Niobium", "Nb" },
		{ "Molybdenum", "Mo" },
		{ "Technetium", "Tc" },
		{ "Ruthenium", "Ru" },
		{ "Rhodium", "Rh" },
		{ "Palladium", "Pd" },
		{ "Silver", "Ag" },
		{ "Cadmium", "Cd" },
		{ "Indium", "In" },
		{ "Tin", "Sn" },
		{ "Antimony", "Sb" },
		{ "Tellurium", "Te" },
		{ "Iodine", "I" },
		{ "Xenon", "Xe" },
		{ "Cesium", "Cs" },
		{ "Barium", "Ba" },
		{ "Lanthanum", "La" },
		{ "Cerium", "Ce" },
		{ "Praseodymium", "Pr" },
		{ "Neodymium", "Nd" },
		{ "Promethium", "Pm" },
		{ "Samarium", "Sm" },
		{ "Europium", "Eu" },
		{ "Gadolinium", "Gd" },
		{ "Terbium", "Tb" },
		{ "Dysprosium", "Dy" },
		{ "Holmium", "Ho" },
		{ "Erbium", "Er" },
		{ "Thulium", "Tm" },
		{ "Ytterbium", "Yb" },
		{ "Lutetium", "Lu" },
		{ "Hafnium", "Hf" },
		{ "Tantalum", "Ta" },
		{ "Tungsten", "W" },
		{ "Rhenium", "Re" },
		{ "Osmium", "Os" },
		{ "Iridium", "Ir" },
		{ "Platinum", "Pt" },
		{ "Gold", "Au" },
		{ "Mercury", "Hg" },
		{ "Thallium", "Tl" },
		{ "Lead", "Pb" },
		{ "Bismuth", "Bi" },
		{ "Polonium", "Po" },
		{ "Astatine", "At" },
		{ "Radon", "Rn" },
		{ "Francium", "Fr" },
		{ "Radium", "Ra" },
		{ "Actinium", "Ac" },
		{ "Thorium", "Th" },
		{ "Protactinium", "Pa" },
		{ "Uranium", "U" },
		{ "Neptunium", "Np" },
		{ "Plutonium", "Pu" },
		{ "Americium", "Am" },
		{ "Curium", "Cm" },
		{ "Berkelium", "Bk" },
		{ "Californium", "Cf" },
		{ "Einsteinium", "Es" },
		{ "Fermium", "Fm" },
		{ "Mendelevium", "Md" },
		{ "Nobelium", "No" },
		{ "Lawrencium", "Lr" },
		{ "Rutherfordium", "Rf" },
		{ "Dubnium", "Db" },
		{ "Seaborgium", "Sg" },
		{ "Bohrium", "Bh" },
		{ "Hassium", "Hs" },
		{ "Meitnerium", "Mt" },
		{ "Darmstadtium", "Ds" },
		{ "Roentgenium", "Rg" },
		{ "Copernicium", "Cn" },
		{ "Nihonium", "Nh" },
		{ "Flerovium", "Fl" },
		{ "Moscovium", "Mc" },
		{ "Livermorium", "Lv" },
		{ "Tennessine", "Ts" },
		{ "Oganesson", "Og" },
		{ "Dummy", "Du" } // Laurentiu
	} };
};