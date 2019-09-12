#include "readAmberInput.hpp"


using namespace std;

int main(int argc, char **argv){

readAmberInput MOL;

//string inpcrdfile = "../../tests_inputs/lin5/ligand.inpcrd";
//string prmtopfile = "../../tests_inputs/lin5/ligand.prmtop";
string prmtopfile = argv[1];
string inpcrdfile = argv[2];

MOL.readAmberFiles(inpcrdfile, prmtopfile);



printf(">>> Check readAMBER input files: %s %s \n\n\n", inpcrdfile.c_str(), prmtopfile.c_str());

printf(">>> Check ATOMS >>> \n\n");

printf("NumberAtoms: %i \n", MOL.getNumberAtoms());
printf("NumberBonds: %i \n", MOL.getNumberBonds());
printf("NumberAngles: %i \n", MOL.getNumberAngles());
printf("NumberDihedrals: %i \n", MOL.getNumberDihedrals());


printf("\n");


printf(" %4s \t %4s \t %4s \t %11s \t %11s \t %11s \t %11s \t %11s \t %11s \t %11s \t %11s \n", "index", "Name", "Name2", "AtomsXcoord", "AtomsYcoord", "AtomsZcoord", "AtomsMass", "AtomsCharge", "AtomsRadii", "AtomsRVdW", "AtomsEpsilon");
for(int i=0; i < MOL.getNumberAtoms(); i++)
{
    printf("%4i \t %4s \t %4s \t %11f \t %11f \t %11f \t %11f \t %11f \t %11f \t %11f \t %11f \n", i, MOL.getAtomsName(i).c_str(), MOL.getAtomsNameAlias(i).c_str(), MOL.getAtomsXcoord(i), MOL.getAtomsYcoord(i), MOL.getAtomsZcoord(i), MOL.getAtomsMass(i), MOL.getAtomsCharge(i), MOL.getAtomsRadii(i), MOL.getAtomsRVdW(i), MOL.getAtomsEpsilon(i));
}

printf("\n\n");

printf(">>> Check BONDS >>> \n\n");

printf("%4s \t %4s \t %4s \t %4s \t %4s \t %11s \t %11s \n", "No", "ind1", "ind2", "name1", "name2", "bond_forceK", "bond_eqval");
for(int i=0; i < MOL.getNumberBonds(); i++)
{
    printf("%4i \t %4i \t %4i \t %4s \t %4s \t %11f \t %11f \n", i, MOL.getBondsAtomsIndex1(i), MOL.getBondsAtomsIndex2(i), MOL.getAtomsName(MOL.getBondsAtomsIndex1(i)).c_str(), MOL.getAtomsName(MOL.getBondsAtomsIndex2(i)).c_str(), MOL.getBondsForceK(i), MOL.getBondsEqval(i));
}

printf("\n\n");



printf(">>> Check ANGLES >>> \n\n");

printf("%4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %11s \t %11s \n", "No", "ind1", "ind2", "ind3", "name1", "name2", "name3", "angle_forceK", "angle_eqval");
for(int i=0; i < MOL.getNumberAngles(); i++)
{
    printf("%4i \t %4i \t %4i \t %4i \t %4s \t %4s \t %4s \t %11f \t %11f \n", i,  MOL.getAnglesAtomsIndex1(i), MOL.getAnglesAtomsIndex2(i), MOL.getAnglesAtomsIndex3(i), MOL.getAtomsName(MOL.getAnglesAtomsIndex1(i)).c_str(), MOL.getAtomsName(MOL.getAnglesAtomsIndex2(i)).c_str(), MOL.getAtomsName(MOL.getAnglesAtomsIndex3(i)).c_str(), MOL.getAnglesForceK(i), MOL.getAnglesEqval(i));
}

printf("\n\n");



printf(">>> Check DIHEDRALS >>> \n\n");

printf("%4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %4s \t %11s \t %11s \t %11s  \n", "No", "ind1", "ind2", "ind3", "ind4", "name1", "name2", "name3", "name4", "dih_forceK", "dih_period", "dih_phase");
for(int i=0; i < MOL.getNumberDihedrals(); i++)
{
    printf("%4i \t %4i \t %4i \t %4i \t %4i \t %4s \t %4s \t %4s \t %4s \t %11f \t %11f \t %11f  \n", i, MOL.getDihedralsAtomsIndex1(i), MOL.getDihedralsAtomsIndex2(i), MOL.getDihedralsAtomsIndex3(i), MOL.getDihedralsAtomsIndex4(i), MOL.getAtomsName(MOL.getDihedralsAtomsIndex1(i)).c_str(), MOL.getAtomsName(MOL.getDihedralsAtomsIndex2(i)).c_str(), MOL.getAtomsName(MOL.getDihedralsAtomsIndex3(i)).c_str(), MOL.getAtomsName(MOL.getDihedralsAtomsIndex4(i)).c_str(), MOL.getDihedralsForceK(i), MOL.getDihedralsPeriod(i), MOL.getDihedralsPhase(i));
}

printf("\n\n");

printf(">>> Check NonBondedAtomsMatrix >>> \n\n");

// MOL.TestNonBondedAtomsMatrix();

for(int i=0; i < MOL.getNumberAtoms(); i++)
{
  for(int j=0; j < MOL.getNumberAtoms(); j++)
  {
    printf(" %i ", MOL.getNonBondedAtomsMatrix(i, j));
  }
  printf("\n");
}

printf("\n\n");

for(int i=0; i < MOL.getNumberAtoms(); i++)
{
  for(int j=0; j < MOL.getNumberAtoms(); j++)
  {
    if( MOL.getNonBondedAtomsMatrix(i, j) == 0)
      printf("EXCLUDE_PAIR %5i %5i \n", i, j);
  }

}
printf("\n\n");

return 0;

}
