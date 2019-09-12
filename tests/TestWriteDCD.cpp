#include "format/PDBObject.hpp"
#include "format/TrajectoryObject.hpp"


using namespace std;



int main(int argc, char *argv[])
{

// This test takes a PDB file and creates a DCD Trajectory file.
//
// 1st argv - pdb input file
// 2nd argv - dcd output file
//


// No unitcell
int unitcell = 0;

// Create a PDBObject from the PDB file
PDBObject MOL;
MOL.readPDBfile(argv[1]);

printf("NumberAtoms %i \n", MOL.NumberAtoms);
printf("AtomsXcoord[0] %f \n", MOL.AtomsXcoord[0]);


// Create a Trajectory object
TrajectoryObject MOL_traj;

// Since it's not an existing trajectory we generate the header
MOL_traj.createTrajectory(argv[2], "dcd", MOL.NumberAtoms, unitcell);

// move the molecule a little bit
int nrsteps = 10;
for(int i=0; i<nrsteps; i++)
{
  // change coordinates:
  for(int j =0; j<MOL_traj.natoms; j++)
  {
     MOL.AtomsXcoord[j] = MOL.AtomsXcoord[j] + i;
     MOL.AtomsYcoord[j] = MOL.AtomsYcoord[j] + i;
     MOL.AtomsZcoord[j] = MOL.AtomsZcoord[j] + i;
  }

  // append new coordinates to a timestep entry
  MOL_traj.appendTimestep("dcd", MOL.AtomsXcoord, MOL.AtomsYcoord, MOL.AtomsZcoord);

}

printf("natoms %i \n", MOL_traj.natoms);
printf("nsets %i \n", MOL_traj.nsets);


return 0;

}
