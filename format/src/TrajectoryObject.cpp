#include "TrajectoryObject.hpp"




// TO DO :
// check mem


// TrajectoryObject::~TrajectoryObject()
// {
//   // close_file_read(v)
//   close_file_write(v);
// }

void TrajectoryObject::createTrajectory(const std::string& TrajectoryFile,
                      const std::string& TrajectoryType,
                      int selectedatoms,
                      int selectedunitcell)
{

try
{
  if (TrajectoryType == "dcd")
  {

    natoms = selectedatoms;
    unitcell = selectedunitcell;
    nsets = 0;

    // printf("selectedatoms %i ; natoms %i \n", selectedatoms, natoms);

    dcdhandle_po = open_dcd_write(TrajectoryFile.c_str(), "dcd", natoms, unitcell);
    if (!dcdhandle_po) {
      fprintf(stderr, "Error in writing DCD header \n");
    }

    dcd = (dcdhandle *)dcdhandle_po;

  }
  else{
    fprintf(stderr, "Only dcd format is supported - for now ...\n");
  }

}
catch(std::exception e){
  std::cout << "Error in generating Trajectory File\n";
  std::cout << e.what() << '\n';
}

}



void TrajectoryObject::openTrajectory(const std::string& TrajectoryFile,
                      const std::string& TrajectoryType)
{

try
{
  if (TrajectoryType == "dcd")
  {
    dcdhandle_po = open_dcd_read(TrajectoryFile.c_str(), "dcd", &natoms, &nsets);
    if (!dcdhandle_po) {
      fprintf(stderr, "Error in opening DCD file \n");
    }

    dcd = (dcdhandle *)dcdhandle_po;
  }
  else{
    fprintf(stderr, "Only dcd format is supported - for now ...\n");
  }

}

catch(std::exception e){
  std::cout << "Error in opening Trajectory File\n";
  std::cout << e.what() << '\n';
}

}


void TrajectoryObject::appendFrame (
                    const std::string& TrajectoryType,
                    const std::vector<TARGET_TYPE>& AtomsXcoord,
                    const std::vector<TARGET_TYPE>& AtomsYcoord,
                    const std::vector<TARGET_TYPE>& AtomsZcoord
                    )
{

try
{
  if (TrajectoryType == "dcd")
  {
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);

    float *pos = timestep.coords;

    for (int j=0; j < natoms; j++) {
      *(pos++) = AtomsXcoord[j];
      *(pos++) = AtomsYcoord[j];
      *(pos++) = AtomsZcoord[j];
    }

    int write_status = write_timestep(dcd, &timestep);

    if (write_status) {
        fprintf(stderr, "error in write_timestep on frame %d\n", nsets);
      }

    nsets++;

  }
  else{
    fprintf(stderr, "Only dcd format is supported - for now ...\n");
  }

}

catch(std::exception e){
  std::cout << "Error in appendng timestep to Trajectory File\n";
  std::cout << e.what() << '\n';
}

}
