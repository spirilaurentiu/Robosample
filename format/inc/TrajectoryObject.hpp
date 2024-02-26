#ifndef __TRAJECTORYOBJECT_HPP__
#define __TRAJECTORYOBJECT_HPP__

#ifndef TARGET_TYPE
#define TARGET_TYPE double
#endif


extern "C" {
#include "mdtraj/dcdplugin.h"
}

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





class TrajectoryObject{

  public:

      // "dcd" - only for now
      std::string TrajectoryType;
      std::string TrajectoryFile;

      int unitcell;
      int natoms;
      int nsets;

      // ~TrajectoryObject();

      // creates trajectory header
      void createTrajectory(const std::string& TrajectoryFile,
                            const std::string& TrajectoryType,
                            int selectedatoms,
                            int selectedunitcell);

      // opens existing trajectory
      void openTrajectory(const std::string& TrajectoryFile, const std::string& TrajectoryType );

      // append timestep to created/opened trajectory
      void appendTimestep (
                          const std::string& TrajectoryType,
                          const std::vector<TARGET_TYPE>& AtomsXcoord,
                          const std::vector<TARGET_TYPE>& AtomsYcoord,
                          const std::vector<TARGET_TYPE>& AtomsZcoord
                          );


      // for DCD type
      void * v;
      dcdhandle *dcd;
      molfile_timestep_t timestep;

};


#endif // __TRAJECTORYOBJECT_HPP__
