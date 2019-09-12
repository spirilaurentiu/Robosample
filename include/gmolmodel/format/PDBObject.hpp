#ifndef __PDBOBJECT_HPP__
#define __PDBOBJECT_HPP__

#ifndef TARGET_TYPE
#define TARGET_TYPE double
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



class PDBObject{


  public:

      // TO DO:
      // header data ?
      // betafactor
      // aminoacids
      // error message

      void readPDBfile(std::string pdbfile);
      void writePDBfile(std::string pdbfile,
                                    std::vector<std::string> AtomsName,
                                    std::vector<std::string> AtomsMolName,
                                    std::vector<std::string> AtomsChain,
                                    std::vector<int> AtomsMolNo,
                                    std::vector<TARGET_TYPE> AtomsXcoord,
                                    std::vector<TARGET_TYPE> AtomsYcoord,
                                    std::vector<TARGET_TYPE> AtomsZcoord,
                                    std::vector<std::string> Atoms1Character);



      int getNumberAtoms();
      std::string getAtomsName(int p);
      std::string getAtomsMolName(int p);
      std::string getAtomsChain(int p);
      std::string getAtoms1Character(int p);
      int getAtomsMolNo(int p);

      TARGET_TYPE getAtomsXcoord(int p);
      TARGET_TYPE getAtomsYcoord(int p);
      TARGET_TYPE getAtomsZcoord(int p);

      std::ifstream pdbinput;

      int NumberAtoms;
      std::vector<std::string> AtomsName;
      std::vector<std::string> AtomsMolName;
      std::vector<std::string> AtomsChain;
      std::vector<int> AtomsMolNo;
      std::vector<TARGET_TYPE> AtomsXcoord;
      std::vector<TARGET_TYPE> AtomsYcoord;
      std::vector<TARGET_TYPE> AtomsZcoord;
      std::vector<std::string> Atoms1Character;


      std::string line;

};


#endif // __PDBOBJECT_HPP__
