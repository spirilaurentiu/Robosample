#include "readAmberInput.hpp"
#include <assert.h>
#include <algorithm>
#include <string>
#include <cctype>

constexpr int ATOM_NAME_SIZE = 4;
constexpr int AMBER_INDEX_MULTIPLIER = 3;
constexpr int AMBER_INDEX_DIFF = 1;

int numres = 0;

void AmberReader::readAmberFiles(const std::string& inpcrdfile, const std::string& prmtopfile)
{

try
{

  readInpcrd(inpcrdfile);

  prmtop.open(prmtopfile.c_str());
  if(!prmtop.is_open())
  {
          printf("Error in prmtop file : %s file not opened \n", prmtopfile.c_str());
          exit(1);
  }

  while(!prmtop.eof())
  {
    getline(prmtop, line);

    if (line.find("FLAG POINTERS") != std::string::npos)
        readPointers();
    else if (line.find("CHARGE") != std::string::npos)
        readAtomsCharge();

    else if (line.find("ATOMIC_NUMBER") != std::string::npos)
        readAtomicNumber();

    else if (line.find("RESIDUE_LABEL") != std::string::npos)
        readResidueLabels();
    else if (line.find("RESIDUE_POINTER") != std::string::npos)
        readResiduePointers();

    else if (line.find("MASS") != std::string::npos)
        readAtomsMass();
    else if (line.find("FLAG RADII") != std::string::npos)
        readAtomsRadii();
    else if (line.find("FLAG SCREEN") != std::string::npos)
        readAtomsScreenGBIS();
    else if (line.find("ATOM_TYPE_INDEX") != std::string::npos)
        readAtomsTypesIndex();
    else if (line.find("ATOM_NAME") != std::string::npos)
        readAtomsName();
    else if (line.find("AMBER_ATOM_TYPE") != std::string::npos)
        readAtomsTypes();
    else if (line.find("NONBONDED_PARM_INDEX") != std::string::npos)
        readNonbondedParmIndex();



    else if (line.find("BOND_FORCE_CONSTANT") != std::string::npos)
        readTempBondsForceK();
    else if (line.find("BOND_EQUIL_VALUE") != std::string::npos)
        readTempBondsEqval();
    else if (line.find("BONDS_INC_HYDROGEN") != std::string::npos)
        readBonds(NumberBondsH);
    else if (line.find("BONDS_WITHOUT_HYDROGEN") != std::string::npos)
        readBonds(NumberBondsNONH);

    else if (line.find("ANGLE_FORCE_CONSTANT") != std::string::npos)
        readTempAnglesForceK();
    else if (line.find("ANGLE_EQUIL_VALUE") != std::string::npos)
        readTempAnglesEqval();
    else if (line.find("ANGLES_INC_HYDROGEN") != std::string::npos)
        readAngles(NumberAnglesH);
    else if (line.find("ANGLES_WITHOUT_HYDROGEN") != std::string::npos)
        readAngles(NumberAnglesNONH);

    else if (line.find("DIHEDRAL_FORCE_CONSTANT") != std::string::npos)
        readTempDihedralsForceK();
    else if (line.find("DIHEDRAL_PERIODICITY") != std::string::npos)
        readTempDihedralsPeriod();
    else if (line.find("DIHEDRAL_PHASE") != std::string::npos)
        readTempDihedralsPhase();
    else if (line.find("SCEE_SCALE_FACTOR") != std::string::npos)
        readTempSCEEScaleFactor();
    else if (line.find("SCNB_SCALE_FACTOR") != std::string::npos)
        readTempSCNBScaleFactor();
    else if (line.find("DIHEDRALS_INC_HYDROGEN") != std::string::npos)
        readDihedrals(NumberDihedralsH);
    else if (line.find("DIHEDRALS_WITHOUT_HYDROGEN") != std::string::npos)
        readDihedrals(NumberDihedralsNONH);


    else if (line.find("NUMBER_EXCLUDED_ATOMS") != std::string::npos)
        readNumberExcludedAtomsList();
    else if (line.find("EXCLUDED_ATOMS_LIST") != std::string::npos)
        readExcludedAtomsList();

    else if (line.find("LENNARD_JONES_ACOEF") != std::string::npos)
        readLennardJonesACOEF();
    else if (line.find("LENNARD_JONES_BCOEF") != std::string::npos)
        readLennardJonesBCOEF();

  }
  
  readLennardJonesRVdWEpsilon();

}
catch(std::exception e){
  std::cout << "Any Error ... \n";
  std::cout << e.what() << '\n';
}

	GeneratePairStartAndLen();

  // populate the interval tree for residue names
  residuePointers.push_back(NumberAtoms + 1);
  for (int i = 0; i < residueLabels.size(); i++) {
    residueLabelTree.insert({
      residuePointers[i] - 1,
      residuePointers[i + 1] - 2,
      residueLabels[i]
    });

    residueIndexTree.insert({
      residuePointers[i] - 1,
      residuePointers[i + 1] - 2,
      i + 1
    });
  }
}

void AmberReader::readInpcrd(const std::string& inpcrdfile) {

  try{
    inpcrd.open(inpcrdfile.c_str());
    if(!inpcrd.is_open())
    {
            printf("Error in inpcrd file : %s file not opened \n", inpcrdfile.c_str());
            exit(1);
    }
      // Ignore the first line
    getline(inpcrd, line);

    // Next line contains one or two values
    // The first is the number of atoms present in this file
    // The second, if present, stores the simulation time
    // We ignore both values
    getline(inpcrd, line);

    // Read the coordinates
    SimTK::Real x, y, z;
    while(inpcrd >> x >> y >> z) {
      AtomsXcoord.push_back(x);
      AtomsYcoord.push_back(y);
      AtomsZcoord.push_back(z);
    }

    inpcrd.close();

  }catch(std::exception e){
    std::cout << "Error reading inpcrd file \n";
    std::cout << e.what() << '\n';
  }

}

void AmberReader::readPointers(){
      // draciile de mai jos au un motiv
      getline(prmtop, line);
      for(global_i = 0; global_i <30 ; global_i++) //  30 used terms for sure...we need only some of them
      {
           prmtop >> FlagPointers[global_i];
      }

      NumberAtoms = FlagPointers[0];
      NumberTypes = FlagPointers[1];

      LennardJonesTypes = (NumberTypes * (NumberTypes+1))/2;

      NumberBondsH = FlagPointers[2];
      NumberBondsNONH = FlagPointers[3];
      NumberBonds = NumberBondsH + NumberBondsNONH;

      NumberAnglesH = FlagPointers[4];
      NumberAnglesNONH = FlagPointers[5];
      NumberAngles = NumberAnglesH + NumberAnglesNONH;

      NumberDihedralsH = FlagPointers[6];
      NumberDihedralsNONH = FlagPointers[7];
      NumberDihedrals = NumberDihedralsH + NumberDihedralsNONH;

      NumberExcludedAtoms = FlagPointers[10];
      numres = FlagPointers[11];

      NumberBondsTypes = FlagPointers[15];
      NumberAnglesTypes = FlagPointers[16];
      NumberDihedralsTypes = FlagPointers[17];

  }

  void AmberReader::readResidueLabels() {
    getline(prmtop, line);
    for(global_i = 0; global_i < numres; global_i++)
    {
      std::string label;
      prmtop >> label;
      residueLabels.push_back(label);
    }
  }

  void AmberReader::readResiduePointers() {
    getline(prmtop, line);
    for(global_i = 0; global_i < numres; global_i++)
    {
      std::string label;
      prmtop >> label;
      residuePointers.push_back(std::stoi(label));
    }
  }

  void AmberReader::readAtomsCharge(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomsCharge.push_back(temp_val);
    }
  }

  void AmberReader::readAtomicNumber(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomicNumbers.push_back(temp_val);
    }
  }

  void AmberReader::readAtomsMass(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomsMass.push_back(temp_val);
    }
  }

  void AmberReader::readAtomsRadii(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomsRadii.push_back(temp_val);
    }
  }

  void AmberReader::readAtomsScreenGBIS(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomsScreenGBIS.push_back(temp_val);
    }
  }

  void AmberReader::readAtomsTypesIndex(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAtoms; global_i++)
    {
         prmtop >> temp_val; AtomsTypesIndex.push_back(temp_val);
    }
  }

  void AmberReader::readAtomsName(){
    getline(prmtop, line); // format flag
    global_i=0;
    while(line.find("FLAG") == std::string::npos && global_i < NumberAtoms) // if stays within same field
    {
      getline(prmtop, line);
      for(unsigned int k=0; k + ATOM_NAME_SIZE <=line.length(); k += ATOM_NAME_SIZE )
      {
        auto name = line.substr(k, ATOM_NAME_SIZE );
        auto it = std::remove_if(name.begin(), name.end(), [](char& c) {return std::isspace<char>(c, std::locale::classic()); });
        name.erase(it, name.end());

        AtomsName.push_back(name);
        global_i++;
      }
    }
  }
  void AmberReader::readAtomsTypes(){
    getline(prmtop, line); // format flag
    global_i=0;
    while(line.find("FLAG") == std::string::npos && global_i < NumberAtoms) // if stays within same field
    {
      getline(prmtop, line);
      for(unsigned int k=0; k + ATOM_NAME_SIZE <=line.length(); k += ATOM_NAME_SIZE )
      {
        auto name = line.substr(k, ATOM_NAME_SIZE );
        auto it = std::remove_if(name.begin(), name.end(), [](char& c) {return std::isspace<char>(c, std::locale::classic()); });
        name.erase(it, name.end());

        AtomsTypes.push_back(name);
        global_i++;
      }
    }
  }


  void AmberReader::readNonbondedParmIndex(){
    getline(prmtop, line);
    for(global_i = 0; global_i < (NumberTypes * NumberTypes); global_i++)
    {
         prmtop >> temp_val; LennardJonesNonbondParmIndex.push_back(temp_val);
         if(temp_val < 0){
            printf("Error in prmtop file : Nonobonded 10-12 interactions are not supported");
            exit(1);
         }
    }
  }

  void AmberReader::readTempBondsForceK(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberBondsTypes; global_i++)
    {
         prmtop >>  temp_val; tempBond_K.push_back(temp_val);
    }
  }

  void AmberReader::readTempBondsEqval(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberBondsTypes; global_i++)
    {
         prmtop >>  temp_val; tempBond_eq.push_back(temp_val);
    }
  }

  void AmberReader::readBonds(int nr){
    getline(prmtop, line);
    for(global_i = 0; global_i < nr; global_i++)
    {
      prmtop >> t1;
      prmtop >> t2;
      prmtop >> t3;

      BondsAtomsIndex.push_back(std::vector<int> {
        t1 / AMBER_INDEX_MULTIPLIER , 
        t2 / AMBER_INDEX_MULTIPLIER});
      BondsForceK.push_back(tempBond_K[ t3 - AMBER_INDEX_DIFF]);
      BondsEqval.push_back(tempBond_eq[ t3 - AMBER_INDEX_DIFF]);
    }
  }


  void AmberReader::readTempAnglesForceK(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAnglesTypes; global_i++)
    {
         prmtop >>  temp_val; tempAngle_K.push_back(temp_val);
    }
  }

  void AmberReader::readTempAnglesEqval(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberAnglesTypes; global_i++)
    {
         prmtop >>  temp_val; tempAngle_eq.push_back(temp_val);
    }
  }

  void AmberReader::readAngles(int nr){
    getline(prmtop, line);
    for(global_i = 0; global_i < nr; global_i++)
    {
      prmtop >> t1;
      prmtop >> t2;
      prmtop >> t3;
      prmtop >> t4;

      AnglesAtomsIndex.push_back(std::vector<int> {
        t1 / AMBER_INDEX_MULTIPLIER ,
        t2 / AMBER_INDEX_MULTIPLIER,
        t3 / AMBER_INDEX_MULTIPLIER});
      AnglesForceK.push_back(tempAngle_K[ t4 - AMBER_INDEX_DIFF]);
      AnglesEqval.push_back(tempAngle_eq[ t4 - AMBER_INDEX_DIFF]);
    }
  }




  void AmberReader::readTempDihedralsForceK(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberDihedralsTypes; global_i++)
    {
         prmtop >>  temp_val; tempDihedral_K.push_back(temp_val);
    }
  }

  void AmberReader::readTempDihedralsPhase(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberDihedralsTypes; global_i++)
    {
         prmtop >>  temp_val; tempDihedral_phase.push_back(temp_val);
    }
  }

  void AmberReader::readTempDihedralsPeriod(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberDihedralsTypes; global_i++)
    {
         prmtop >>  temp_val; tempDihedral_per.push_back(temp_val);
    }
  }
  void AmberReader::readTempSCEEScaleFactor(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberDihedralsTypes; global_i++)
    {
         prmtop >>  temp_val; tempDihedral_SCEE.push_back(temp_val);
    }
  }
  void AmberReader::readTempSCNBScaleFactor(){
    getline(prmtop, line);
    for(global_i = 0; global_i < NumberDihedralsTypes; global_i++)
    {
         prmtop >>  temp_val; tempDihedral_SCNB.push_back(temp_val);
    }
  }


  void AmberReader::readDihedrals(int nr){
    getline(prmtop, line);
    for(global_i = 0; global_i < nr; global_i++)
    {
      prmtop >> t1;
      prmtop >> t2;
      prmtop >> t3;
      prmtop >> t4;
      prmtop >> t5;

      DihedralsAtomsIndex.push_back(std::vector<int> {
        std::abs(t1 / AMBER_INDEX_MULTIPLIER), // AmberIndexMultiplier=3
        std::abs(t2 / AMBER_INDEX_MULTIPLIER),
        std::abs(t3 / AMBER_INDEX_MULTIPLIER), 
        std::abs(t4 / AMBER_INDEX_MULTIPLIER)});

      DihedralsForceK.push_back(
        tempDihedral_K[ t5 - AMBER_INDEX_DIFF]); // AmberIndexDiff=1
      DihedralsPeriod.push_back(
        tempDihedral_per[ t5 - AMBER_INDEX_DIFF]);
      DihedralsPhase.push_back(
        tempDihedral_phase[ t5 - AMBER_INDEX_DIFF]);
    }

    /* std::cout << "readAmberInput::readDihedrals " << nr << "\n" ;
    for(int i = 0; i < DihedralsAtomsIndex.size(); i++){

    //   for(int j = 0; j < DihedralsAtomsIndex[i].size(); j++){
    //     std::cout << DihedralsAtomsIndex[i][j] << " ";
    //   }
    
    //   std::cout << DihedralsForceK[i] << " " << DihedralsPeriod[i] << " " << DihedralsPhase[i] << std::endl;

    } */

  }





    void AmberReader::readLennardJonesACOEF(){
      getline(prmtop, line);
      for(global_i = 0; global_i < LennardJonesTypes; global_i++)
      {
           prmtop >>  temp_val; tempLJONES_ACOEFF.push_back(temp_val);
      }
    }
    void AmberReader::readLennardJonesBCOEF(){
      getline(prmtop, line);
      for(global_i = 0; global_i < LennardJonesTypes; global_i++)
      {
           prmtop >>  temp_val; tempLJONES_BCOEFF.push_back(temp_val);
      }
    }

    void AmberReader::readLennardJonesRVdWEpsilon(){
    // adapted from OpenMM - amber_file_parser.py

    SimTK::Real typeRVdW[NumberTypes];
    SimTK::Real typeEpsilon[NumberTypes];

    SimTK::Real rmin;
    SimTK::Real eps;


      for(global_i = 0; global_i < NumberAtoms; global_i++)
      {

          // int index = LennardJonesNonbondParmIndex[ ((NumberTypes + 1)*(AtomsTypesIndex[global_i] - 1)) ] -1;

          int ind = NumberTypes*(AtomsTypesIndex[global_i] - 1) + AtomsTypesIndex[global_i] - 1;
          int index = LennardJonesNonbondParmIndex[ ind ] - 1;

          if( index < 0)
          {
              printf("Error in prmtop file : 10-12 interactions are not supported \n");
              exit(1);
          }

          acoef = tempLJONES_ACOEFF[ index ];
          bcoef = tempLJONES_BCOEFF[ index ];

          if(bcoef != 0 && acoef !=0 ){
            rmin = pow((2.0f*acoef/bcoef), (1.0f/6.0f)) / 2.0f;
            eps = (0.25f*bcoef*bcoef) /acoef ;
          }
          else{
              rmin = 0.5f;
              eps = 0.0f;
          }

          AtomsRVdW.push_back(rmin);
          AtomsEpsilon.push_back(eps);

          typeRVdW[ AtomsTypesIndex[global_i] - 1 ] = rmin;
          typeEpsilon[ AtomsTypesIndex[global_i] - 1 ] = eps;

      }

      // Check if off-diagonal LJ terms
      for(global_i = 0; global_i < NumberTypes; global_i++)
      {
        for(int j = 0; j < NumberTypes; j++)
        {
            int index = LennardJonesNonbondParmIndex[ NumberTypes * global_i + j ] - 1;
            if (index < 0 )
            {
              SimTK::Real rij = typeRVdW[global_i] + typeRVdW[j];
              SimTK::Real wij = sqrt( typeEpsilon[global_i] * typeEpsilon[j] );
              SimTK::Real pairA = tempLJONES_ACOEFF[index];
              SimTK::Real pairB = tempLJONES_BCOEFF[index];

              // Check for 0 terms
              if( ( pairA * pairB ==0 && ( pairA + pairB != 0 || wij * rij !=0 )))
              {
                printf("Error in prmtop file : Off-diagonal terms detected in LJ parameters - null \n");
                exit(1);
              }
              // Check for small terms

              SimTK::Real termA = ( pairA - (wij * pow(rij, 12) ) ) / pairA;
              SimTK::Real termB = ( pairB - ( 2* wij * pow(rij, 6) ) ) / pairB;

              if ( std::abs(termA) > pow(10, -6) || std::abs(termB) > pow(10, -6) )
              {
                printf("Error in prmtop file : Off-diagonal terms detected in LJ parameters - small \n");
                exit(1);
              }
            }
        }
      }
    }

    void AmberReader::readNumberExcludedAtomsList(){
      getline(prmtop, line);
      for(global_i = 0; global_i < NumberAtoms; global_i++)
      {
           prmtop >> temp_val; NumberExcludedAtomsList.push_back(temp_val);
      }
    }

    void AmberReader::readExcludedAtomsList(){
      getline(prmtop, line);


      // create NonBondedAtomsMatrix
      for(global_i = 0; global_i < NumberAtoms; global_i++)
      {
        std::vector<bool> row;
        for( int j = 0; j < NumberAtoms; j++)
        {
          row.push_back(global_i!=j);
        }
        NonBondedAtomsMatrix.push_back(row);
      }

      global_i = 0;
      while ( global_i < NumberAtoms )
      {
           // eliminates pairs involve in bonded interactions
           int count = 1;
           while( count <= NumberExcludedAtomsList[global_i])
           {
             prmtop >> temp_val;

             if(temp_val != 0)
             {
                NonBondedAtomsMatrix[ global_i ][ temp_val - 1 ] = 0;
                NonBondedAtomsMatrix[ temp_val - 1 ][ global_i ] = 0;
             }

             count++;
           }
           global_i++;
      }
    }






// GET FUNCTIONS


int AmberReader::getNumberAtoms() const {
  return NumberAtoms;
}

int AmberReader::getNumberBonds() const {
  return NumberBonds;
}

int AmberReader::getNumberAngles() const {
  return NumberAngles;
}

int AmberReader::getNumberDihedrals() const {
  return NumberDihedrals;
}


const std::string& AmberReader::getAtomsName(int p) const {
  return AtomsName[p];
}

const std::string& AmberReader::getAtomsType(int p) const {
  return AtomsTypes[p];
}


SimTK::Real AmberReader::getAtomsXcoord(int p) const {
  return AtomsXcoord[p];
}
SimTK::Real AmberReader::getAtomsYcoord(int p) const {
  return AtomsYcoord[p];
}
SimTK::Real AmberReader::getAtomsZcoord(int p) const {
  return AtomsZcoord[p];
}
SimTK::Real AmberReader::getAtomsMass(int p) const {
  return AtomsMass[p];
}
SimTK::Real AmberReader::getAtomsCharge(int p) const {
  return AtomsCharge[p];
}
SimTK::Real AmberReader::getAtomsRadii(int p) const {
  return AtomsRadii[p];
}
SimTK::Real AmberReader::getAtomsRVdW(int p) const {
  return AtomsRVdW[p];
}
SimTK::Real AmberReader::getAtomsEpsilon(int p) const {
  return AtomsEpsilon[p];
}

SimTK::Real AmberReader::getAtomicNumber(int p) const {
  return AtomicNumbers[p];
}

int AmberReader::getBondsAtomsIndex1(int bond) const {
  return BondsAtomsIndex[bond][0];
}
int AmberReader::getBondsAtomsIndex2(int bond) const {
  return BondsAtomsIndex[bond][1];
}
SimTK::Real AmberReader::getBondsForceK(int bond) const {
  return BondsForceK[bond];
}
SimTK::Real AmberReader::getBondsEqval(int bond) const {
  return BondsEqval[bond];
}




int AmberReader::getAnglesAtomsIndex1(int angle) const {
  return AnglesAtomsIndex[angle][0];
}
int AmberReader::getAnglesAtomsIndex2(int angle) const {
  return AnglesAtomsIndex[angle][1];
}
int AmberReader::getAnglesAtomsIndex3(int angle) const {
  return AnglesAtomsIndex[angle][2];
}
SimTK::Real AmberReader::getAnglesForceK(int angle) const {
  return AnglesForceK[angle];
}
SimTK::Real AmberReader::getAnglesEqval(int angle) const {
  return AnglesEqval[angle];
}



int AmberReader::getDihedralsAtomsIndex1(int dihedral) const {
  return DihedralsAtomsIndex[dihedral][0];
}
int AmberReader::getDihedralsAtomsIndex2(int dihedral) const {
  return DihedralsAtomsIndex[dihedral][1];
}
int AmberReader::getDihedralsAtomsIndex3(int dihedral) const {
  return DihedralsAtomsIndex[dihedral][2];
}
int AmberReader::getDihedralsAtomsIndex4(int dihedral) const {
  return DihedralsAtomsIndex[dihedral][3];
}

int AmberReader::getDihedralsAtomsIndex(int dihIndex, int atomIndx)
{
	assert( atomIndx > 0 && atomIndx < 4 );
	return DihedralsAtomsIndex[ dihIndex ][ atomIndx ];
}

//HOREA
/* void readAmberInput::GeneratePairStartAndLen()
{
  if( DihedralsAtomsIndex.size() ){
    std::vector<int> currentDihedralIndices = DihedralsAtomsIndex[0];
    PairStartAndLen.push_back( std::make_pair(0,1) );

    for(unsigned int idx = 1 ; idx < DihedralsAtomsIndex.size() ; idx++ )
    {
      if( currentDihedralIndices[0] == DihedralsAtomsIndex[idx][0] &&
          currentDihedralIndices[1] == DihedralsAtomsIndex[idx][1] &&
          currentDihedralIndices[2] == DihedralsAtomsIndex[idx][2] &&
          currentDihedralIndices[3] == DihedralsAtomsIndex[idx][3] )
      {
        PairStartAndLen[ PairStartAndLen.size() -1 ].second ++;
      }
      else {
        int nextIndx = PairStartAndLen.back().first;
        PairStartAndLen.push_back( std::make_pair( nextIndx + 1, 1) );
      }
    }
  }
	
} */

// Laurentiu
void AmberReader::GeneratePairStartAndLen()
{
  if( DihedralsAtomsIndex.size() ){

    std::vector<int> currentDihedralIndices = DihedralsAtomsIndex[0];
    PairStartAndLen.push_back( std::make_pair(0,0) );
    int idx = 0;

    while(idx < DihedralsAtomsIndex.size())
    {

      if( currentDihedralIndices[0] == DihedralsAtomsIndex[idx][0] &&
          currentDihedralIndices[1] == DihedralsAtomsIndex[idx][1] &&
          currentDihedralIndices[2] == DihedralsAtomsIndex[idx][2] &&
          currentDihedralIndices[3] == DihedralsAtomsIndex[idx][3] )
      {
        PairStartAndLen[ PairStartAndLen.size() - 1 ].second ++;
      }
      else {
        //int nextIndx = PairStartAndLen.back().first;
        int nextIndx = idx;
        PairStartAndLen.push_back( std::make_pair( nextIndx, 1) );

        currentDihedralIndices[0] = DihedralsAtomsIndex[idx][0];
        currentDihedralIndices[1] = DihedralsAtomsIndex[idx][1];
        currentDihedralIndices[2] = DihedralsAtomsIndex[idx][2];
        currentDihedralIndices[3] = DihedralsAtomsIndex[idx][3];

      }

/*
      std::cout << "GeneratePairStartAndLen " 
        << " " << DihedralsAtomsIndex[idx][0] << " " << DihedralsAtomsIndex[idx][1]
        << " " << DihedralsAtomsIndex[idx][2] << " " << DihedralsAtomsIndex[idx][3]
        << " || " << PairStartAndLen[ PairStartAndLen.size() - 1 ].first
        << " " << PairStartAndLen[ PairStartAndLen.size() - 1 ].second
        << std::endl << std::flush;
*/
      idx++;

    }
  }
	
}

std::vector < std::pair<int, int> > AmberReader::getPairStartAndLen() const
{
	return PairStartAndLen;
}


SimTK::Real AmberReader::getDihedralsForceK(int dihedral) const {
  return DihedralsForceK[dihedral];
}
SimTK::Real AmberReader::getDihedralsPhase(int dihedral) const {
  return DihedralsPhase[dihedral];
}
SimTK::Real AmberReader::getDihedralsPeriod(int dihedral) const {
  return DihedralsPeriod[dihedral];
}

bool AmberReader::getNonBondedAtomsMatrix(int at1, int at2) const {
  return NonBondedAtomsMatrix[at1][at2];
}

const std::string& AmberReader::getResidueLabel(int i) const {
  return residueLabelTree.get(i);
}

int AmberReader::getResidueIndex(int i) const {
  return residueIndexTree.get(i);
}
