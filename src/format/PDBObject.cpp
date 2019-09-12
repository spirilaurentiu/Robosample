#include "format/PDBObject.hpp"

// Uses standard PDB formatting reccomandations
void PDBObject::readPDBfile(std::string pdbfile)
{

try
{

  pdbinput.open(pdbfile.c_str());
  if(!pdbinput.is_open())
  {
          printf("Error in reading PDB file : file not opened \n");
          exit(1);
  }

  NumberAtoms = 0;
  int lineno = 0;

  while(!pdbinput.eof())
  {
      getline(pdbinput, line);
      lineno++;
      if(line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
      {
        try
        {
          AtomsName.push_back( line.substr(12, 4) );
          AtomsMolName.push_back( line.substr(17, 3) );
          AtomsChain.push_back( line.substr(21, 1) );
          AtomsMolNo.push_back( std::stoi( line.substr(22, 4) ) );

          AtomsXcoord.push_back( std::stod( line.substr(30, 8) ) );
          AtomsYcoord.push_back( std::stod( line.substr(38, 8) ) );
          AtomsZcoord.push_back( std::stod( line.substr(46, 8) ) );

          std::string temp = line.substr(12, 4);
          for(int it=0; it<4; it++)
          {
            if (temp.substr(it, 1) != " ")
            {
              Atoms1Character.push_back( temp.substr(it, 1) );
              break;
            }
          }

          NumberAtoms++;
        }
        catch(std::exception e){
          std::cout << "Error in reading PDB input file : ATOM/HETATM section data format at line: " << lineno << "\n" ;
          std::cout << e.what() << '\n';
        }
      }
  }

  pdbinput.close();


}
catch(std::exception e){
  std::cout << "Error in reading PDB file - Unknown Error\n";
  std::cout << e.what() << '\n';
}

}


// Uses standard PDB formatting reccomandations
// TO DO: occupancy and temperature factor
void PDBObject::writePDBfile(std::string pdbfile,
                              std::vector<std::string> AtomsName,
                              std::vector<std::string> AtomsMolName,
                              std::vector<std::string> AtomsChain,
                              std::vector<int> AtomsMolNo,
                              std::vector<TARGET_TYPE> AtomsXcoord,
                              std::vector<TARGET_TYPE> AtomsYcoord,
                              std::vector<TARGET_TYPE> AtomsZcoord,
                              std::vector<std::string> Atoms1Character)
{

try
{
  FILE * pdboutput;
  pdboutput = fopen(pdbfile.c_str(), "w");

  for(int it = 0; it < NumberAtoms; it++)
  {
    fprintf(pdboutput, "ATOM  ");
    fprintf(pdboutput, "%5i ", it);
    fprintf(pdboutput, "%4s ", AtomsName[it].c_str() );
    fprintf(pdboutput, "%3s ", AtomsMolName[it].c_str() );
    fprintf(pdboutput, "%1s", AtomsChain[it].c_str() );
    fprintf(pdboutput, "%4i    ", AtomsMolNo[it] );
    fprintf(pdboutput, "%8.3f", AtomsXcoord[it] );
    fprintf(pdboutput, "%8.3f", AtomsYcoord[it] );
    fprintf(pdboutput, "%8.3f", AtomsZcoord[it] );

    // No occupancy or temperature factor
    fprintf(pdboutput, "%6.2f", 1.0 );
    fprintf(pdboutput, "%6.2f", 0.0 );

    fprintf(pdboutput, "      %1s", AtomsChain[it].c_str());
    fprintf(pdboutput, "    %1s", Atoms1Character[it].c_str());

    fprintf(pdboutput, "\n");

  }


  fprintf(pdboutput, "END");
  fclose(pdboutput);

}
catch(std::exception e){
  std::cout << "Error in writing PDB file \n";
  std::cout << e.what() << '\n';
}

}


// GET FUNCTIONS

int PDBObject::getNumberAtoms(){
  return NumberAtoms;
}
int PDBObject::getAtomsMolNo(int p){
  return AtomsMolNo[p];
}


std::string PDBObject::getAtomsName(int p){
  return AtomsName[p];
}
std::string PDBObject::getAtomsMolName(int p){
  return AtomsMolName[p];
}
std::string PDBObject::getAtomsChain(int p){
  return AtomsChain[p];
}
std::string PDBObject::getAtoms1Character(int p){
  return Atoms1Character[p];
}




TARGET_TYPE PDBObject::getAtomsXcoord(int p){
  return AtomsXcoord[p];
}
TARGET_TYPE PDBObject::getAtomsYcoord(int p){
  return AtomsYcoord[p];
}
TARGET_TYPE PDBObject::getAtomsZcoord(int p){
  return AtomsZcoord[p];
}
