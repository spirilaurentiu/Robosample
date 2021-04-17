#!/bin/bash

<<<<<<< Updated upstream
# notes for the future
# make -j$((`nproc`*2)+1) -> c++: fatal error: Killed signal terminated program cc1plus (i.e too much memory consumed)
# make -j$(`nproc`) is too slow on my machine

# remove openmm previous openmm compilation
sudo rm /usr/local/openmm -rf
sudo find /usr/local/lib -iname '*openmm*' -delete
rm openmm/build-debug -rf
rm openmm/build-release -rf

# delete previous simbody compilation from /usr/
sudo rm /usr/local/lib/cmake/simbody -rf
sudo rm /usr/local/lib/simbody -rf
sudo rm /usr/local/lib/pkgconfig/simbody.pc
sudo find /usr/local/lib -name 'libSimTK*' -delete

# compile openmm
cd openmm
git checkout master && git pull
mkdir build-release
cd build-release/
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install
cd ../../
=======
# # notes for the future
# # make -j$((`nproc`*2)+1) -> c++: fatal error: Killed signal terminated program cc1plus (i.e too much memory consumed)
# # make -j$(`nproc`) is too slow on my machine
>>>>>>> Stashed changes

# # remove openmm previous openmm compilation
# sudo rm /usr/local/openmm -rf
# sudo find /usr/local/lib -iname '*openmm*' -delete
# rm openmm/build-debug -rf
# rm openmm/build-release -rf

# # delete previous simbody compilation from /usr/
# sudo rm /usr/local/lib/cmake/simbody -rf
# sudo rm /usr/local/lib/simbody -rf
# sudo rm /usr/local/lib/pkgconfig/simbody.pc
# sudo find /usr/local/lib -name 'libSimTK*' -delete

# # compile openmm
# cd openmm
# mkdir build-release
# cd build-release/
# cmake -DCMAKE_BUILD_TYPE=Release ..
# make -j$((`nproc`*2))
# sudo make install
# cd ../../

<<<<<<< Updated upstream
# compile project (Simbody)
cd Molmodel/Simbody01
git checkout master && git pull
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install

# compile Molmodel
cd ../../
git checkout master && git pull
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install
=======
# # delete previous compilation from project (Simbody)
# cd Molmodel/Simbody01
# rm build-debug -rf
# rm build-release -rf

# # remove Molmodel
# cd ../
# rm build-debug -rf
# rm build-release -rf
>>>>>>> Stashed changes

# # remove Robosample
# cd ../
# rm build-debug -rf
# rm build-release -rf

<<<<<<< Updated upstream
# compile Robosample
cd ../../
git checkout master && git pull
=======
# # compile project (Simbody)
# cd Molmodel/Simbody01
# mkdir build-release
# cd build-release
# cmake -DCMAKE_BUILD_TYPE=Release ..
# make -j$((`nproc`*2))
# sudo make install

# # compile Molmodel
# cd ../../
# mkdir build-release
# cd build-release
# cmake -DCMAKE_BUILD_TYPE=Release ..
# make -j$((`nproc`*2))
# sudo make install

# # ensure this is correctly installed
# sudo mkdir -p /usr/local/lib/plugins/ && sudo cp -f libOpenMMPlugin.so $_
# rm /usr/local/lib/plugins/libOpenMMPlugin_d.so

# # compile Robosample
# cd ../../
rm build-release -rf
>>>>>>> Stashed changes
mkdir build-release
cd build-release
mkdir pgo
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2)) 2> out.txt
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# # compile Robosample
# cd ../../
# mkdir build-release
# cd build-release
# mkdir pgo
# cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Generate -DROBO_PGO_PATH="$(pwd)/pgo/" ..
# make -j$((`nproc`*2))
# sudo /sbin/ldconfig
# # add test input files
# cp -ri ../tests_inputs/* .
# mkdir temp
# mkdir temp/pdbs
# # run PGO on Robosample and recompile
# ./tests/Robosample inp.dock.test
# cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Use -DROBO_PGO_PATH="$(pwd)/pgo/" ..
# make -j$((`nproc`*2))
# sudo /sbin/ldconfig

# add tests
cp ../tests/test-memory.sh test-memory.sh


World 1, NU 287:
       pe_o 13846.48659, pe_n nan, pe_nB nan
       ke_prop 810.30598, ke_n -nan
       fix_o 8842.72392, fix_n -nan
       logSineSqrGamma2_o -0.57732, logSineSqrGamma2_n -0.57732
       ts 0.00110, exp(bdE) -nan
       etot_n nan, etot_proposed 23501.01653
       MSD= 0.0000000000, RRdot= -nan
CONTACT INFO: #forces= 0 dissEnergy= 0.0000000000 hasDefaultForceGenerator= 1 #mobods= 289 ctsNofSurfaces= 0
Parameter 3 to routine DGEBAL was incorrect
Parameter 2 to routine DGEHRD  was incorrect
Parameter 2 to routine DORGHR DORGQR was incorrect
Parameter 4 to routine DHSEQR was incorrect
terminate called after throwing an instance of 'SimTK::Exception::IllegalLapackArg'
 what():  SimTK Exception thrown at LapackInterface.cpp:598:
 SimTK internal error: dgeev called with an illegal value to argument #4.
Please report this at SimTK.org.
de obicei cand primeam nan imi respingea sample-ul
acum nu mi-l mai respinge, imi da eroare
daca ai vreo idee, please help
o sa ma uit si eu 
MHAcceptProbability 
trebuie sa puna si conditia cu nan-uri cum era inainte
cand ai timp, da-mi de veste poate vorbim online
