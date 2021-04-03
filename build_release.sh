#!/bin/bash

# notes for the future
# make -j$((`nproc`*2)+1) -> c++: fatal error: Killed signal terminated program cc1plus (i.e too much memory consumed)

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
mkdir build-release
cd build-release/
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install
cd ../../

# delete previous compilation from project (Simbody)
cd Molmodel/Simbody01
rm build-debug -rf
rm build-release -rf

# remove Molmodel
cd ../
rm build-debug -rf
rm build-release -rf

# remove Robosample
cd ../
rm build-debug -rf
rm build-release -rf

# compile project (Simbody)
cd Molmodel/Simbody01
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install

# compile Molmodel
cd ../../
mkdir build-release
cd build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$((`nproc`*2))
sudo make install

# ensure this is correctly installed
sudo mkdir -p /usr/local/lib/plugins/ && sudo cp -f libOpenMMPlugin.so $_
rm /usr/local/lib/plugins/libOpenMMPlugin_d.so

# compile Robosample
cd ../../
mkdir build-release
cd build-release
mkdir pgo
cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Generate -DROBO_PGO_PATH="$(pwd)/pgo/" ..
make -j$((`nproc`*2))
sudo /sbin/ldconfig

# add test input files
cp -ri ../tests_inputs/* .
mkdir temp
mkdir temp/pdbs

# run PGO on Robosample and recompile
./tests/Robosample inp.dock.test
cmake -DCMAKE_BUILD_TYPE=Release -DROBO_PGO=Use -DROBO_PGO_PATH="$(pwd)/pgo/" ..
make -j$((`nproc`*2))
sudo /sbin/ldconfig

# add tests
cp ../tests/test-memory.sh test-memory.sh


# ###
# mkdir build-release
# cd build-release
# mkdir pgo
# cmake -DCMAKE_BUILD_TYPE=Release ..
# make -j$((`nproc`*2))
# cp -ri ../tests_inputs/* .
# mkdir temp
# mkdir temp/pdbs
# sudo /sbin/ldconfig
# ./tests/Robosample inp
# ###

# Avem date care arata bine. So what?
# in prezentare, pot sa fie cine vreau (extrovertit, actor, inginer etc).
# e vina mea daca merge ceva prost in prezentare - asa zice audienta (laptop, proiector stricat etc)
# cei care spun glume bune in prezentare au glumele respective pregatite
# 99% din carti folosesc 7 scenarii de baza (vezi poza)
# hero defeating a monster: firma mea face mult mai mult bine fata de competitori (praslea cel voinic, terminator, filmele marvel) - frecvent fol in industrie
# rags to riches: metoda mea taie costurile (alladin, punguta cu doi bani)
# quest for a treasure: un grup de pers (sus era doar una): chemare init, calatoria, frustrare la sfarsit, incerca ceva -> final fericit (vrajitorul in oz, indiana jones) - rar fol in afaceri
# prezentare ft buna: richard turere - my invetion that made peace with the lions; info spuse complementare fata de slide (deci nu repeta)
# folosesti o gama larga de stimuli: explic sunetul, mirosul etc -> activezi imaginatia audientei
# leslie morgen steiner domestic violence
# poveste personala: cum a reusit op sa salveze cu tehnologia gemenele cu care era gravida sotia lui
# spune o poveste predictibila, iar apoi surprinde audienta
# star moment: smth people will remember
# bill gates vb despre malarie: a adus un borcan cu tantari si l-a desfacut pe scena (nu doar saracii trebuie sa sufere)
# info trebuie sa fie relateable, sa inteleaga oricine
# data: start with eyecatcher, be simple
# do not use 3d maps (unless you are cool, but it does not work neaparat)
# the joy of stats 200 countries 200 years in 4 minutes