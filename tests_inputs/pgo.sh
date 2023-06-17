./robosample.pgo.train inp.2but
mv profile-data/default.profraw profile-data/2but.profraw
# mv profile-data/default*.profraw profile-data/2but.profraw

./robosample.pgo.train inp.ala10
mv profile-data/default.profraw profile-data/ala10.profraw
# mv profile-data/default*.profraw profile-data/ala10.profraw

./robosample.pgo.train inp.aper
mv profile-data/default.profraw profile-data/aper.profraw
# mv profile-data/default*.profraw profile-data/aper.profraw

./robosample.pgo.train inp.e2un
mv profile-data/default.profraw profile-data/e2un.profraw
# mv profile-data/default*.profraw profile-data/e2un.profraw

./robosample.pgo.train inp.mla
mv profile-data/default.profraw profile-data/mla.profraw
# mv profile-data/default*.profraw profile-data/mla.profraw

# does not work with pgo
# ./robosample.pgo.train inp.2mols
# mv profile-data/default.profraw profile-data/2mols.profraw

llvm-profdata-14 merge -o profile-data/default.profdata $(ls profile-data/*.profraw)