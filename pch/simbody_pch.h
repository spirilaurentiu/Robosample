#pragma once

#include <algorithm>
// #include <altivec.h>
#include <array>
#include <atomic>
// #include <Availability.h>
#include <bitset>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <chrono>
#include <climits>
#include <cmath>
#include <complex>
#include <condition_variable>
// #include <cpodes/cpodes_band.h>
// #include <cpodes/cpodes_bandpre.h>
// #include <cpodes/cpodes_bbdpre.h>
// #include <cpodes/cpodes_dense.h>
// #include <cpodes/cpodes_direct.h>
// #include <cpodes/cpodes.h>
// #include <cpodes/cpodes_lapack_exports.h>
// #include <cpodes/cpodes_lapack.h>
// #include <cpodes/cpodes_spbcgs.h>
// #include <cpodes/cpodes_spgmr.h>
// #include <cpodes/cpodes_spils.h>
// #include <cpodes/cpodes_sptfqmr.h>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cxxabi.h>
#include <deque>
// #include <direct.h>
#include <dirent.h>
#include <dlfcn.h>
#include <emmintrin.h>
#include <exception>
#include <fcntl.h>
#include <float.h>
#include <fstream>
#include <future>
// #include <GL/glext.h>
// #include <GL/gl.h>
// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glu.h>
// #include <GL/glut.h>
// #include <GLUT/glut.h>
#include <initializer_list>
#include <inttypes.h>
// #include <io.h>
#include <iomanip>
// #include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
// #include <mach/mach.h>
// #include <mach/mach_time.h>
// #include <mach-o/dyld.h>
#include <map>
#include <math.h>
#include <memory>
// #include <mpi.h>
#include <mutex>
// #include <nvector/nvector_serial.h>
#include <ostream>
// #include <process.h>
#include <queue>
#include <regex>
#include <set>
#include <sstream>
#include <stdarg.h>
#include <stddef.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
// #include <sundials/sundials_band.h>
// #include <sundials/sundials_config.h>
// #include <sundials/sundials_dense.h>
// #include <sundials/sundials_direct.h>
// #include <sundials/sundials_fnvector.h>
// #include <sundials/sundials_iterative.h>
// #include <sundials/sundials_lapack.h>
// #include <sundials/sundials_math.h>
// #include <sundials/sundials_nvector.h>
// #include <sundials/sundials_spbcgs.h>
// #include <sundials/sundials_spgmr.h>
// #include <sundials/sundials_sptfqmr.h>
// #include <sundials/sundials_types.h>
// #include <sys/stat.h>
// #include <sys/sysctl.h>
// #include <sys/time.h>
#include <thread>
#include <time.h>
#include <typeinfo>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>
// #include <windows.h>




// #include "SimTKlapack.h"


#include "SimTKcommon/basics.h"
#include "SimTKcommon/Constants.h"
// #include "SimTKcommon/SmallMatrix.h"
// #include "SimTKcommon/TemplatizedLapack.h"
// #include "SimTKcommon/Testing.h"


#include "SimTKcommon/Scalar.h"


#include "SimTKcommon/internal/common.h"
// #include "SimTKcommon/Orientation.h"
// #include "SimTKcommon/Mechanics.h"
