#ifndef __ROBO_HPP
#define __ROBO_HPP

#include "Simbody.h"
#include "Molmodel.h"

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//#include <Eigen/LU>
//#include <Eigen/QR>

#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cctype>
#include <deque>
#include <memory>
#include <thread>
#include <bitset>

#include <random> // FOR RANDOM BLOCKS
#include <functional>

// bits/stdc++.h (clang compatibility)
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#include "bgeneral.hpp"
#include "readAmberInput.hpp"
#include "SetupReader.hpp"
#include "trim.hpp"

//#ifndef DEBUG_ROBO
//#define DEBUG_ROBO
//#endif

#ifdef DEBUG_ROBO
#define TRACE(STR) printf("%s", STR);
#else
#define TRACE(STR)
#endif


#endif // __ROBO_HPP
