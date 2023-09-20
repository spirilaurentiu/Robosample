#include "Context.hpp"

struct DISTANCE {
    int i = -1,
        j = -1;
    SimTK::Real distance = -1;
};

constexpr std::array<DISTANCE, 14> distances = {{
    { 0, 1, 0.97 },
    { 0, 2, 1.42 },
    { 2, 11, 1.52 },
    { 2, 3, 1.1 },
    { 11, 12, 1.1 },
    { 11, 13, 1.1 },
    { 11, 14, 1.1 },
    { 2, 4, 1.53 },
    { 4, 5, 1.1 },
    { 4, 6, 1.1 },
    { 4, 7, 1.52 },
    { 7, 8, 1.09 },
    { 7, 9, 1.09 },
    { 7, 10, 1.09 },
}};

int main() {
    Context c;
    if (!c.initializeFromFile("inp.2but")) {
        std::cout << "Failed to initialize context from inp.2but" << std::endl;
        return 1;
    }

    for (std::size_t i = 0; i < distances.size(); i++) {
        const auto d = c.Distance(0, 0, 0, distances[i].i, distances[i].j) * 10;
        std::cout << distances[i].i << " " << distances[i].j << " " << distances[i].distance << " " << d << std::endl;
    }

    // c.getWorld(0)->topologies[0]->

    return 0;
}