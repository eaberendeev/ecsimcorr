// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>

#include "simulation.h"

// Main function simply hands off control to the Simulation class
int main(int argc, char **argv) {
    
    Simulation simulation(argc, argv);

    simulation.make_all();

    return 0;
}