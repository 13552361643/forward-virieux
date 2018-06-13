//
// Created by Lars Gebraad on 30/05/18.
//

#include <iostream>
#include <armadillo>

#define _USE_MATH_DEFINES

#include <cmath>

// Own includes
#include "forward-psv/fwiExperiment.h"
#include "forward-psv/fwiPropagator.h"
#include "forward-psv/fwiShot.h"
#include "misc/functions.h"

using namespace arma;
using namespace std;
typedef vector<double> stdvec;

int main() {

    // Loading sources and receivers
    string experimentFolder = "Original600x600model";
    string receiversFile = experimentFolder + string("/receivers.txt");
    string sourcesFile = experimentFolder + string("/sources.txt");
    string dimensionFile = experimentFolder + string("/dimensions.txt");
    imat receivers;
    imat sources;
    vec dimensions;
    receivers.load(receiversFile);
    sources.load(sourcesFile);
    dimensions.load(dimensionFile);

    // Load dimensions
    uword nx = static_cast<uword>(dimensions[0]);
    double dx = dimensions[1];
    uword nz = static_cast<uword>(dimensions[2]);
    double dz = dimensions[3];
    uword np = static_cast<uword>(dimensions[4]);
    double np_f = dimensions[5];

    // Create stf
    vec sourcefunction;
    double dt = dimensions[6];
    unsigned int nt = static_cast<unsigned int>(dimensions[7]);
    double freq = dimensions[8];
    sourcefunction = generateRicker(dt, nt, freq);

    // Loading set up
    fwiExperiment experiment(dx, dz, nx, nz, np, np_f, receivers, sources, sourcefunction, dt, nt, fwiShot::momentSource);

    // Loading material parameters
    mat rho;
    mat c11;
    mat c13;
    mat c33;
    mat c55;

    rho.load(experimentFolder + string("/M.txt"));
    c11.load(experimentFolder + string("/C1111.txt"));
    c13.load(experimentFolder + string("/C1122.txt"));
    c33.load(experimentFolder + string("/C2222.txt"));
    c55.load(experimentFolder + string("/C1212.txt"));

    std::cout <<  sqrt(max(max(c11))/ min(min(rho)));

//    return EXIT_SUCCESS;

    experiment.updateAnisotropicElasticity(rho.t(), c11.t(), c13.t(), c33.t(), c55.t());

    experiment.exportSnapshots = true;
//    experiment.useRamSnapshots = true;
    experiment.performFWI = false;

    for (unsigned int i = 0; i < nt; i += 10) {
        experiment.snapshots.emplace_back(i);
    }

    experiment.forwardData();

    return EXIT_SUCCESS;
}