//
// Created by Lars Gebraad on 17.03.18.
//

#include <omp.h>
#include "fwiExperiment.h"
#include "fwiPropagator.h"

using namespace arma;

// Constructor
fwiExperiment::fwiExperiment(imat _receivers, imat _sources, vec _sourceFunction,
                             double _samplingTime, double _samplingTimestep, int _samplingAmount) {
    // Create a fwiExperiment
    receivers = std::move(_receivers);
    sources = std::move(_sources);
    sourceFunction = std::move(_sourceFunction);
    samplingTime = _samplingTime;
    samplingTimestep = _samplingTimestep;
    samplingAmount = _samplingAmount;
    currentModel.setTime(samplingTimestep, samplingAmount, samplingTime);

    muKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    densityKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    lambdaKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);

    // Check for positions
    for (auto &&yPosReceiver : receivers.col(1)) {
        if (yPosReceiver >= static_cast<int>(currentModel.nz_domain)) {
            throw std::invalid_argument("Invalid y position for receiver (in or beyond the Gaussian taper).");
        }
    }
    for (auto &&yPosSource : sources.col(1)) {
        if (yPosSource >= static_cast<int>(currentModel.nz_domain)) {
            throw std::invalid_argument("Invalid y position for receiver (in or beyond the Gaussian taper).");
        }
    }

    // This way all shots have the same receivers and sources
    for (uword ishot = 0; ishot < sources.n_rows; ++ishot) {
        shots.emplace_back(
                fwiShot(sources.row(ishot), receivers, sourceFunction, samplingAmount, samplingTimestep, samplingTime,
                        ishot, snapshotInterval));
    }

}

void fwiExperiment::forwardData() {
    for (uword iShot = 0; iShot < sources.n_rows; ++iShot) {
        fwiPropagator::propagateForward(currentModel, shots[iShot]);
    }
}

void fwiExperiment::writeShots(file_type type, std::string &_folder) {
    for (auto &&shot : shots) {
        shot.writeShot(type, _folder);
    }
}

void fwiExperiment::computeKernel() {
    calculateAdjointSourcesL2();
    muKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    densityKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    lambdaKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    backwardAdjoint();
}

void fwiExperiment::calculateMisfit() {
    misfit = 0;
    for (auto &&shot : shots) {
        misfit += 0.5 * shot.samplingTimestep * accu(square(shot.seismogramObs_ux - shot.seismogramSyn_ux));
        misfit += 0.5 * shot.samplingTimestep * accu(square(shot.seismogramObs_uz - shot.seismogramSyn_uz));
    }
}

void fwiExperiment::calculateAdjointSourcesL2() {
    for (auto &&shot : shots) {
        shot.calculateAdjointSources();
    }
}

void fwiExperiment::backwardAdjoint() {
    for (uword iShot = 0; iShot < sources.n_rows; ++iShot) {
        fwiPropagator::propagateAdjoint(currentModel, shots[iShot], densityKernel_par1, muKernel_par1, lambdaKernel_par1);
    }
}

void fwiExperiment::loadShots(std::string &_folder) {
    for (auto &&shot : shots) {
        shot.loadShot(_folder);
    }
}

fwiExperiment::fwiExperiment() {
    receivers = imat();
    sources = imat();
    sourceFunction = vec();
    muKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    densityKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    lambdaKernel_par1 = zeros(currentModel.nx_domain, currentModel.nz_domain);
    shots = std::vector<fwiShot>();
}

void fwiExperiment::mapKernels() {
    densityKernel_par2 = densityKernel_par1 + (square(currentModel.vp) - 2 * square(currentModel.vs)) % lambdaKernel_par1
            + square(currentModel.vs) % muKernel_par1;

    vpKernel_par2 = 2 * currentModel.vp % lambdaKernel_par1 / currentModel.b_vx;

    // TODO validate next line, the vs * Klambda is a bit weird
    vsKernel_par2 = 2 * currentModel.vs % muKernel_par1 / currentModel.b_vx; - 4 * currentModel.vs % lambdaKernel_par1 / currentModel.b_vx;
}