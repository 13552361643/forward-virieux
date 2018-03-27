//
// Created by lars on 17.03.18.
//
#include <armadillo>
#include "model.h"
#include "propagator.h"
#include "experiment.h"

void propagator::propagateForward(model &_currentModel, shot &_shot, bool storeWavefieldBoundary) {

    // Some standard output
    std::cout << "    Stability number: " <<
              sqrt(_currentModel.lm.max() * _currentModel.b_vx.max()) * _shot.dt *
              sqrt(1.0 / (_currentModel.dx * _currentModel.dx) + 1.0 / (_currentModel.dz * _currentModel.dz))
              << std::endl;

    // Loading simulation parameters
    double dx = _currentModel.dx;
    const arma::uword nx = _currentModel.nx;
    double dz = _currentModel.dz;
    const arma::uword nz = _currentModel.nz;

    // Create dynamic fields
    arma::mat vx = arma::zeros(_currentModel.nx, _currentModel.nz);
    arma::mat vz = arma::zeros(_currentModel.nx, _currentModel.nz);
    arma::mat txx = arma::zeros(_currentModel.nx, _currentModel.nz);
    arma::mat tzz = arma::zeros(_currentModel.nx, _currentModel.nz);
    arma::mat txz = arma::zeros(_currentModel.nx, _currentModel.nz);

    // Create taper matrix
    arma::mat taper = _currentModel.np_boundary * arma::ones(nx, nz);
    for (int iTaper = 0; iTaper < _currentModel.np_boundary; ++iTaper) {
        taper.submat(iTaper, 0, nx - iTaper - 1, nz - iTaper - 1) =
                1 + iTaper * arma::ones(nx - 2 * iTaper, nz - iTaper);
    }
    taper = arma::exp(-arma::square(_currentModel.np_factor * (_currentModel.np_boundary - taper)));

    // Create indices for injection
    auto horSpanX = arma::span(_currentModel.np_boundary, nx - _currentModel.np_boundary - 1);
    auto topSpanY = arma::span(1, 1);
    auto bottomSpanY = arma::span(_currentModel.nz_domain - 1, _currentModel.nz_domain - 1);
    auto verSpanY = arma::span(0 + 1, _currentModel.nz_domain - 2);
    auto leftSpanX = arma::span(_currentModel.np_boundary, _currentModel.np_boundary);
    auto rightSpanX = arma::span(nx - _currentModel.np_boundary - 1, nx - _currentModel.np_boundary - 1);

    // Time marching through all time levels
    for (int it = 0; it < _shot.nt; ++it) {
        // Inject source
        for (int source = 0; source < _shot.source.n_rows; ++source) {
            int ix = _shot.source.row(source)[0];
            int iz = _shot.source.row(source)[1];

            txx(ix + _currentModel.np_boundary, iz) += 0.5 * _shot.dt * _shot.sourceFunction[it];
            tzz(ix + _currentModel.np_boundary, iz) += 0.5 * _shot.dt * _shot.sourceFunction[it];
        }

        // Record wavefield at receivers
        for (int receiver = 0; receiver < _shot.receivers.n_rows; ++receiver) {
            int ix = _shot.receivers.row(receiver)[0];
            int iz = _shot.receivers.row(receiver)[1];

            _shot.forwardData_vx(receiver, it) = vx(ix + _currentModel.np_boundary, iz);
            _shot.forwardData_vx(receiver, it) = vz(ix + _currentModel.np_boundary, iz);
        }

        bool pragmaMethod = true;
        if (storeWavefieldBoundary) {
            if (pragmaMethod) {
                _shot.boundaryRecVxTop.col(it) = vx(horSpanX, topSpanY);
                _shot.boundaryRecVzTop.col(it) = vz(horSpanX, topSpanY);
                _shot.boundaryRecTxxTop.col(it) = txx(horSpanX, topSpanY);
                _shot.boundaryRecTzzTop.col(it) = tzz(horSpanX, topSpanY);
                _shot.boundaryRecTxzTop.col(it) = txz(horSpanX, topSpanY);

                _shot.boundaryRecVxBottom.col(it) = vx(horSpanX, bottomSpanY);
                _shot.boundaryRecVzBottom.col(it) = vz(horSpanX, bottomSpanY);
                _shot.boundaryRecTxxBottom.col(it) = txx(horSpanX, bottomSpanY);
                _shot.boundaryRecTzzBottom.col(it) = tzz(horSpanX, bottomSpanY);
                _shot.boundaryRecTxzBottom.col(it) = txz(horSpanX, bottomSpanY);

                _shot.boundaryRecVxLeft.row(it) = vx(leftSpanX, verSpanY);
                _shot.boundaryRecVzLeft.row(it) = vz(leftSpanX, verSpanY);
                _shot.boundaryRecTxxLeft.row(it) = txx(leftSpanX, verSpanY);
                _shot.boundaryRecTzzLeft.row(it) = tzz(leftSpanX, verSpanY);
                _shot.boundaryRecTxzLeft.row(it) = txz(leftSpanX, verSpanY);

                _shot.boundaryRecVxRight.row(it) = vx(rightSpanX, verSpanY);
                _shot.boundaryRecVzRight.row(it) = vz(rightSpanX, verSpanY);
                _shot.boundaryRecTxxRight.row(it) = txx(rightSpanX, verSpanY);
                _shot.boundaryRecTzzRight.row(it) = tzz(rightSpanX, verSpanY);
                _shot.boundaryRecTxzRight.row(it) = txz(rightSpanX, verSpanY);
            } else {
                // Record wavefield on top and bottom boundary
#pragma omp parallel
#pragma omp for
                for (auto ix = (int) _currentModel.np_boundary;
                     ix < _currentModel.nx - _currentModel.np_boundary; ++ix) {
                    _shot.boundaryRecVxTop(ix - _currentModel.np_boundary, it) = vx(ix, 0);
                    _shot.boundaryRecVzTop(ix - _currentModel.np_boundary, it) = vz(ix, 0);
                    _shot.boundaryRecTxxTop(ix - _currentModel.np_boundary, it) = txx(ix, 0);
                    _shot.boundaryRecTzzTop(ix - _currentModel.np_boundary, it) = tzz(ix, 0);
                    _shot.boundaryRecTxzTop(ix - _currentModel.np_boundary, it) = txz(ix, 0);

                    _shot.boundaryRecVxBottom(ix - _currentModel.np_boundary, it) =
                            vx(ix, _currentModel.nz_domain - 1);
                    _shot.boundaryRecVzBottom(ix - _currentModel.np_boundary, it) =
                            vz(ix, _currentModel.nz_domain - 1);
                    _shot.boundaryRecTxxBottom(ix - _currentModel.np_boundary, it) =
                            txx(ix, _currentModel.nz_domain - 1);
                    _shot.boundaryRecTzzBottom(ix - _currentModel.np_boundary, it) =
                            tzz(ix, _currentModel.nz_domain - 1);
                    _shot.boundaryRecTxzBottom(ix - _currentModel.np_boundary, it) =
                            txz(ix, _currentModel.nz_domain - 1);
                }


#pragma omp parallel
#pragma omp for
                for (int iz = 1; iz < _currentModel.nz - _currentModel.np_boundary - 1; ++iz) {
                    _shot.boundaryRecVxLeft(it, iz) = vx(_currentModel.np_boundary, iz);
                    _shot.boundaryRecVzLeft(it, iz) = vz(_currentModel.np_boundary, iz);
                    _shot.boundaryRecTxxLeft(it, iz) = txx(_currentModel.np_boundary, iz);
                    _shot.boundaryRecTzzLeft(it, iz) = tzz(_currentModel.np_boundary, iz);
                    _shot.boundaryRecTxzLeft(it, iz) = txz(_currentModel.np_boundary, iz);

                    _shot.boundaryRecVxRight(it, iz) = vx(_currentModel.nx - _currentModel.np_boundary, iz);
                    _shot.boundaryRecVzRight(it, iz) = vz(_currentModel.nx - _currentModel.np_boundary, iz);
                    _shot.boundaryRecTxxRight(it, iz) = txx(_currentModel.nx - _currentModel.np_boundary, iz);
                    _shot.boundaryRecTzzRight(it, iz) = tzz(_currentModel.nx - _currentModel.np_boundary, iz);
                    _shot.boundaryRecTxzRight(it, iz) = txz(_currentModel.nx - _currentModel.np_boundary, iz);
                }
            }
        }

        if (it == _shot.nt - 1) {
            _shot.lastVx = vx;
            _shot.lastVz = vz;
            _shot.lastTxx = txx;
            _shot.lastTzz = tzz;
            _shot.lastTxz = txz;
            break;
        }
        // After this point is only integration, which doesn't have to be done at the last time level

        // Time integrate stress
#pragma omp parallel
#pragma omp for
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                if (iz > 1 and ix > 1 and ix < nx - 1 and iz < nz - 1) {
                    txx(ix, iz) = taper(ix, iz) *
                                  (txx(ix, iz) +
                                   (_shot.dt * _currentModel.lm(ix, iz) * (
                                           -coeff2 * vx(ix + 1, iz - 1) + coeff1 * vx(ix, iz - 1)
                                           - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 2, iz - 1)
                                   ) / dx +
                                    (_shot.dt * _currentModel.la(ix, iz)) * (
                                            -coeff2 * vz(ix - 1, iz + 1) + coeff1 * vz(ix - 1, iz)
                                            - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 1, iz - 2)
                                    ) / dz));
                    tzz(ix, iz) = taper(ix, iz) *
                                  (tzz(ix, iz) +
                                   (_shot.dt * _currentModel.la(ix, iz) * (
                                           -coeff2 * vx(ix + 1, iz - 1) + coeff1 * vx(ix, iz - 1)
                                           - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 2, iz - 1)
                                   ) / dx +
                                    (_shot.dt * _currentModel.lm(ix, iz)) * (
                                            -coeff2 * vz(ix - 1, iz + 1) + coeff1 * vz(ix - 1, iz)
                                            - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 1, iz - 2)
                                    ) / dz));
                    txz(ix, iz) = taper(ix, iz) *
                                  (txz(ix, iz) + _shot.dt * _currentModel.mu(ix, iz) * (
                                          (
                                                  -coeff2 * vx(ix - 1, iz + 1) + coeff1 * vx(ix - 1, iz)
                                                  - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 1, iz - 2)
                                          ) / dz +
                                          (
                                                  -coeff2 * vz(ix + 1, iz - 1) + coeff1 * vz(ix, iz - 1)
                                                  - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 2, iz - 1)
                                          ) / dx));
                } else {
                    txx(ix, iz) = txx(ix, iz) * taper(ix, iz);
                    txz(ix, iz) = txz(ix, iz) * taper(ix, iz);
                    tzz(ix, iz) = tzz(ix, iz) * taper(ix, iz);
                }
            }
        }

        // Time integrate velocity
#pragma omp parallel
#pragma omp for
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                if (iz < nz - 2 and ix < nx - 2 and ix > 0 and iz > 0) {
                    vx(ix, iz) =
                            taper(ix, iz) *
                            (vx(ix, iz) + _currentModel.b_vx(ix, iz) *
                                          (_shot.dt * (
                                                  -coeff2 * txx(ix + 2, iz + 1) + coeff1 * txx(ix + 1, iz + 1)
                                                  - coeff1 * txx(ix, iz + 1) + coeff2 * txx(ix - 1, iz + 1)
                                          ) / dx +
                                           _shot.dt * (
                                                   -coeff2 * txz(ix + 1, iz + 2) + coeff1 * txz(ix + 1, iz + 1)
                                                   - coeff1 * txz(ix + 1, iz) + coeff2 * txz(ix + 1, iz - 1)
                                           ) / dz));
                    vz(ix, iz) =
                            taper(ix, iz) *
                            (vz(ix, iz) + _currentModel.b_vz(ix, iz) *
                                          (_shot.dt * (
                                                  -coeff2 * txz(ix + 2, iz + 1) + coeff1 * txz(ix + 1, iz + 1)
                                                  - coeff1 * txz(ix, iz + 1) + coeff2 * txz(ix - 1, iz + 1)
                                          ) / dx +
                                           _shot.dt * (
                                                   -coeff2 * tzz(ix + 1, iz + 2) + coeff1 * tzz(ix + 1, iz + 1)
                                                   - coeff1 * tzz(ix + 1, iz) + coeff2 * tzz(ix + 1, iz - 1)
                                           ) / dz));
                } else {
                    vx(ix, iz) = vx(ix, iz) * taper(ix, iz);
                    vz(ix, iz) = vz(ix, iz) * taper(ix, iz);
                }
            }
        }
    }
}

void propagator::propagateBackward(model &_currentModel, shot &_shot) {
    // Some standard output
    std::cout << "    Stability number: " <<
              sqrt(_currentModel.lm.max() * _currentModel.b_vx.max()) * _shot.dt *
              sqrt(1.0 / (_currentModel.dx * _currentModel.dx) + 1.0 / (_currentModel.dz * _currentModel.dz))
              << std::endl;

    // Loading simulation parameters
    double dx = _currentModel.dx;
    const arma::uword nx = _currentModel.nx;
    double dz = _currentModel.dz;
    const arma::uword nz = _currentModel.nz;

    // Create dynamic fields
    arma::mat vx = -_shot.lastVx;
    arma::mat vz = -_shot.lastVz;
    arma::mat txx = _shot.lastTxx;
    arma::mat tzz = _shot.lastTzz;
    arma::mat txz = _shot.lastTxz;

    // Create taper matrix
    arma::mat taper = _currentModel.np_boundary * arma::ones(nx, nz);
    for (int iTaper = 0; iTaper < _currentModel.np_boundary; ++iTaper) {
        taper.submat(iTaper, 0, nx - iTaper - 1, nz - iTaper - 1) =
                1 + iTaper * arma::ones(nx - 2 * iTaper, nz - iTaper);
    }
    taper = arma::exp(-arma::square(_currentModel.np_factor * (_currentModel.np_boundary - taper)));

    // Create indices for injection
    auto horSpanX = arma::span(_currentModel.np_boundary, nx - _currentModel.np_boundary - 1);
    auto topSpanY = arma::span(1, 1);
    auto bottomSpanY = arma::span(_currentModel.nz_domain - 1, _currentModel.nz_domain - 1);
    auto verSpanY = arma::span(0 + 1, _currentModel.nz_domain - 2);
    auto leftSpanX = arma::span(_currentModel.np_boundary, _currentModel.np_boundary);
    auto rightSpanX = arma::span(nx - _currentModel.np_boundary - 1, nx - _currentModel.np_boundary - 1);

    arma::cube acc(nx, nz, static_cast<const arma::uword>(_shot.nt));

    // Time marching through all time levels
    for (int it = _shot.nt - 1; it >= 0; --it) {

        vx(horSpanX, topSpanY) = -_shot.boundaryRecVxTop.col(it);
        vz(horSpanX, topSpanY) = -_shot.boundaryRecVzTop.col(it);
        txx(horSpanX, topSpanY) = _shot.boundaryRecTxxTop.col(it);
        tzz(horSpanX, topSpanY) = _shot.boundaryRecTzzTop.col(it);
        txz(horSpanX, topSpanY) = _shot.boundaryRecTxzTop.col(it);

        vx(horSpanX, bottomSpanY) = -_shot.boundaryRecVxBottom.col(it);
        vz(horSpanX, bottomSpanY) = -_shot.boundaryRecVzBottom.col(it);
        txx(horSpanX, bottomSpanY) = _shot.boundaryRecTxxBottom.col(it);
        tzz(horSpanX, bottomSpanY) = _shot.boundaryRecTzzBottom.col(it);
        txz(horSpanX, bottomSpanY) = _shot.boundaryRecTxzBottom.col(it);

        vx(leftSpanX, verSpanY) = -_shot.boundaryRecVxLeft.row(it);
        vz(leftSpanX, verSpanY) = -_shot.boundaryRecVzLeft.row(it);
        txx(leftSpanX, verSpanY) = _shot.boundaryRecTxxLeft.row(it);
        tzz(leftSpanX, verSpanY) = _shot.boundaryRecTzzLeft.row(it);
        txz(leftSpanX, verSpanY) = _shot.boundaryRecTxzLeft.row(it);

        vx(rightSpanX, verSpanY) = -_shot.boundaryRecVxRight.row(it);
        vz(rightSpanX, verSpanY) = -_shot.boundaryRecVzRight.row(it);
        txx(rightSpanX, verSpanY) = _shot.boundaryRecTxxRight.row(it);
        tzz(rightSpanX, verSpanY) = _shot.boundaryRecTzzRight.row(it);
        txz(rightSpanX, verSpanY) = _shot.boundaryRecTxzRight.row(it);

        if(true){
            acc.slice(it) = vx; // takes a lot of ram
        }


        if (it == 0) {
            break;
        }
        // After this point is only integration, which doesn't have to be done at the last time level

        // Time integrate stress
#pragma omp parallel
#pragma omp for
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                if (iz > 1 and ix > 1 and ix < nx - 1 and iz < nz - 1) {
                    txx(ix, iz) = taper(ix, iz) *
                                  (txx(ix, iz) +
                                   (_shot.dt * _currentModel.lm(ix, iz) * (
                                           -coeff2 * vx(ix + 1, iz - 1) + coeff1 * vx(ix, iz - 1)
                                           - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 2, iz - 1)
                                   ) / dx +
                                    (_shot.dt * _currentModel.la(ix, iz)) * (
                                            -coeff2 * vz(ix - 1, iz + 1) + coeff1 * vz(ix - 1, iz)
                                            - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 1, iz - 2)
                                    ) / dz));
                    tzz(ix, iz) = taper(ix, iz) *
                                  (tzz(ix, iz) +
                                   (_shot.dt * _currentModel.la(ix, iz) * (
                                           -coeff2 * vx(ix + 1, iz - 1) + coeff1 * vx(ix, iz - 1)
                                           - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 2, iz - 1)
                                   ) / dx +
                                    (_shot.dt * _currentModel.lm(ix, iz)) * (
                                            -coeff2 * vz(ix - 1, iz + 1) + coeff1 * vz(ix - 1, iz)
                                            - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 1, iz - 2)
                                    ) / dz));
                    txz(ix, iz) = taper(ix, iz) *
                                  (txz(ix, iz) + _shot.dt * _currentModel.mu(ix, iz) * (
                                          (
                                                  -coeff2 * vx(ix - 1, iz + 1) + coeff1 * vx(ix - 1, iz)
                                                  - coeff1 * vx(ix - 1, iz - 1) + coeff2 * vx(ix - 1, iz - 2)
                                          ) / dz +
                                          (
                                                  -coeff2 * vz(ix + 1, iz - 1) + coeff1 * vz(ix, iz - 1)
                                                  - coeff1 * vz(ix - 1, iz - 1) + coeff2 * vz(ix - 2, iz - 1)
                                          ) / dx));
                } else {
                    txx(ix, iz) = txx(ix, iz) * taper(ix, iz);
                    txz(ix, iz) = txz(ix, iz) * taper(ix, iz);
                    tzz(ix, iz) = tzz(ix, iz) * taper(ix, iz);
                }
            }
        }

        // Time integrate velocity
#pragma omp parallel
#pragma omp for
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                if (iz < nz - 2 and ix < nx - 2 and ix > 0 and iz > 0) {
                    vx(ix, iz) =
                            taper(ix, iz) *
                            (vx(ix, iz) + _currentModel.b_vx(ix, iz) *
                                          (_shot.dt * (
                                                  -coeff2 * txx(ix + 2, iz + 1) + coeff1 * txx(ix + 1, iz + 1)
                                                  - coeff1 * txx(ix, iz + 1) + coeff2 * txx(ix - 1, iz + 1)
                                          ) / dx +
                                           _shot.dt * (
                                                   -coeff2 * txz(ix + 1, iz + 2) + coeff1 * txz(ix + 1, iz + 1)
                                                   - coeff1 * txz(ix + 1, iz) + coeff2 * txz(ix + 1, iz - 1)
                                           ) / dz));
                    vz(ix, iz) =
                            taper(ix, iz) *
                            (vz(ix, iz) + _currentModel.b_vz(ix, iz) *
                                          (_shot.dt * (
                                                  -coeff2 * txz(ix + 2, iz + 1) + coeff1 * txz(ix + 1, iz + 1)
                                                  - coeff1 * txz(ix, iz + 1) + coeff2 * txz(ix - 1, iz + 1)
                                          ) / dx +
                                           _shot.dt * (
                                                   -coeff2 * tzz(ix + 1, iz + 2) + coeff1 * tzz(ix + 1, iz + 1)
                                                   - coeff1 * tzz(ix + 1, iz) + coeff2 * tzz(ix + 1, iz - 1)
                                           ) / dz));
                } else {
                    vx(ix, iz) = vx(ix, iz) * taper(ix, iz);
                    vz(ix, iz) = vz(ix, iz) * taper(ix, iz);
                }
            }
        }
    }

#pragma omp parallel
#pragma omp for
    for (int it = 0; it < _shot.nt; ++it) { // takes a lot of time
        char filename[1024];
        sprintf(filename, "output/vx%i.txt", it);
        acc.slice(it).save(filename, arma::raw_ascii);
    }
}