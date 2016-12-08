#include "cstdlib"
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "iostream"
using namespace std;

#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(BlackScholesModel *model, Option *option, int nbSamples, double fdStep,double rank) {
    mod_ = model;
    opt_ = option;
    rng_ = pnl_rng_dcmt_create_id(rank,PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng_,time(NULL));
    //rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
    //pnl_rng_sseed(rng_, time(NULL));
    fdStep_ = (fdStep == 0) ? DEFAULT_VALUE_FDSTEP : fdStep;
    nbSamples_ = nbSamples;
    path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
    pathShifted = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
}



// Price at t=0
void MonteCarlo::price(double &prix, double &ic) {

    double T = opt_->T_;
    double r = mod_->r_;

    double discountFactor = exp(-r*T);

    double espEstimation = 0;
    double varEstimateur = 0;
    for (int j = 0; j < nbSamples_; ++j) {
        double estimation = payOffSimulation();
        espEstimation += estimation;
        varEstimateur += estimation * estimation;
    }

    espEstimation /= (double)nbSamples_;
    varEstimateur /= (double)nbSamples_;
    varEstimateur = exp(-2*r*T)*fabs(varEstimateur - espEstimation * espEstimation);

    prix = discountFactor * espEstimation;
    ic = (1.96) * sqrt(varEstimateur/(double)nbSamples_);

}

void MonteCarlo::price_master(double &prix, double &stdDev, double &varEstimateur,double &espEstimation) {

    double T = opt_->T_;
    double r = mod_->r_;
    double discountFactor = exp(-r*T);
    prix = discountFactor * espEstimation;

    espEstimation /= (double)nbSamples_;
    varEstimateur /= (double)nbSamples_;
    varEstimateur = exp(-2*r*T)*fabs(varEstimateur - espEstimation * espEstimation);

    prix = discountFactor * espEstimation;
    stdDev = sqrt(varEstimateur/(double)nbSamples_);


}

void MonteCarlo::price_slave(double &espEstimation, double &varEstimateur,int nbSamples) {

    double T = opt_->T_;
    double r = mod_->r_;


    espEstimation = 0;
    varEstimateur = 0;
    for (int j = 0; j < nbSamples; ++j) {
        double estimation = payOffSimulation();
        espEstimation += estimation;
        varEstimateur += estimation * estimation;
    }

}




double MonteCarlo::payOffSimulation() {
    double payOffSimu = 0;
    double T = opt_->T_;
    int nbTimeSteps = opt_->nbTimeSteps_;

    // Simulation
    mod_->asset(path, T, nbTimeSteps, rng_);

    payOffSimu = opt_->payoff(path);

    return payOffSimu;
}



MonteCarlo::~MonteCarlo() {
    // delete mod_;
    // delete opt_;
    pnl_rng_free(&rng_);
    pnl_mat_free(&path);
    pnl_mat_free(&pathShifted);
}
