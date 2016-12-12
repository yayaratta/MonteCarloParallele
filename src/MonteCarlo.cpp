
//using namespace std;

#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(BlackScholesModel *model, Option *option, int nbSamples, double rank) {
    mod_ = model;
    opt_ = option;
    rng_ = pnl_rng_dcmt_create_id(rank,PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng_,time(NULL));
    nbSamples_ = nbSamples;
    path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
}


void MonteCarlo::price_master(double &prix, double &stdDev, double varEstimateur,double espEstimation,double nbSamples) {
    double T = opt_->T_;
    double r = mod_->r_;
    double discountFactor = exp(-r*T);

    espEstimation /= nbSamples;
    varEstimateur /= nbSamples;
    varEstimateur = exp(-2*r*T)*fabs(varEstimateur - espEstimation * espEstimation);
    prix = discountFactor * espEstimation;
    stdDev = sqrt(varEstimateur/nbSamples);
}

void MonteCarlo::price_slave(double &espEstimation, double &varEstimateur,int nbSamples) {

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
    pnl_rng_free(&rng_);
    pnl_mat_free(&path);
}
