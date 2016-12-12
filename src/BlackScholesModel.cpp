#include <iostream>
#include <stdexcept>
#include "cstdlib"


using namespace std;


#include "BlackScholesModel.hpp"

#define DEFAULT_VALUE_FOR_TREND 0.03

//////////////// PROTOTYPES //////////////////
/**
 * Permit to get Cholesky matrix from correlation matrix Cor = (rho)*ones + (1-rho)Id
 */
PnlMat *getCholeskyFromRho(int size, double rho);


//////////////// IMPLEMENTATION OF BLACKSCHOLESMODEL CLASS ///////////////////////////

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot, PnlVect *trend, int H,double T) :
        size_(size), r_(r), rho_(rho), sigma_(sigma), spot_(spot), H_(H),T_(T){
    if (sigma->size != size)
        throw new std::invalid_argument("Size of sigma is different of the number of shares");
    if (spot->size != size)
        throw new std::invalid_argument("Size of spot is different of the number of shares");
    if (rho >= 1 || rho <= (-(double)1/(size-1)))
        throw new std::invalid_argument("Correlation not in ]-1/(D-1);1[");
    if (trend == NULL) {
        cout << "## WARNING : Trend is null, default value " << DEFAULT_VALUE_FOR_TREND << " is used";
        trend_ = pnl_vect_create_from_scalar(size,DEFAULT_VALUE_FOR_TREND);
    }
    else
        trend_ = trend;
    if (trend_->size != size)
        throw new std::invalid_argument("Size of trend is different of the number of shares");

    // Cholesky initialisation (computed just one time and not for each asset call)
    L = getCholeskyFromRho(size,rho);
    Gi_ = pnl_vect_new();
    LGi_ = pnl_vect_new();
}



void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng) {
    // Case multidimensional with correlation
    // Step initialisation
    double step = T/(double)nbTimeSteps;
    double sqrtStep = sqrt(step);
    double sigma_d;
    double Sd_tiMinus1;
    double LdGi;
    double Sd_ti;
    // Spot initialisation
    for (int d = 0; d < path->n; ++d)
        PNL_MSET(path, 0, d, GET(spot_,d));

    for (int i = 1; i < path->m; ++i) {
        pnl_vect_rng_normal_d(Gi_,path->n,rng); // Gi gaussian vector
        LGi_ =  pnl_mat_mult_vect(L,Gi_); // All the LdGi
        // For each asset
        for (int d = 0; d < path->n; ++d) {
            sigma_d = GET(sigma_,d);
            Sd_tiMinus1 = PNL_MGET(path, (i-1), d);
            LdGi = GET(LGi_,d);
            Sd_ti = Sd_tiMinus1 * exp((r_ - sigma_d * sigma_d / 2) * step + sigma_d * sqrtStep * LdGi);
            PNL_MSET(path,i,d,Sd_ti);
        }
    }
}


BlackScholesModel::~BlackScholesModel() {
    pnl_mat_free(&L);
    pnl_vect_free(&trend_);
    pnl_vect_free(&spot_);
    pnl_vect_free(&sigma_);
    pnl_vect_free(&Gi_);
    pnl_vect_free(&LGi_);

}

/////////////////////// IMPLEMENTATION OF FUNCTIONS /////////////////////////

PnlMat *getCholeskyFromRho(int size, double rho) {
    // Correlation matrix
    PnlMat *res = pnl_mat_create_from_scalar(size,size,rho);
    for (int i = 0; i < res->n; ++i)
        PNL_MSET(res,i,i,1);
    // Cholesky of the correlation matrix
    pnl_mat_chol(res);
    return res;
}
