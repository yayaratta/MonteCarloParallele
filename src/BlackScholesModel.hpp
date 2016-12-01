#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
private:
    PnlMat *L; /// matrice de cholesky du modèle
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *trend_; /// vecteur de tendances instantanées en probabilité historique
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales du sous-jacent
    PnlVect *Gi_;
    PnlVect *LGi_;
    int H_;
    double T_;

    BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot, PnlVect *trend = NULL, int H = 100,double T = 3);

/**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (N+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);



    /**
     *
     *Destructeur de BlackScholesModel
     *
     */
    ~BlackScholesModel();

};


