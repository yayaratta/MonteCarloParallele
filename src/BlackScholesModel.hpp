#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

class BlackScholesModel
{
private:
    PnlMat *L; /// matrice de cholesky du modèle
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales du sous-jacent
    PnlVect *divid_; /// valeurs des dividendes
    PnlVect *Gi_;
    PnlVect *LGi_;
    double T_;

    /**
     *
     * Constructeur de BloackScholesModel
     *
     */

    BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot, PnlVect *divid, double T = 3);

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


