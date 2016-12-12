#pragma once

#include "Options/Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "cstdlib"
#include <math.h>
#include <time.h>
#include "iostream"
#define DEFAULT_VALUE_NBSAMPLES 10000

class MonteCarlo
{
public:
    BlackScholesModel *mod_; /*! pointeur vers le modèle */
    Option *opt_; /*! pointeur sur l'option */
    PnlRng *rng_; /*! pointeur sur le générateur */
    int nbSamples_; /*! nombre de tirages Monte Carlo */
    PnlMat *path; /*! path for simulation : initialized one times only*/

    /**
     * Constructor
     *
     * @param[in] *model : Le modèle de Black Scholes à utiliser
     * @param[in] *option : L'option à pricer
     * @param[in] nbSamples : Le nombre de tirage de Monte Carlo
     */
    MonteCarlo(BlackScholesModel *model,Option *option, int nbSamples = DEFAULT_VALUE_NBSAMPLES, double rank=0);

    /**
     * Calcule le prix de l'option à la date 0
     *
     * @param[out] prix valeur de l'estimateur Monte Carlo
     * @param[out] ic largeur de l'intervalle de confiance
     */


    /**
     * Price master qui compute les résultats !
     *
     * @param[out] &prix: Calcul du prix à partir du discount factor et de la somme des payoffs simulés.
     * @param[out] &stdDev : Calcul de l'écart type
     * @param[in] nbSamples : Le nombre de tirage de Monte Carlo courant ! ( utile pour la partie avec précision )
     * @param[in] varEstimation : Somme des carrés des payoffs simulés
     * @param[in] espEstimation : Somme des prix simulé
     *
     */

    void price_master(double &prix, double &stdDev,double varEstimation,double espEstimation,double nbSamples);


    void price_slave(double &espEstimation, double &varEstimation, int nbSamples);

    /**
     *
     *Destructeur de MonteCarlo
     *
     */
    ~MonteCarlo();

private:
    /**
     * Permet une simulation du pay-off de l'option à partir du modèle
     *
     * @return Le pay-off de l'option avec une simulation du sous-jacent
     */

    double payOffSimulation();

};