#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include "MonteCarlo.hpp"
#include "BlackScholesModel.hpp"
#include "Options/Option.hpp"
#include <mpi.h>
#include "../src/parser.hpp"
#include "Options/OptionBasket.hpp"
#include "Options/AsianOption.hpp"
#include "Options/PerformanceOption.hpp"
#include <string>

using namespace std;

//////////////////////////////////// PROTOTYPE ///////////////////////////////////////////

/**
 * parserData is a structure due to store datas when parsing input file.
 */
typedef struct {
    string type;
    Option *option; /// Option
    double r; /// Free risk rate of the model
    double rho; /// Correlation between assets
    int nbSamples; /// Number of samples
    int nbHedging; /// Number of hedging dates
    double fdstep; /// fdstep number
    PnlVect* sigma; /// Asset volatilities
    PnlVect* spot; /// Asset spots
    PnlVect* trend; /// Asset trends
    PnlVect* divid; /// divid
} ParserDatas;

/**
 * Permet de pricer à t=0 une option pour différents paramètres
 *
 * @param[in] ParserDatas contenant l'option et les paramètres pour le pricing
 */
void priceAtZero(ParserDatas *datas);
/**
 * Permet de calculer le PnL d'une option avec différents paramètres
 *
 * @param[in] ParserDatas contenant l'option et les paramètres pour le pricing
 */
void computePnL(ParserDatas *datas);
/**
 * Permet de parser le fichier en input
 *
 * @param[in] infile est le chemin vers le fichier d'input
 * @return ParserDatas contenant les données nécessaires à la création du pricer
 */
ParserDatas *parseInputFile(char *infile);
/**
 * Permet de checker la bonne conformité des paramètres
 *
 * @param[in] les datas parsées au préalable
 * @thrown Des exceptions s'il y a des erreurs de conformité
 */
void checkParameters(ParserDatas *datas);

void priceAtZeroWithPrecision(ParserDatas *datas,double precision);

//////////////////////////////////// MAIN ///////////////////////////////////////////

int main(int argc, char **argv){
    MPI_Init(&argc,&argv);    
    ParserDatas *datas;
    char *infile;
    double test = 0;
    try {
        switch (argc) {
            case 2:
                infile = argv[1];
                datas = parseInputFile(infile);
                checkParameters(datas);
                priceAtZero(datas);
                break;
            case 3:
                infile = argv[1];
                datas = parseInputFile(infile);
                checkParameters(datas);

                test = atof(argv[2]);
                priceAtZeroWithPrecision(datas,test);
                break;
            default:
                throw invalid_argument(
                        "(Invalid Arguement) Please call pricer executable as follows : \n./pricer fichier_input \n./pricer -c fichier_input \n");
        }


    }catch(exception const& e){
        cerr << "[ERREUR] : " << e.what();
        return EXIT_FAILURE;
    }
    // Free
    free(datas);
    MPI_Finalize();
    return EXIT_SUCCESS;

}

//////////////////////////////////// IMPLEMENTATION ///////////////////////////////////////////

void displayParameters(ParserDatas *datas);
void priceAtZero(ParserDatas *datas) {
    double start = MPI_Wtime();
    // Model initialisation
    BlackScholesModel* model = new BlackScholesModel(
            datas->option->size_,datas->r,datas->rho,
            datas->sigma,datas->spot, datas->divid, datas->option->T_);

    // MonteCarlo initialisation
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    double mean = 0;
    double var = 0;
    double price = 0;
    double stdDev = 0;
    double ic = 0;
    MonteCarlo* monteCarlo = new MonteCarlo(model,datas->option,datas->nbSamples,rank);
    double espEstimation = 0;
    double varToAgregate = 0;
    if(rank != 0){

        // compute price
        //double price;

        int nbSamples_Slave;
        int temp = datas->nbSamples/(size - 1);
        nbSamples_Slave = (rank <= (datas->nbSamples % (size - 1))) ? (temp + 1) : temp;
        monteCarlo->price_slave(espEstimation,varToAgregate,nbSamples_Slave);

        // Display results


    }

    MPI_Reduce(&espEstimation,&mean,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&varToAgregate,&var,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);


    if(rank == 0){

        monteCarlo->price_master(price,stdDev,var,mean,datas->nbSamples);
        ic = 1.96*stdDev;
        displayParameters(datas);
        cout << "\n -----> Price [ " << price << " ]\n";
        cout << "\n -----> IC [ " << price - ic << " ; "<< price + ic << " ]" << endl;
        cout << "\n------> Standard Deviation : " << stdDev << endl;
        cout << "\n------> Number of Samples : " << monteCarlo->nbSamples_ << endl;
        cout << "\n------> Time of calculation : " << MPI_Wtime() - start << "seconds" << endl;


    }

    delete model;
    delete monteCarlo;

}


void priceAtZeroWithPrecision(ParserDatas *datas,double precision) {
    double start = MPI_Wtime();
    // Model initialisation
    BlackScholesModel* model = new BlackScholesModel(
            datas->option->size_,datas->r,datas->rho,
            datas->sigma,datas->spot, datas->divid,datas->option->T_);

    // MonteCarlo initialisation
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    double mean = 0;
    double var  = 0;
    double price = 0;
    double stdDev = 1000;
    double ic = 0;
    double nbSamples_Slave = 0;
    double nbSamples = 0;
    MonteCarlo* monteCarlo = new MonteCarlo(model,datas->option,datas->nbSamples,rank);
    double espEstimation = 0;
    double varToAgregate = 0;
    double nbSamples_tmp = 0;
    double var_tmp = 0;
    double mean_tmp = 0;

    while(stdDev > precision) {
        nbSamples_tmp = 0;
        mean_tmp = 0;
        var_tmp = 0;
        if (rank != 0) {
            nbSamples_Slave = 1;
            // compute price
            //double price;
            //monteCarlo->price(price,ic);
            monteCarlo->price_slave(espEstimation, varToAgregate, nbSamples_Slave);


        }
    //PB : MAJ de var et mean

        MPI_Reduce(&nbSamples_Slave, &nbSamples_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&espEstimation, &mean_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&varToAgregate, &var_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            nbSamples += nbSamples_tmp;
            var += var_tmp;
            mean += mean_tmp;
            price = 0;
            monteCarlo->price_master(price, stdDev, var, mean,nbSamples);

        }
        MPI_Bcast(&stdDev,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //Attetion on incrémente pas le nombre d'itération a faire

    }

    if (rank == 0) {
        ic = 1.96 * stdDev;

        displayParameters(datas);
        cout << "\n -----> Price [ " << price << " ]\n";
        cout << "\n -----> IC [ " << price - ic << " ; " << price + ic << " ]" << endl;
        cout << "\n------> Standard Deviation : " << stdDev << endl;
        cout << "\n------> Number of Samples : " << nbSamples << endl;
        cout << "\n------> Time of calculation : " << MPI_Wtime() - start << "seconds" << endl;


    }

    delete model;
    delete monteCarlo;

}

ParserDatas *parseInputFile(char *infile)
{
    double T, r, strike, rho, fdStep;
    PnlVect *spot, *sigma, *divid, *weights, *trend;
    string type;
    int size;
    int nbTimeStep;
    int H;
    size_t n_samples;

    Param *P = new Parser(infile);

    //Option params
    P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("strike", strike);
    P->extract("timestep number",nbTimeStep);
    P->extract("payoff coefficients",weights, size);
    //Other params
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("interest rate", r);
    P->extract("correlation",rho);
    P->extract("sample number", n_samples);
    P->extract("hedging dates number",H);
    P->extract("fd step", fdStep);

    if (!P->extract("dividend rate", divid, size))
        divid = pnl_vect_create_from_zero(size);

    // Previous compute of trend with previous data in back/ folder
    /*double trend_d = r + 0.01;
    PnlVect *trend = pnl_vect_create_from_scalar(size,trend_d);*/

    P->extract("trend", trend, size);

    //Filling the structure parserData
    ParserDatas *data = new ParserDatas();
    data -> type = type;
    data -> r = r;
    data -> rho = rho;
    data -> sigma = sigma;
    data -> spot = spot;
    data -> nbSamples = (int) n_samples;
    data -> trend = trend;
    if(H > 0) {
        data->nbHedging = H;
    }else{
        data->nbHedging = 3*nbTimeStep;
    }

    data-> fdstep = fdStep;
    //Creation de l'option basket
    if (type == "basket")
    {
        Option *basketOption = new OptionBasket(T,nbTimeStep,size,weights,strike);
        data->option = basketOption;
    }
        //Creation de l'option asian
    else if (type == "asian") {
        Option *asianOption = new AsianOption(T,nbTimeStep,size,weights,strike);
        data->option = asianOption;
    }
        //Creation de l'option performance
    else {
        Option *performanceOption = new PerformanceOption(T, nbTimeStep, size, weights);
        data->option = performanceOption;
    }

    return data;
}

void checkParameters(ParserDatas *datas){
    // Lever des erreurs de cette façon :
    // if (test[datas.x] is false)
    //      throw invalid_argument("message")

    //nbSamples
    if ( datas->nbSamples <= 0 )
        throw invalid_argument("nbSamples must be > 0");

    //spot
    if (pnl_vect_min(datas -> spot) < 0)
    {
        throw invalid_argument("Invalid negative spot in spot vector");
    }
    //frr
    if(datas -> r < 0)
    {
        throw invalid_argument("Invalid negative free risk rate ");
    }
    //volatility
    if (pnl_vect_min(datas -> sigma) < 0)
    {
        throw invalid_argument("Invalid negative volatility in volatility vector");
    }
    //nbHedging
    if (datas->nbHedging <= 0)
        throw invalid_argument("nbHedging must be > 0");
    if (datas->nbHedging <= datas->option->nbTimeSteps_)
        throw invalid_argument("nbHedging must be > nbTimeSteps");

}

void displayParameters(ParserDatas *datas){
    cout << "\n**** DATAS ****";
    cout << "\n* Free risk rate : " << datas->r;
    cout << "\n* Correlation : " << datas->rho;
    cout << "\n* Volatility (first) : " << GET(datas->sigma,0);
    cout << "\n* Samples number : " << datas->nbSamples;
    cout << "\n* Hedging dates number : " << datas->nbHedging;
    cout << "\n* Spots (first) : " << GET(datas->spot,0);
    //pnl_vect_print(datas->spot);
    if (datas->trend != NULL)
        cout << "\n* Trends (first) : " << GET(datas->trend,0);
    //pnl_vect_print(datas->trend);
    cout << "\n**** OPTION ****";
    cout << "\n* Type : " << datas->type;
    cout << "\n* Size : " << datas->option->size_;
    cout << "\n* Maturity : " << datas->option->T_;
    cout << "\n* Time step number : " << datas->option->nbTimeSteps_;
    cout << "\n*****************\n";
}
