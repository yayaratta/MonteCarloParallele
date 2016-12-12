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
    double strike; /// strike for packing
    int nbSamples; /// Number of samples
    int nbTimeStep; /// Time Step Number
    PnlVect* sigma; /// Asset volatilities
    PnlVect* spot; /// Asset spots
    PnlVect* weights; /// weights option for packing
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
void displayParameters(ParserDatas *datas);
//////////////////////////////////// MAIN ///////////////////////////////////////////

int main(int argc, char **argv){
    MPI_Init(&argc,&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ParserDatas  *datas;
    char *infile;
    double test = 0;
    try {
        switch (argc) {
            case 2:case 3:

                //infile = argv[1];
                //datas = parseInputFile(infile);
                //checkParameters(datas);

                if(rank == 0){

                    //PACKING
                    infile = argv[1];
                    datas = parseInputFile(infile);
                    checkParameters(datas);

                    double T, r, strike, rho, fdStep;
                    PnlVect *spot, *sigma, *weights;
                    int typeOption;
                    int size;
                    int nbSamples;


                    if(datas->type == "basket"){
                        typeOption = 0;
                    }else if(datas->type == "asian"){
                        typeOption = 1;
                    }else if(datas->type == "performance"){
                        typeOption = 2;
                    }else{
                        throw invalid_argument(
                                "Type Option invalid \n");
                    }



                    // packing des paramètres du problème : 1ère étape, récuperer la taille du buffer
                    char *buffer;
                    int info, count, buffersize=0, pos=0;
                    info = MPI_Pack_size(4, MPI_INT, MPI_COMM_WORLD, &count);
                    if (info) return info;
                    buffersize += count;
                    int pnlSize = (3 * datas->sigma->size);
                    info = MPI_Pack_size(4 + pnlSize , MPI_DOUBLE, MPI_COMM_WORLD, &count);
                    if (info) return info;
                    buffersize += count;


                    buffer = (char *)malloc(buffersize);

                    // Packing !!!!
                    cout << "ResPack : " << datas->nbTimeStep << std::endl;
                    cout << "Impression de weights : " << std::endl;
                    pnl_vect_print(datas->weights);
                    cout << endl;
                    info = MPI_Pack(&datas->sigma->size, 1, MPI_INT, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->nbSamples, 1, MPI_INT, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&typeOption, 1, MPI_INT, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->nbTimeStep, 1, MPI_INT, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->option->T_, 1, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->r, 1, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->strike, 1, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(&datas->rho, 1, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(datas->weights->array, datas->sigma->size, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(datas->spot->array, datas->sigma->size, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;
                    info = MPI_Pack(datas->sigma->array, datas->sigma->size, MPI_DOUBLE, buffer, buffersize, &pos, MPI_COMM_WORLD);
                    if(info) return info;


                    MPI_Bcast(&buffersize, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Bcast(buffer, buffersize, MPI_PACKED, 0, MPI_COMM_WORLD);

                }else{
                    datas = new ParserDatas();
                    char *buffer;
                    int info, buffersize,pos = 0,size=0,nbSamples,nbTimeStep;
                    double T = 0, r=0,strike = 0,rho;
                    int type = 0;
                    PnlVect * weights, *spot, *sigma;
                    info = MPI_Bcast(&buffersize, 1, MPI_INT, 0, MPI_COMM_WORLD);

                    if (info) return info;
                    buffer = (char *) malloc((unsigned) buffersize);
                    info = MPI_Bcast(buffer, buffersize, MPI_PACKED, 0, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &size, 1, MPI_INT, MPI_COMM_WORLD);
                    if (info) return info;

                    weights = pnl_vect_create_from_zero(size);
                    spot = pnl_vect_create_from_zero(size);
                    sigma = pnl_vect_create_from_zero(size);






                    info = MPI_Unpack(buffer, buffersize, &pos, &nbSamples, 1, MPI_INT, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &type, 1, MPI_INT, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &nbTimeStep, 1, MPI_INT, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &T, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &r, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &strike, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, &rho, 1, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, weights->array, size, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, spot->array, size, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;
                    info = MPI_Unpack(buffer, buffersize, &pos, sigma->array, size, MPI_DOUBLE, MPI_COMM_WORLD);
                    if (info) return info;



                    if(type == 0){
                        datas ->type = "basket";
                        Option *basketOption = new OptionBasket(T,nbTimeStep,size,weights,strike);
                        datas->option = basketOption;
                    }else if(type == 1){
                        datas->type = "asian";
                        Option *asianOption = new AsianOption(T,nbTimeStep,size,weights,strike);
                        datas->option = asianOption;

                    }else{
                        datas->type = "performance";
                        Option *performanceOption = new PerformanceOption(T, nbTimeStep, size, weights);
                        datas->option = performanceOption;
                    }

                    datas->r = r;
                    datas->nbSamples = nbSamples;
                    datas->nbTimeStep = nbTimeStep;
                    datas->strike = strike;
                    datas->rho = rho;
                    datas->weights = weights;
                    datas->spot = spot;
                    datas->sigma = sigma;



                    displayParameters(datas);
                    /*

                    cout << "SIZE : " << size << std::endl;

                    cout << "T : " << T << std::endl;

                    cout << "r : " << r << std::endl;

                    cout << "Nb Samples : " << datas->nbSamples << std::endl;

                    cout << "Type option: " << datas->type << std::endl;

                    cout << "Strike : " << datas->strike << std::endl;

                    cout << "nbTimeStep : " << datas->nbTimeStep << std::endl;

                    cout << "rho : " << datas->rho << std::endl;
                    cout << "Impression de weights : " << std::endl;
                    pnl_vect_print(datas->weights);
                    cout << endl;
                    cout << "Impression de sigma : " << std::endl;
                    pnl_vect_print(datas->sigma);
                    cout << endl;
                    cout << "Impression de spot : " << std::endl;
                    pnl_vect_print(datas->spot);
                    cout << endl;
*/
                    //UNPACKING

                }




                if(argc == 2)
                    priceAtZero(datas);
                else {
                    test = atof(argv[2]);
                    priceAtZeroWithPrecision(datas, test);
                }



                break;
            default:
                throw invalid_argument(
                        "(Invalid Arguement) Please call pricer executable as follows : \n./pricer fichier_input \n./pricer -c fichier_input \n");
        }


    }catch(exception const& e){
        cerr << "[ERREUR] : " << e.what();
        return EXIT_FAILURE;
    }


    // PENSER A FAIRE TOUS LES FREEEEEEE
    free(datas);

    // Free
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
            datas->sigma,datas->spot);

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
            datas->sigma,datas->spot);

    // MonteCarlo initialisation
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    double mean = 0;
    double var  = 0;
    double price = 0;
    double stdDev = 10000;
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
    data -> strike = strike;
    data -> nbTimeStep = nbTimeStep;
    data -> weights = weights;
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


}

void displayParameters(ParserDatas *datas){
    cout << "\n**** DATAS ****";
    cout << "\n* Free risk rate : " << datas->r;
    cout << "\n* Correlation : " << datas->rho;
    cout << "\n* Volatility (first) : " << GET(datas->sigma,0);
    cout << "\n* Samples number : " << datas->nbSamples;
    cout << "\n* Spots (first) : " << GET(datas->spot,0);

    cout << "\n**** OPTION ****";
    cout << "\n* Type : " << datas->type;
    cout << "\n* Size : " << datas->option->size_;
    cout << "\n* Maturity : " << datas->option->T_;
    cout << "\n* Time step number : " << datas->option->nbTimeSteps_;
    cout << "\n*****************\n";
}
