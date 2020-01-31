#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include "Solver.hpp"
#include "../quadrature/gausslegendre.hpp"

Solver::Solver(const Params& params, Mesh& mesh, ProblemType problemType) :
params(params), mesh(mesh), problemType(problemType), currentTimeStep(params.maxTimeStep)
{
    if(problemType == INCOMPRESSIBLE_PSPG)
    {
        matrices.resize(5);
        vectors.resize(2);
        sideVariables.resize(1);
    }
    else
        throw std::runtime_error("Unsupported problem type.");
}

Solver::~Solver()
{

}

void Solver::solveProblem()
{
    auto startTime = std::chrono::high_resolution_clock::now();
    this->mesh.remesh();
    auto endTime = std::chrono::high_resolution_clock::now();

    auto ellapsedTime =
        std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "Ellapsed time for remeshing: "
              << static_cast<double>(ellapsedTime.count())/1000.0
              << " s" << std::endl;
    this->mesh.computeDetInvJ();

    _computeSideVariables();
    _buildMatrices({true, true, true, true, true}, {true, true});
}

void Solver::_buildMatrices(std::vector<bool> matricesToBuild, std::vector<bool> vectorsToBuild)
{
    assert(matricesToBuild.size() == this->matrices.size());
    assert(vectorsToBuild.size() == this->vectors.size());

    std::vector<Eigen::MatrixXd> N = this->mesh.getN();
    std::vector<std::vector<Eigen::MatrixXd>> B;

    for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
    {
        B.push_back(this->mesh.getB(elm));
    }

    assert(N.size() == 3);

    const Eigen::MatrixXd NrightP = (Eigen::MatrixXd(6,3) << 1, 0, 0,
                                                             0, 1, 0,
                                                             0, 0, 1,
                                                             0, 0, 0,
                                                             0, 0, 0,
                                                             0, 0, 0).finished();

    const Eigen::MatrixXd NleftP = (Eigen::MatrixXd(1,2) << 1, 0).finished();

    const Eigen::MatrixXd BrightP = (Eigen::MatrixXd(6,3) << 1, 0, 0,
                                                             0, 1, 0,
                                                             0, 0, 1,
                                                             1, 0, 0,
                                                             0, 1, 0,
                                                             0, 0, 1).finished();

    const Eigen::MatrixXd BleftP = (Eigen::MatrixXd(2,3) << 1, 0, 0,
                                                            0, 1, 0).finished();

    /********************************************************************************
                                        Build M
    ********************************************************************************/
    if(matricesToBuild[0] == true)
    {
        this->matrices[0].resize(2*this->mesh.nodesList.size(), 2*this->mesh.nodesList.size());

        Eigen::MatrixXd Me(6, 6); Me.setZero();
        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Me = Me + N[k].transpose()*N[k]*GP2Dweight<double>[k];

        Me = this->params.fluidParameters[0]*Me;

        std::vector<Eigen::Triplet<double>> index;
        index.reserve(2*3*3*this->mesh.getElementNumber());

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j],
                                                           Me(i,j)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i] + this->mesh.nodesList.size(),
                                                           this->mesh.getElement(elm)[j] + this->mesh.nodesList.size(),
                                                           Me(i,j)*this->mesh.getDetJ(elm)));
                }
            }
        }

        this->matrices[0].setFromTriplets(index.begin(), index.end(),
                                          [](const double& a,const double& b){
                                              return a + b;
                                          });
    }

    /********************************************************************************
                                        Build K
    ********************************************************************************/
    if(matricesToBuild[1] == true)
    {
        this->matrices[1].resize(2*this->mesh.nodesList.size(), 2*this->mesh.nodesList.size());

        const Eigen::DiagonalMatrix<double, 3> ddev(2, 2, 1);

        std::vector<Eigen::Triplet<double>> index;
        index.reserve(3*3*this->mesh.getElementNumber());

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            Eigen::MatrixXd Ke = this->params.fluidParameters[1]*B[elm][0].transpose()*ddev*B[elm][0]; //same matrices for all gauss point ^^

            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j],
                                                           Ke(i,j)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j] + this->mesh.nodesList.size(),
                                                           Ke(i,j+3)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i] + this->mesh.nodesList.size(),
                                                           this->mesh.getElement(elm)[j],
                                                           Ke(i+3,j)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i] + this->mesh.nodesList.size(),
                                                           this->mesh.getElement(elm)[j] + this->mesh.nodesList.size(),
                                                           Ke(i+3,j+3)*this->mesh.getDetJ(elm)));
                }
            }
        }

        this->matrices[1].setFromTriplets(index.begin(), index.end(),
                                          [](const double& a,const double& b){
                                              return a + b;
                                          });
    }

    /********************************************************************************
                                        Build D
    ********************************************************************************/
    if(matricesToBuild[2] == true)
    {
        this->matrices[2].resize(this->mesh.nodesList.size(), 2*this->mesh.nodesList.size());

        const Eigen::Vector3d m(1, 1, 0);

        std::vector<Eigen::Triplet<double>> index;
        index.reserve(3*6*this->mesh.getElementNumber());

        Eigen::MatrixXd Dprev(3,1); Dprev.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Dprev = Dprev + (NleftP*N[k]*NrightP).transpose()*GP2Dweight<double>[k];

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            Eigen::MatrixXd De = Dprev*m.transpose()*B[elm][0]; //same matrices for all gauss point ^^

            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j],
                                                           De(i,j)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j] + this->mesh.nodesList.size(),
                                                           De(i,j+3)*this->mesh.getDetJ(elm)));
                }
            }
        }

        this->matrices[2].setFromTriplets(index.begin(), index.end(),
                                          [](const double& a,const double& b){
                                              return a + b;
                                          });
    }

    /********************************************************************************
                                        Build C
    ********************************************************************************/
    if(matricesToBuild[3] == true)
    {
        this->matrices[3].resize(this->mesh.nodesList.size(), 2*this->mesh.nodesList.size());

        std::vector<Eigen::Triplet<double>> index;
        index.reserve(3*6*this->mesh.getElementNumber());

        Eigen::MatrixXd Cprev(2,6); Cprev.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Cprev = Cprev + N[k]*GP2Dweight<double>[k];

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            Eigen::MatrixXd Ce = this->sideVariables[0][elm]*(BleftP*B[elm][0]*BrightP).transpose()*Cprev;

            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j],
                                                           Ce(i,j)*this->mesh.getDetJ(elm)));

                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j] + this->mesh.nodesList.size(),
                                                           Ce(i,j+3)*this->mesh.getDetJ(elm)));
                }
            }
        }

        this->matrices[3].setFromTriplets(index.begin(), index.end(),
                                          [](const double& a,const double& b){
                                              return a + b;
                                          });
    }

    /********************************************************************************
                                        Build L
    ********************************************************************************/
    if(matricesToBuild[4] == true)
    {
        this->matrices[4].resize(this->mesh.nodesList.size(), this->mesh.nodesList.size());

        std::vector<Eigen::Triplet<double>> index;
        index.reserve(3*3*this->mesh.getElementNumber());

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            Eigen::MatrixXd Le = (this->sideVariables[0][elm]/this->params.fluidParameters[0])
                                 *(BleftP*B[elm][0]*BrightP).transpose()*(BleftP*B[elm][0]*BrightP);

            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                for(unsigned short j = 0 ; j < 3 ; ++j)
                {
                    index.push_back(Eigen::Triplet<double>(this->mesh.getElement(elm)[i],
                                                           this->mesh.getElement(elm)[j],
                                                           Le(i,j)*this->mesh.getDetJ(elm)));
                }
            }
        }

        this->matrices[4].setFromTriplets(index.begin(), index.end(),
                                          [](const double& a,const double& b){
                                              return a + b;
                                          });
    }

    /********************************************************************************
                                        Build f
    ********************************************************************************/
    if(vectorsToBuild[0] == true)
    {
        this->vectors[0].setZero();
        this->vectors[0].resize(2*this->mesh.nodesList.size());

        const Eigen::Vector2d b(0, -this->params.gravity);

        Eigen::MatrixXd Fprev(6,2); Fprev.setZero();

        for (unsigned short k = 0 ; k < N.size() ; ++k)
            Fprev = Fprev + N[k].transpose()*GP2Dweight<double>[k];

        Eigen::VectorXd Fe = this->params.fluidParameters[0]*Fprev*b;

        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                this->vectors[0](this->mesh.getElement(elm)[i]) =
                this->vectors[0](this->mesh.getElement(elm)[i]) + Fe(i)*this->mesh.getDetJ(elm);

                this->vectors[0](this->mesh.getElement(elm)[i] + this->mesh.nodesList.size()) =
                this->vectors[0](this->mesh.getElement(elm)[i] + this->mesh.nodesList.size()) + Fe(i + 3)*this->mesh.getDetJ(elm);
            }
        }
    }

    /********************************************************************************
                                        Build h
    ********************************************************************************/
    if(vectorsToBuild[1] == true)
    {
        this->vectors[1].setZero();
        this->vectors[1].resize(this->mesh.nodesList.size());

        const Eigen::Vector2d b(0, -this->params.gravity);


        for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
        {
            Eigen::VectorXd He = this->sideVariables[0][elm]*(BleftP*B[elm][0]*BrightP).transpose()*b;

            for(unsigned short i = 0 ; i < 3 ; ++i)
            {
                this->vectors[1](this->mesh.getElement(elm)[i]) =
                this->vectors[1](this->mesh.getElement(elm)[i]) + He(i)*this->mesh.getDetJ(elm);
            }
        }
    }
}

void Solver::_computeSideVariables()
{
    sideVariables[0].resize(this->mesh.getElementNumber());

    double rho = this->params.fluidParameters[0];
    double mu = this->params.fluidParameters[1];

    double U = 0;
    std::size_t trueNnodes = 0;
    for(std::size_t n = 0 ; n < this->mesh.nodesList.size() ; ++n)
    {
        if(this->mesh.nodesList[n].isBound == false ||
           (this->mesh.nodesList[n].isBound == true && this->mesh.nodesList[n].isFree == false))
        {
            U += std::sqrt(this->mesh.nodesList[n].states[0]*
               this->mesh.nodesList[n].states[0] +
               this->mesh.nodesList[n].states[1]*
               this->mesh.nodesList[n].states[1]);

            trueNnodes++;
        }
    }
    U /= static_cast<double>(trueNnodes);

    for(std::size_t elm = 0 ; elm < this->mesh.getElementNumber() ; ++elm)
    {
        double h = std::sqrt(2*this->mesh.getDetJ(elm)/M_PI);

        sideVariables[0][elm] = std::sqrt((2/this->currentTimeStep)*(2/this->currentTimeStep)
                                 + (2*U/h)*(2*U/h)
                                 + 9*(4*mu/(h*h*rho))*(4*mu/(h*h*rho)));
    }
}
