#include "Mesh.hpp"

inline double Mesh::getDetJ(std::size_t elm) const
{
    assert(elm < elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!this->detJ.empty());
    return this->detJ[elm];
}

inline std::vector<std::size_t> Mesh::getElement(std::size_t elm) const
{
    assert(elm < elementList.size() && "elm should be between 0 and size - 1 !");
    return this->elementList[elm];
}

inline std::size_t Mesh::getElementNumber() const
{
    return this->elementList.size();
}

inline std::vector<bool> Mesh::getIndices() const
{
    std::vector<bool> indices(this->nodesList.size());

    for (std::size_t i = 0 ; i < this->nodesList.size() ; ++i)
    {
        if(this->nodesList[i].isBound || this->nodesList[i].isFree)
            indices[i] = true;
    }

    return indices;
}

inline double Mesh::getInvJ(std::size_t elm, unsigned short i, unsigned short j) const
{
    assert(elm < elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!this->invJ.empty());
    assert(i < 2 && j < 2 && "i and j should be 0 or 1 in 2D");
    return this->invJ[elm](i, j);
}

inline std::vector<Eigen::MatrixXd> Mesh::getN() const
{
    std::vector<Eigen::MatrixXd> Ns;

    for(auto point: GP2Dpoints<double>)
    {
        Eigen::MatrixXd N(2, 6); N.setZero();

        N(0,0) = N(1,3) =  1 - point[0] - point[1];
        N(0,1) = N(1,4) =  1 - point[0];
        N(0,2) = N(1,5) =  1 - point[1];

        Ns.push_back(N);
    }

    return Ns;
}

inline std::vector<Eigen::MatrixXd> Mesh::getB(std::size_t elm) const
{
    assert(elm < elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!this->invJ.empty());

    std::vector<Eigen::MatrixXd> Bs;

    for(auto point: GP2Dpoints<double>)
    {
        Eigen::MatrixXd B(3, 6); B.setZero();

        B(0,0) = B(2,0) =  - this->invJ[elm](0,0) - this->invJ[elm](0,1);
        B(0,1) = B(2,1) =  - this->invJ[elm](0,0);
        B(0,2) = B(2,2) =  - this->invJ[elm](0,1);

        B(1,3) = B(2,3) =  - this->invJ[elm](1,0) - this->invJ[elm](1,1);
        B(1,4) = B(2,4) =  - this->invJ[elm](1,0);
        B(1,5) = B(2,5) =  - this->invJ[elm](1,1);

        Bs.push_back(B) ;
    }

    return Bs;
}
