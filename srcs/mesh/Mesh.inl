#include "Mesh.hpp"

inline std::vector<Eigen::MatrixXd> Mesh::getB(std::size_t elm) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_invJ.empty());

    std::vector<Eigen::MatrixXd> Bs;

    for(auto point: GP2Dpoints<double>)
    {
        Eigen::MatrixXd B(3, 6); B.setZero();

        B(0,0) = B(2,0) =  - m_invJ[elm][0][0] - m_invJ[elm][1][0];
        B(0,1) = B(2,1) =  m_invJ[elm][0][0];
        B(0,2) = B(2,2) =  m_invJ[elm][1][0];

        B(1,3) = B(2,3) =  - m_invJ[elm][0][1] - m_invJ[elm][1][1];
        B(1,4) = B(2,4) =  m_invJ[elm][0][1];
        B(1,5) = B(2,5) =  m_invJ[elm][1][1];

        Bs.push_back(B) ;
    }

    return Bs;
}

inline double Mesh::getDetJ(std::size_t elm) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_detJ.empty());
    return m_detJ[elm];
}

inline std::vector<std::size_t> Mesh::getElement(std::size_t elm) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    return m_elementList[elm];
}

inline std::size_t Mesh::getElementNumber() const
{
    return m_elementList.size();
}

inline double Mesh::getInvJ(std::size_t elm, unsigned short i, unsigned short j) const
{
    assert(elm < m_elementList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_invJ.empty());
    assert(i < 2 && j < 2 && "i and j should be 0 or 1 in 2D");
    return m_invJ[elm][i][j];
}

inline unsigned short Mesh::getMeshDim() const
{
    return m_dim;
}

inline std::vector<Eigen::MatrixXd> Mesh::getN() const
{
    std::vector<Eigen::MatrixXd> Ns;

    for(auto point: GP2Dpoints<double>)
    {
        Eigen::MatrixXd N(2, 6); N.setZero();

        N(0,0) = N(1,3) =  1 - point[0] - point[1];
        N(0,1) = N(1,4) =  point[0];
        N(0,2) = N(1,5) =  point[1];

        Ns.push_back(N);
    }

    return Ns;
}

inline std::size_t Mesh::getNodesNumber() const
{
    return m_nodesList.size();
}

inline double Mesh::getNodePosition(std::size_t nodeIndex, unsigned short coordinate) const
{
    return m_nodesList[nodeIndex].position[coordinate];
}

inline double Mesh::getNodeState(std::size_t nodeIndex, unsigned short state) const
{
    return m_nodesList[nodeIndex].states[state];
}

inline bool Mesh::isNodeFree(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isFree;
}

inline bool Mesh::isNodeBound(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isBound;
}

inline bool Mesh::isNodeFluidInput(std::size_t nodeIndex) const
{
    return m_nodesList[nodeIndex].isFluidInput;
}

inline void Mesh::setNodeState(std::size_t nodeIndex, unsigned short stateIndex, double state)
{
    m_nodesList[nodeIndex].states[stateIndex] = state;
}
