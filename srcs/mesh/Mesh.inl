#include "Mesh.hpp"

inline double Mesh::getAlpha() const
{
    return m_alpha;
}

inline std::vector<std::size_t> Mesh::getElement(std::size_t elm) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    return m_elementsList[elm];
}

inline double Mesh::getElementDetJ(std::size_t elm) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_elementsDetJ.empty());
    return m_elementsDetJ[elm];
}

inline double Mesh::getElementInvJ(std::size_t elm, unsigned short i, unsigned short j) const
{
    assert(elm < m_elementsList.size() && "elm should be between 0 and size - 1 !");
    assert(!m_elementsInvJ.empty());
    assert(i < m_dim && j < m_dim && "i and j should be 0 or dim - 1");
    return m_elementsInvJ[elm][i][j];
}

inline std::size_t Mesh::getElementsNumber() const
{
    return m_elementsList.size();
}

//inline double Mesh::getFreeSurfaceDetJ(std::size_t edge) const
//{
//    assert(edge < m_freeSurfaceEdgesList.size() && "edge should be between 0 and size - 1 !");
//    return m_freeSurfaceEdgeDetJ[edge];
//}
//
//inline std::vector<std::size_t>  Mesh::getFreeSurfaceEdge(std::size_t edge) const
//{
//    return m_freeSurfaceEdgesList[edge];
//}
//
//inline std::size_t  Mesh::getFreeSurfaceEdgesNumber() const
//{
//    return m_freeSurfaceEdgesList.size();
//}

inline double Mesh::getGamma() const
{
    return m_gamma;
}

inline double Mesh::getGaussPoints(unsigned short point, unsigned short coordinate) const
{
    assert(coordinate < m_dim);

    constexpr std::array<std::array<double, 2>, 3> GL2Dpoints = {{{0.16666666666666666666666666666666, 0.16666666666666666666666666666666},
                                                                  {0.66666666666666666666666666666666, 0.16666666666666666666666666666666},
                                                                  {0.16666666666666666666666666666666, 0.66666666666666666666666666666666}}};

    constexpr std::array<std::array<double, 3>, 4> GL3Dpoints = {{{0.1381966011250105, 0.1381966011250105, 0.1381966011250105},
                                                                  {0.5854101966249685, 0.1381966011250105, 0.1381966011250105},
                                                                  {0.1381966011250105, 0.5854101966249685, 0.1381966011250105},
                                                                  {0.1381966011250105, 0.1381966011250105, 0.5854101966249685}}};
    switch(m_dim)
    {
        case 2:
            return GL2Dpoints[point][coordinate];
        case 3:
            return GL3Dpoints[point][coordinate];
        default:
            throw std::runtime_error("Unsupported dimension for gauss legendre quadrature!");
    }
}

inline unsigned short Mesh::getGaussPointsNumber() const
{
    switch(m_dim)
    {
        case 2:
            return 3;
        case 3:
            return 4;
        default:
            throw std::runtime_error("Unsupported dimension for gauss legendre quadrature!");
    }
}

inline double Mesh::getGaussWeight(unsigned short point) const
{
    constexpr std::array<double, 3> GL2Dweights = {0.33333333333333333333333333333333,
                                                   0.33333333333333333333333333333333,
                                                   0.33333333333333333333333333333333};

    constexpr std::array<double, 4> GL3Dweights = {0.25,
                                                   0.25,
                                                   0.25,
                                                   0.25};

    switch(m_dim)
    {
        case 2:
            return GL2Dweights[point];
        case 3:
            return GL3Dweights[point];
        default:
            throw std::runtime_error("Unsupported dimension for gauss legendre quadrature!");
    }
}

inline double Mesh::getHchar() const
{
    return m_hchar;
}

inline unsigned short Mesh::getMeshDim() const
{
    return m_dim;
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

inline double Mesh::getOmega() const
{
    return m_omega;
}

inline double Mesh::getRefElementSize() const
{
    switch(m_dim)
    {
        case 2:
            return 0.5;
        case 3:
            return 0.16666666666666666666666666666667;
        default:
            throw std::runtime_error("Unsupported dimension for gauss legendre quadrature!");
    }
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

void Mesh::setStatesNumber(unsigned short statesNumber)
{
    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        m_nodesList[n].states.resize(statesNumber);
    }
}
