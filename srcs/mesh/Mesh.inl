#include "Mesh.hpp"

#include <array>
#include <cassert>
#include <iostream>

inline double Mesh::getAlpha() const noexcept
{
    return m_alpha;
}

inline unsigned short Mesh::getDim() const noexcept
{
    return m_dim;
}

inline std::vector<IndexType> Mesh::getElement(IndexType elm) const noexcept
{
    assert(elm < m_elementsList.size());

    return m_elementsList[elm].nodesIndexes;
}

inline double Mesh::getElementDetJ(IndexType elm) const noexcept
{
    assert(elm < m_elementsList.size());

    return m_elementsList[elm].detJ;
}

inline double Mesh::getElementInvJ(IndexType elm, unsigned short i, unsigned short j) const noexcept
{
    assert(elm < m_elementsList.size());
    assert(i < m_dim + 1 && j < m_dim + 1);

    return m_elementsList[elm].invJ[i][j];
}

inline IndexType Mesh::getElementsNumber() const noexcept
{
    return m_elementsList.size();
}

//inline double Mesh::getFreeSurfaceDetJ(IndexType edge) const
//{
//    assert(edge < m_freeSurfaceEdgesList.size() && "edge should be between 0 and size - 1 !");
//    return m_freeSurfaceEdgeDetJ[edge];
//}
//
//inline std::vector<IndexType>  Mesh::getFreeSurfaceEdge(IndexType edge) const
//{
//    return m_freeSurfaceEdgesList[edge];
//}
//
//inline IndexType  Mesh::getFreeSurfaceEdgesNumber() const
//{
//    return m_freeSurfaceEdgesList.size();
//}

inline double Mesh::getGamma() const noexcept
{
    return m_gamma;
}

inline double Mesh::getGaussPoints(unsigned short point, unsigned short coordinate) const noexcept
{
    assert(coordinate < m_dim);

    constexpr std::array<std::array<double, 2>, 3> GL2Dpoints {{{0.16666666666666666666666666666666, 0.16666666666666666666666666666666},
                                                                {0.66666666666666666666666666666666, 0.16666666666666666666666666666666},
                                                                {0.16666666666666666666666666666666, 0.66666666666666666666666666666666}}};

    constexpr std::array<std::array<double, 3>, 4> GL3Dpoints {{{0.1381966011250105, 0.1381966011250105, 0.1381966011250105},
                                                                {0.5854101966249685, 0.1381966011250105, 0.1381966011250105},
                                                                {0.1381966011250105, 0.5854101966249685, 0.1381966011250105},
                                                                {0.1381966011250105, 0.1381966011250105, 0.5854101966249685}}};

    assert(point < m_dim + 1);

    switch(m_dim)
    {
        case 2:
            return GL2Dpoints[point][coordinate];
        default:
            return GL3Dpoints[point][coordinate];
    }
}

inline unsigned short Mesh::getGaussPointsNumber() const noexcept
{
    switch(m_dim)
    {
        case 2:
            return 3;
        default:
            return 4;
    }
}

inline double Mesh::getGaussWeight(unsigned short point) const noexcept
{
    constexpr std::array<double, 3> GL2Dweights {0.33333333333333333333333333333333,
                                                 0.33333333333333333333333333333333,
                                                 0.33333333333333333333333333333333};

    constexpr std::array<double, 4> GL3Dweights {0.25,
                                                 0.25,
                                                 0.25,
                                                 0.25};

    switch(m_dim)
    {
        case 2:
            return GL2Dweights[point];
        default:
            return GL3Dweights[point];
    }
}

inline double Mesh::getHchar() const noexcept
{
    return m_hchar;
}

inline IndexType Mesh::getNodesNumber() const noexcept
{
    return m_nodesList.size();
}

inline std::vector<double> Mesh::getNodeInitialPosition(IndexType nodeIndex) const
{
    assert(nodeIndex < m_nodesList.size());

    if(!m_nodesList[nodeIndex].isBound)
        throw std::runtime_error("initial position is only available for boundary nodes!");

    return m_nodesList[nodeIndex].initialPosition;

}

inline std::vector<double> Mesh::getNodePosition(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());

    return m_nodesList[nodeIndex].position;
}

inline double Mesh::getNodePosition(IndexType nodeIndex, unsigned short coordinate) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    assert(coordinate < m_dim);

    return m_nodesList[nodeIndex].position[coordinate];
}

inline double Mesh::getNodeState(IndexType nodeIndex, unsigned short state) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    assert(state < m_nodesList[nodeIndex].states.size());

    return m_nodesList[nodeIndex].states[state];
}

inline std::string Mesh::getNodeType(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());

    return m_tagNames[m_nodesList[nodeIndex].tag];
}

inline double Mesh::getOmega() const noexcept
{
    return m_omega;
}

inline double Mesh::getRefElementSize() const noexcept
{
    switch(m_dim)
    {
        case 2:
            return 0.5;
        default:
            return 0.16666666666666666666666666666667;
    }
}

inline bool Mesh::isNodeFree(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    return m_nodesList[nodeIndex].isFree;
}

inline bool Mesh::isNodeBound(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    return m_nodesList[nodeIndex].isBound;
}

inline bool Mesh::isNodeDirichlet(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    return m_nodesList[nodeIndex].isDirichlet;
}

inline bool Mesh::isNodeOnFreeSurface(IndexType nodeIndex) const noexcept
{
    assert(nodeIndex < m_nodesList.size());
    return m_nodesList[nodeIndex].isOnFreeSurface;
}

inline void Mesh::setNodeIsDirichlet(IndexType nodeIndex, bool isDirichlet) noexcept
{
    assert(nodeIndex < m_nodesList.size());
    m_nodesList[nodeIndex].isDirichlet = isDirichlet;
}

inline void Mesh::setNodeState(IndexType nodeIndex, unsigned short stateIndex, double state) noexcept
{
    assert(nodeIndex < m_nodesList.size());
    m_nodesList[nodeIndex].states[stateIndex] = state;
}

void Mesh::setStatesNumber(unsigned short statesNumber)
{
    #pragma omp parallel for default(shared)
    for(IndexType n = 0 ; n < m_nodesList.size() ; ++n)
    {
        m_nodesList[n].states.resize(statesNumber);
    }
}
