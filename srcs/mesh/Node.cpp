#include "Node.hpp"

#include <cmath>
#include "Mesh.hpp"
#include "Element.hpp"
#include "Facet.hpp"

Node::Node(Mesh& mesh):
m_pMesh(&mesh)
{
    m_neighbourNodes.clear();
    m_elements.clear();
    m_facets.clear();

    m_isOnBoundary = false;
}

const Element& Node::getElement(unsigned int elementIndex) const noexcept
{
    return m_pMesh->getElement(m_elements[elementIndex]);
}

const Facet& Node::getFacet(unsigned int facetIndex) const noexcept
{
    return m_pMesh->getFacet(m_facets[facetIndex]);
}

bool Node::isContact() const noexcept
{
    for(std::size_t elm = 0 ; elm < getElementCount() ; ++elm)
    {
        const Element& element = getElement(elm);
        if(element.isContact())
            return true;
    }

    return false;
}

double Node::distance(const Node& n0, const Node& n1)
{
    double dist2 = 0;
    for(unsigned int i = 0 ; i < 3 ; ++i)
    {
        dist2 += (n0.getCoordinate(i) - n1.getCoordinate(i))*
                 (n0.getCoordinate(i) - n1.getCoordinate(i));
    }

    return std::sqrt(dist2);
}

void Node::updateNaturalMeshSize()
{
    double nb_adjacent_elements = double(m_elements.size());
    double element_size_sum = 0.;
    //std::cout << "nb_adjacent_elements = "<<nb_adjacent_elements << "\n";
    for (int i=0;i < nb_adjacent_elements;++i)
    {
        Element el = getElement(i);
        element_size_sum += el.getNaturalMeshSize();
    }
    //std::coutelement_size_sum
    if (nb_adjacent_elements != 0)
        m_natural_mesh_size = element_size_sum / nb_adjacent_elements;
    else
        m_natural_mesh_size = 0.;
}
double Node::getWeight(double alphaRatio, double minTargetMeshSize)
{
    //std::cout << "target_mesh_size= " << m_target_mesh_size << "\n";
    //std::cout << "natural_mesh_size= " << m_natural_mesh_size<< "\n";
    return std::pow(alphaRatio, 2.) * (std::pow(m_target_mesh_size, 2.) - std::pow(minTargetMeshSize, 2.));
}

