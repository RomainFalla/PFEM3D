#include "Node.hpp"

#include <cmath>
#include "Mesh.hpp"
#include "Element.hpp"
#include "Facet.hpp"

Node::Node(Mesh& mesh):
m_pMesh(&mesh)
{

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
