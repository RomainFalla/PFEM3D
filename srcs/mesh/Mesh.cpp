#include "Mesh.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include <gmsh.h>


Mesh::Mesh(const MeshCreateInfo& meshInfos) :
m_hchar(meshInfos.hchar),
m_alpha(meshInfos.alpha),
m_omega(meshInfos.omega),
m_gamma(meshInfos.gamma),
m_boundingBox(meshInfos.boundingBox),
m_useMeshRefinement(meshInfos.useMeshRefinement),
m_computeNormalCurvature(true)
{
    if (m_useMeshRefinement)
    {
        m_alphaRatio = meshInfos.refInfo.alphaRatio;
        m_gamma2 = meshInfos.refInfo.gamma2;
        m_minTargetMeshSize = meshInfos.refInfo.minTargetMeshSize;
        m_maxTargetMeshSize = meshInfos.refInfo.maxTargetMeshSize;
        m_maxProgressionFactor = meshInfos.refInfo.maxProgressionFactor;
        m_ImposedTargetMeshSizeBCtags = meshInfos.refInfo.imposedTargetMeshSizeBCtags;
        m_ImposedTargetMeshSizeValues = meshInfos.refInfo.imposedTargetMeshSizeValues;
    }
    loadFromFile(meshInfos.mshFile);
}
bool Mesh::addNodes(bool verboseOutput)
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh!");

    double m_dim_d = double(m_dim);
    double minLocalMeshSizeAfterRefinement = 0.;// for non-uniform refinement only
    if (m_useMeshRefinement)
    {
        std::size_t bulkFacetCount = m_inBulkFacetList.size();
        std::size_t boundaryFacetCount = m_facetsList.size();

        double ratio;

        /*for (std::size_t facetIndex = 0; facetIndex < bulkFacetCount; ++facetIndex) 
        {
            std::size_t indexElm1 = m_inBulkFacetList[facetIndex].m_elementIndexes[0];
            std::size_t indexElm2 = m_inBulkFacetList[facetIndex].m_elementIndexes[1];

            double L = m_inBulkFacetList[facetIndex].getLocalMeshSize();

            double meanMesure = (1. / 2.) * (m_elementsList[indexElm1].getSize() + m_elementsList[indexElm2].getSize()); // mesure means area or volume in 3D
            ratio = std::pow(meanMesure, (m_dim_d - 1.) / m_dim_d) / m_inBulkFacetList[facetIndex].getSize(); 

            minLocalMeshSizeAfterRefinement = m_refinementFactor * std::pow(L, m_dim_d);

            if (meanMesure / minLocalMeshSizeAfterRefinement < m_refinementLevels[0]) continue;

            if (meanMesure / minLocalMeshSizeAfterRefinement > m_refinementLevels[1]) continue;

            if (ratio > 1. / m_dim_d) continue;
            //it is difficult to understand why, it just results from small calculus on paper. The idea is that the surface or length of the facet must be high enough compare to the adjacent elements

            Node newNode(*this);

            newNode.m_target_mesh_size = L;
            newNode.m_natural_mesh_size = std::pow(meanMesure, 1. / m_dim_d);

            for (unsigned short k = 0; k < m_dim; ++k)
            {
                newNode.m_position[k] = 0;
                for (unsigned short d = 0; d <= m_dim-1; ++d)
                {
                    assert(m_inBulkFacetList[facetIndex].m_nodesIndexes[d] < bulkFacetCount);
                    newNode.m_position[k] += m_nodesList[m_inBulkFacetList[facetIndex].m_nodesIndexes[d]].m_position[k];
                }
                newNode.m_position[k] /= m_dim;
            }

            newNode.m_states.resize(m_nodesList[0].m_states.size());
            for (unsigned short k = 0; k < m_nodesList[0].m_states.size(); ++k)
            {
                newNode.m_states[k] = 0;
                for (unsigned short d = 0; d <= m_dim-1; ++d)
                {
                    newNode.m_states[k] += m_nodesList[m_inBulkFacetList[facetIndex].m_nodesIndexes[d]].m_states[k];
                }
                newNode.m_states[k] /= m_dim ;
            }

            m_nodesList.push_back(std::move(newNode));
        }

        for (std::size_t facetIndex = 0; facetIndex < boundaryFacetCount; ++facetIndex)
        {
            if (m_facetsList[facetIndex].isBound()) continue;

            std::size_t indexElm = m_facetsList[facetIndex].m_elementIndexes[0];

            double L = m_facetsList[facetIndex].getLocalMeshSize();

            double mesure = m_elementsList[indexElm].getSize(); // mesure means area or volume in 3D
            ratio = std::pow(mesure, (m_dim_d - 1.) / m_dim_d) / m_facetsList[facetIndex].getSize();

            minLocalMeshSizeAfterRefinement = m_refinementFactor * std::pow(L, m_dim_d);

            if (mesure / minLocalMeshSizeAfterRefinement < m_refinementLevels[0]) continue;

            if (mesure / minLocalMeshSizeAfterRefinement > m_refinementLevels[1]) continue;

            if (ratio > 1. / m_dim_d) continue;

            Node newNode(*this);

            newNode.m_target_mesh_size = L;
            newNode.m_natural_mesh_size = std::pow(mesure, 1. / m_dim_d);

            for (unsigned short k = 0; k < m_dim; ++k)
            {
                newNode.m_position[k] = 0;
                for (unsigned short d = 0; d <= m_dim-1; ++d)
                {
                    assert(m_facetsList[facetIndex].m_nodesIndexes[d] < boundaryFacetCount);
                    newNode.m_position[k] += m_nodesList[m_facetsList[facetIndex].m_nodesIndexes[d]].m_position[k];
                }
                newNode.m_position[k] /= m_dim;
            }

            newNode.m_states.resize(m_nodesList[0].m_states.size());
            for (unsigned short k = 0; k < m_nodesList[0].m_states.size(); ++k)
            {
                newNode.m_states[k] = 0;
                for (unsigned short d = 0; d <= m_dim-1; ++d)
                {
                    newNode.m_states[k] += m_nodesList[m_facetsList[facetIndex].m_nodesIndexes[d]].m_states[k];
                }
                newNode.m_states[k] /= m_dim;
            }

            newNode.m_isBound = false;
            newNode.m_isOnBoundary = true;
            newNode.m_isOnFreeSurface = true;

            m_nodesList.push_back(std::move(newNode));
        }*/
    }

    bool addedNodes = false;
    //bool refineMore = false;

    std::vector<bool> toBeDeletedElements(m_elementsList.size(), false);
    //std::vector<bool> toBeDeletedFacets(m_facetsList.size(), false);

    double limitSize = 0.;

    limitSize = m_omega * std::pow(m_hchar, m_dim);
    std::size_t elementCount = m_elementsList.size(); //avoid problems becuse the follwinĝ for loop will increase the element number

    for(std::size_t elm = 0 ; elm < elementCount ; ++elm)
    {
        double currentRefinementRatio; // only for non-uniform meshes
        double L = m_elementsList[elm].getLocalMeshSize(); // only for non-uniform meshes

        if (m_useMeshRefinement)
            minLocalMeshSizeAfterRefinement = m_refinementFactor * std::pow(L, double(m_dim));

        if (m_elementsList[elm].m_forcedRefinement) // m_forcedRefinement is (should be) always false for each element in uniform meshes.
        {
            Node newNode(*this);

            newNode.m_target_mesh_size = L;
            newNode.m_natural_mesh_size = m_elementsList[elm].getNaturalMeshSize(true);

            bool wallNode = true;
            int tag = -1;
            std::bitset<8> userDefFlags;

            for (unsigned short k = 0; k < m_dim; ++k)
            {
                newNode.m_position[k] = 0;
                for (unsigned short d = 0; d <= m_dim; ++d)
                {
                    assert(m_elementsList[elm].m_nodesIndexes[d] < m_nodesList.size());
                    if (m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].isOnBoundary())
                    {
                        newNode.m_position[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_position[k];
                        wallNode = wallNode && m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_isBound; // this will be false as soon as a boundary node is not flagged as a wall node
                        tag = m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_tag;
                        userDefFlags = m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_userDefFlags;
                    }
                }
                newNode.m_position[k] /= m_dim;
            }

            if (wallNode) continue; // we do not add a node on a solid wall for the moment due to some inconsistencies with the implementation of the initial BCs. R. Falla

            newNode.m_states.resize(m_nodesList[0].m_states.size());
            for (unsigned short k = 0; k < m_nodesList[0].m_states.size(); ++k)
            {
                newNode.m_states[k] = 0;
                for (unsigned short d = 0; d <= m_dim; ++d)
                {
                    if (m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].isOnBoundary())
                        newNode.m_states[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_states[k];
                }
                newNode.m_states[k] /= m_dim;
            }
            newNode.m_isOnBoundary = true;
            newNode.m_isOnFreeSurface = !wallNode; // true
            newNode.m_isBound = wallNode; // false
            if (wallNode)
            {
                newNode.m_tag = tag;
                newNode.m_userDefFlags = userDefFlags;
                std::cout << "tag = " << tag<<"\n" ;
            }
            m_nodesList.push_back(std::move(newNode));
            //m_boundaryInitialPos[m_nodesList.size() - 1] = m_nodesList.back().m_position;
        }

        if (m_useMeshRefinement) 
        {
            currentRefinementRatio = m_elementsList[elm].getSize() / minLocalMeshSizeAfterRefinement;
            unsigned int refinementLevel = 0;
            while (currentRefinementRatio > m_refinementLevels[refinementLevel + 1] && refinementLevel < 5)
                refinementLevel++;

            if (refinementLevel == 1)
                goto addOneNodeInTheMiddle; // get back the classic refinement used for uniform meshes, except that the if condition no more matters

            if (refinementLevel > 1)
                divideElement(m_elementsList[elm], refinementLevel);
        }
        else
        {
            //If an element is too big, we add a node at his center
            if (/**!needToRefineMore &&*/ m_elementsList[elm].getSize() > limitSize)
            {
                addOneNodeInTheMiddle:
                Node newNode(*this);

                if (m_useMeshRefinement)
                {
                    newNode.m_target_mesh_size = L;
                    newNode.m_natural_mesh_size = m_elementsList[elm].getNaturalMeshSize(true);
                }

                /*for (auto it = m_elementsList[elm].m_facets.begin(); it != m_elementsList[elm].m_facets.end(); ++it)
                {
                    Facet* facet = *it;
                    facet->m_elements.erase(&(m_elementsList[elm]));
                }*/

                for (unsigned short k = 0; k < m_dim; ++k)
                {
                    newNode.m_position[k] = 0;
                    for (unsigned short d = 0; d <= m_dim; ++d)
                    {
                        assert(m_elementsList[elm].m_nodesIndexes[d] < m_nodesList.size());
                        newNode.m_position[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_position[k];
                    }
                    newNode.m_position[k] /= (m_dim + 1);
                }

                newNode.m_states.resize(m_nodesList[0].m_states.size());
                for (unsigned short k = 0; k < m_nodesList[0].m_states.size(); ++k)
                {
                    newNode.m_states[k] = 0;
                    for (unsigned short d = 0; d <= m_dim; ++d)
                    {
                        newNode.m_states[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_states[k];
                    }
                    newNode.m_states[k] /= (m_dim + 1);
                }

                m_nodesList.push_back(std::move(newNode));

                for (unsigned short d = 0; d <= m_dim; ++d)
                {
                    Element element(*this);

                    switch (m_dim)
                    {
                    case 2:
                        if (d == m_dim)
                        {
                            element.m_nodesIndexes = { m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_nodesList.size() - 1 };
                        }
                        else
                        {
                            element.m_nodesIndexes = { m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_nodesList.size() - 1 };
                        }
                        break;

                    default:
                        if (d == m_dim)
                        {
                            element.m_nodesIndexes = { m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_elementsList[elm].m_nodesIndexes[1],
                                                      m_nodesList.size() - 1 };
                        }
                        else if (d == m_dim - 1)
                        {
                            element.m_nodesIndexes = { m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_nodesList.size() - 1 };
                        }
                        else
                        {
                            element.m_nodesIndexes = { m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_elementsList[elm].m_nodesIndexes[d + 2],
                                                      m_nodesList.size() - 1 };
                        }
                        break;
                    }

                    //No recompute of J, detJ and ivJ, everything is recomputed in triangulateAndAlphaShape
                    element.m_index = m_elementsList.size();
                    m_elementsList.push_back(std::move(element));
                    //if (element.getSize() > limitSize) refineMore = true;
                }
                toBeDeletedElements[elm] = true;

                /*if (verboseOutput)
                {
                    std::cout << "Adding node (" << "(";
                    for (unsigned short d = 0; d < m_dim; ++d)
                    {
                        std::cout << m_nodesList[m_nodesList.size() - 1].m_position[d];
                        if (d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cout << std::endl;
                }*/
                addedNodes = true;
            }
        }
        /**if (m_useMeshRefinement) // in that case the limit size depend on the local target mesh size
        {
            Facet* largestFacet;
            getLargestFacet(&m_elementsList[elm], largestFacet);

            if (largestFacet->m_elements.size() != 2) // doesn't work for boundary facets
                continue;

            if (toBeDeletedElementslargestFacet->m_elements[0]->m_index] == true || toBeDeletedElements[largestFacet->m_elements[1]->m_index] == true
                || largestFacet->m_elements[0]->getSize() > limitSize || largestFacet->m_elements[1]->getSize() > limitSize)
                continue; // check that the element considered haven't been already refined.

            double meanMesure = (largestFacet->m_elements[0]->getSize() + largestFacet->m_elements[1]->getSize()) / 2.; // this criteria lead to a smoother node density variation

            // The elements around one largest facet are too big.
            if (meanMesure > limitSize2)
            {
                std::vector<Element*> childElements;
                splitElements(largestFacet, childElements);

                for (auto it = childElements.begin(); it != childElements.end(); ++it)
                {
                    Element* el = *it;
                    el->m_index = m_elementsList.size();
                    m_elementsList.push_back(std::move(*el));
                    el->computeJ();
                    el->computeDetJ();
                    el->computeInvJ();
                    if (el->getSize() > limitSize) refineMore = true;
                }

                toBeDeletedElements[largestFacet->m_elements[0]->m_index] = true;
                toBeDeletedElements[largestFacet->m_elements[1]->m_index] = true;

                if (verboseOutput)
                {
                    std::cout << "Adding node (" << "(";
                    for (unsigned short d = 0; d < m_dim; ++d)
                    {
                        std::cout << m_nodesList[m_nodesList.size() - 1].m_position[d];
                        if (d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cout << std::endl;
                }
                addedNodes = true;
            }
        }
        needToRefineMore = refineMore;*/
    }
    

    m_elementsList.erase(std::remove_if(m_elementsList.begin(), m_elementsList.end(), [this, &toBeDeletedElements](const Element& element){
        if(toBeDeletedElements[&element - &*std::begin(m_elementsList)])
            return true;
        else
            return false;
    }), m_elementsList.end());

    return addedNodes;
}

void Mesh::divideElement(Element& el,unsigned int level) 
{
    std::cout << "Refine element with refinement level == "<< level <<"!!\n";
    if (level < 2)
    {
        std::cout<<"Warning: Divide element has no effect because the level is less than 2!!\n";
        return;
    }
    unsigned int sum;
    double level_d = double(level);

    std::vector<unsigned int> charges;
    for (unsigned int i = 0; i <= m_dim; i++)
    {
        if (el.m_chargedNode[i])
            charges.push_back(1);
        else
            charges.push_back(0);
    }

    const Node& n0 = getNode(el.m_nodesIndexes[0]);
    const Node& n1 = getNode(el.m_nodesIndexes[1]);
    const Node& n2 = getNode(el.m_nodesIndexes[2]);

    double f0, f1, f2;

    if (m_dim == 2)
    {     
        for (unsigned int i = charges[0]; i < level; i++)
        {
            f0 = double(i) / level_d;
            for (unsigned int j = charges[1]; j < level; j++)
            {
                sum = i+j;
                if (sum > level) break;
                f1 = double(j) / level_d;
                for (unsigned int k = charges[2]; k < level; k++)
                {
                    sum = i+j+k;
                    if (sum > level) break;
                    f2 = double(k) / level_d;

                    if (sum != level) continue;

                    Node newNode(*this);
                    newNode.m_target_mesh_size = el.getLocalMeshSize();
                    newNode.m_natural_mesh_size = el.getNaturalMeshSize(true);

                    for (unsigned short m = 0; m < m_dim; ++m)
                        newNode.m_position[m] = f0 * n0.m_position[m] + f1 * n1.m_position[m] + f2 * n2.m_position[m];


                    newNode.m_states.resize(m_nodesList[0].m_states.size());
                    for (unsigned short m = 0; m < m_nodesList[0].m_states.size(); ++m)
                        newNode.m_states[m] = f0 * n0.m_states[m] + f1 * n1.m_states[m] + f2 * n2.m_states[m];

                    m_nodesList.push_back(std::move(newNode));
                }

            }
        }
    }
    else //dim == 3
    {
        const Node& n3 = getNode(el.m_nodesIndexes[3]);
        double f3;

        for (unsigned int i = charges[0]; i < level; i++)
        {
            f0 = double(i) / level_d;
            for (unsigned int j = charges[1]; j < level; j++)
            {
                sum = i+j;
                if (sum > level) break;
                f1 = double(j) / level_d;
                for (unsigned int k = charges[2]; k < level; k++)
                {
                    sum = i+j+k;
                    if (sum> level) break;
                    f2 = double(k) / level_d;
                    for (unsigned int l = charges[3]; l < level; l++)
                    {
                        sum = i+j+k+l;
                        if (sum > level) break;
                        f3 = double(l) / level_d;

                        if (sum != level) continue;

                        Node newNode(*this);
                        newNode.m_target_mesh_size = el.getLocalMeshSize();
                        newNode.m_natural_mesh_size = el.getNaturalMeshSize(true);

                        for (unsigned short m = 0; m < m_dim; ++m)
                            newNode.m_position[m] = f0 * n0.m_position[m] + f1 * n1.m_position[m] + f2 * n2.m_position[m] + f3 * n3.m_position[m];


                        newNode.m_states.resize(m_nodesList[0].m_states.size());
                        for (unsigned short m = 0; m < m_nodesList[0].m_states.size(); ++m)
                            newNode.m_states[m] = f0 * n0.m_states[m] + f1 * n1.m_states[m] + f2 * n2.m_states[m] + f3 * n3.m_states[m];

                        m_nodesList.push_back(std::move(newNode));
                    }
                }

            }
        }
     
    }

    //loop to avoid the adjacent elements to add nodes on the same facets
    for (auto it = el.m_facets.begin(); it != el.m_facets.end(); ++it)
    {
        std::size_t i = *it;
        if (i > m_inBulkFacetList.size())
            continue;

        std::size_t otherElemIndex = m_inBulkFacetList[i].m_elementIndexes[0];
        std::size_t otherNodeIndex = m_inBulkFacetList[i].m_outNodeIndexes[0];

        if (otherElemIndex == el.m_index) 
        {
            otherElemIndex = m_inBulkFacetList[i].m_elementIndexes[1];
            otherNodeIndex = m_inBulkFacetList[i].m_outNodeIndexes[1];
        }

        std::size_t inc = 0;
        for (std::size_t j : m_elementsList[otherElemIndex].m_nodesIndexes)
        {
            if (j == otherNodeIndex) 
            {
                m_elementsList[otherElemIndex].m_chargedNode[inc] = true;
            }
            inc++;
        }
    }
}

bool Mesh::checkBoundingBox(bool verboseOutput) noexcept
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    std::vector<bool> toBeDeletedElement(m_elementsList.size(), false);   //Should the element be deleted
    std::vector<bool> toBeDeletedNodes(m_nodesList.size(), false);   //Should the node be deleted

    std::vector<std::size_t> nodesIndexesDeleted = {};
    std::vector<std::size_t> elementIndexesDeleted = {};

    bool outofBBNodes = false;

    auto isNodeOutOfBB = [this](const Node& node) -> bool {
        for(unsigned short d = 0 ; d < m_dim ; ++d)
        {
            if(node.m_position[d] < m_boundingBox[d] ||
               node.m_position[d] > m_boundingBox[d + m_dim])
            {
                return true;
            }
        }

        return false;
    };

    //If the whole element is out of the bounding box, we delete it.
    // Bounding box fromat: [xmin, ymin, zmin, xmax, ymax, zmax]
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(isNodeOutOfBB(m_nodesList[n]))
        {
            toBeDeletedNodes[n] = true;

            if(!m_nodesList[n].isFree())
            {
                //There are elements to delete
                for(std::size_t i = 0 ; i < m_nodesList[n].m_elements.size() ; ++i)
                {
                    toBeDeletedElement[m_nodesList[n].m_elements[i]] = true;
                }
            }
        }
    }

    for(std::size_t n = 0 ; n < toBeDeletedNodes.size() ; ++n)
        if(toBeDeletedNodes[n]) nodesIndexesDeleted.push_back(n);

    for(std::size_t elm = 0 ; elm < toBeDeletedElement.size() ; ++elm)
        if(toBeDeletedElement[elm]) elementIndexesDeleted.push_back(elm);

    m_elementsList.erase(
    std::remove_if(m_elementsList.begin(), m_elementsList.end(), [this, &toBeDeletedElement, verboseOutput](const Element& element){
        if(toBeDeletedElement[&element - &*std::begin(m_elementsList)])
        {
            /*if(verboseOutput)
            {
                std::cout << "Removing out of bounding box element: [";
                for(unsigned short n= 0 ; n < m_dim + 1 ; ++n)
                {
                    std::cout << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cout << this->m_nodesList[element.m_nodesIndexes[n]].m_position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    if(n == m_dim)
                            std::cout << "]";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }*/

            return true;
        }
        else
            return false;
    }), m_elementsList.end());

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this, &toBeDeletedNodes, verboseOutput](const Node& node)
    {
        if(toBeDeletedNodes[&node - &*std::begin(m_nodesList)])
        {
            /*if(verboseOutput)
            {
                std::cout << "Removing out of bounding box node (";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }*/

            return true;
        }
        else
           return false;
    }), m_nodesList.end());

    for(std::size_t i = 0 ; i < nodesIndexesDeleted.size() ; ++i)
    {
        for(Element& element : m_elementsList)
        {
            for(std::size_t& n : element.m_nodesIndexes)
            {
                if(n > nodesIndexesDeleted[i])
                    n--;
            }
        }

        for(Node& node : m_nodesList)
        {
            for(std::size_t& n : node.m_neighbourNodes)
            {
                if(n > nodesIndexesDeleted[i])
                    n--;
            }
        }

        for(std::size_t j = i + 1 ; j < nodesIndexesDeleted.size() ; ++j)
            nodesIndexesDeleted[j]--;
    }

    for(std::size_t i = 0 ; i < elementIndexesDeleted.size() ; ++i)
    {
        for(Node& node : m_nodesList)
        {
            for(std::size_t& elm : node.m_elements)
            {
                if(elm > elementIndexesDeleted[i])
                    elm--;
            }
        }

        for(std::size_t j = i + 1 ; j < elementIndexesDeleted.size() ; ++j)
            elementIndexesDeleted[j]--;
    }

    return outofBBNodes;
}

void Mesh::smoothingLcLoc()
{
    std::vector<Node*> seedingNodes;
    std::vector<Node*> nodesToCheckNextIteration;
    std::vector<Node*> freshlyRefinedNodes;

    double TargetMeshSizeRatio;
    double TargetMeshVsNeighboursRatio;
    double realProg;

    // loop to find the nodes which have a too high LcLocRatio and "seed" a
    // refinement front from there to propagate until smooth
    for (auto it = m_nodesList.begin(); it != m_nodesList.end(); ++it)
    {
        Node* nod1 = &(*it);
        for (auto it2 = nod1->m_neighbourNodes.begin(); it2 != nod1->m_neighbourNodes.end(); ++it2)
        {
            Node* nod2 = &(m_nodesList[*it2]);
            TargetMeshSizeRatio = nod2->m_target_mesh_size/ nod1->m_target_mesh_size; // = (localProgression + 1) / 2 = (LcOtherNode/LcThisNode + 1) / 2
            TargetMeshVsNeighboursRatio = nod1->m_natural_mesh_size / nod1->m_target_mesh_size;
            realProg = 1. + (m_maxProgressionFactor - 1.) * TargetMeshVsNeighboursRatio;
            //f = (realProg + 1.) / 2.; // = max allowed edgeLcLoc/nodeLcLoc
            //more complex but more rigorous, take into account the real size of the elements to update the target mesh size of the neighbouring elements

            if (TargetMeshSizeRatio > realProg)
            {
                seedingNodes.push_back(nod1);
                break;
            }
        }
    }

    //std::cout << "seedingNodes size = " << seedingNodes.size() << "\n";

    // use nod->propagate() on these nodes, to propagate a "wave front" of
    // refinement over sevral neighboring nodes. Each node that is checked and
    // if found to need refinement will trigger its neighboring nodes to be
    // checked also and so on, until there are no more harsh increases of LcLoc
    int inc = 0;
    while (true) // we are sure to go out because we keep increasing the LcLoc to check with and the LcLoc are bounded. Therefore, seedingNodes or nodesToCheckNextIteration must be empty at some point.
    {
        for (auto it = seedingNodes.begin(); it != seedingNodes.end(); ++it)
        {
            TargetMeshVsNeighboursRatio = (*it)->m_natural_mesh_size / (*it)->m_target_mesh_size;
            realProg = 1. + (m_maxProgressionFactor - 1.) * TargetMeshVsNeighboursRatio;
            freshlyRefinedNodes.clear();
            freshlyRefinedNodes = (*it)->propagate(realProg * (*it)->m_target_mesh_size);
            for (auto it2 = freshlyRefinedNodes.begin(); it2 != freshlyRefinedNodes.end(); ++it2)
            {
                Node* nod = *it2;
                nodesToCheckNextIteration.push_back(nod); // add all neighbouring nodes that were just refined
            }
            //std::cout << "realProg = " << realProg << "\n";
            //std::cout << "TargetMeshVsNeighboursRatio  = " << TargetMeshVsNeighboursRatio << "\n";
        }
        seedingNodes.clear();
        if (nodesToCheckNextIteration.empty()) // no more nodes were refined in this iteration -> the "wave front" ends
            break;

        // prepare next iteration
        seedingNodes = nodesToCheckNextIteration;
        nodesToCheckNextIteration.clear();

        inc++;
        if (inc == 100)
        {
            throw std::runtime_error("while loop in the smoothLcLoc function goes forever...");
        }
    }
    //throw std::runtime_error("I want the program to stop here 4! ");
    //std::cout << "number of iteration in smoothing = " << inc << " \n";
}


/**void Mesh::linkElementsAndFacets()
{
    for (auto it = m_facetsList.begin(); it != m_facetsList.end(); ++it)
    {
        Facet* f = &(*it);
        f->m_elements.clear();
        if (m_dim == 2)
        {
            Node n1 = m_nodesList[f->m_nodesIndexes[0]];
            Node n2 = m_nodesList[f->m_nodesIndexes[1]];
            // find all triangles that both nodes have in common
            for (std::size_t i = 0; i != m_dim; ++i)
            {
                for (std::size_t j = 0; j != m_dim; ++j)
                {
                    if (n1.m_elements[i] == n2.m_elements[j]) // found a common triangle
                    {
                        Element* elm = &(m_elementsList[n1.m_elements[i]]);
                        f->m_elements.push_back(elm);
                        elm->m_facets.push_back(f);
                    }
                }
            }
        }
        else 
        {
            Node n1 = m_nodesList[f->m_nodesIndexes[0]];
            Node n2 = m_nodesList[f->m_nodesIndexes[1]];
            Node n3 = m_nodesList[f->m_nodesIndexes[2]];

            // much time consuming than for 2D... this should be parallelized at least
            for (std::size_t i = 0; i != m_dim; ++i)
            {
                for (std::size_t j = 0; j != m_dim; ++j)
                {
                    if (n1.m_elements[i] != n2.m_elements[j]) continue;

                    for (std::size_t k = 0; k != m_dim; ++k)
                    {
                        if (n1.m_elements[i] == n3.m_elements[k]) // found a common tetrahedron
                        {
                            Element* elm = &(m_elementsList[n1.m_elements[i]]);
                            f->m_elements.push_back(elm);
                            elm->m_facets.push_back(f);
                        }
                    }
                }
            }
        }
    }
}*/

void Mesh::computeMeshDim()
{
    int elementDim = -1;

    // loop over the dimension i to get the maximum element dimension in the mesh
    for(unsigned short i = 2 ; i <= 3 ; ++i)
    {
        std::vector<int> eleTypes;
        gmsh::model::mesh::getElementTypes(eleTypes, i);

        switch(eleTypes.size())
        {
            case 0:
                break;
            case 1:
                elementDim = i;
                break;
            default:
                throw std::runtime_error("do not use hybrid meshes for PFEM simulations!");
        }
    }

    if(elementDim == -1)
        throw std::runtime_error("there is no suitable elements in the .msh file!");
    else
        m_dim = static_cast<unsigned short>(elementDim);
}

void Mesh::computeFSNormalCurvature()
{
    m_freeSurfaceCurvature.clear();
    m_boundFSNormal.clear();

    switch(m_dim)
    {
        case 2:
            computeFSNormalCurvature2D();
            break;
        default:

            computeFSNormalCurvature3D();
            break;
    }
}

void Mesh::displayToConsole() const noexcept
{
    std::cout << "Mesh dimension: " << m_dim << "D\n";
    std::cout << "hchar: " << m_hchar << "\n";
    std::cout << "alpha: " << m_alpha << "\n";
    std::cout << "omega: " << m_omega << "\n";
    std::cout << "gamma: " << m_gamma << std::endl;
}

std::vector<std::array<double, 3>> Mesh::getGaussPoints(unsigned int dimension, unsigned int n) const
{
    switch(dimension)
    {
        case 1:
        {
            switch(n)
            {
                case 1:
                    return { {0.0, 0.0, 0.0} };
                case 2:
                    return { {1.0/std::sqrt(3.0), 0, 0},
                             {-1.0/std::sqrt(3.0), 0, 0} };
                case 3:
                    return { {-std::sqrt(3.0/5.0), 0, 0},
                             {0, 0, 0},
                             {std::sqrt(3.0/5.0), 0, 0} };
                case 4:
                    return { {-0.861136311594053, 0, 0},
                             {-0.339981043584856, 0, 0},
                             {0.339981043584856, 0, 0},
                             {0.861136311594053, 0, 0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 1D: " + std::to_string(n));
            }
        }

        case 2:
        {
            switch(n)
            {
                case 1:
                    return { {1.0/3.0, 1.0/3.0, 0.0} };
                case 3:
                    return { {1.0/6.0, 1.0/6.0, 0.0},
                             {1.0/6.0, 2.0/3.0, 0.0},
                             {2.0/3.0, 1.0/6.0, 0.0} };
                case 4:
                    return { {1.0/3.0, 1.0/3.0, 0.0},
                             {0.6, 0.2, 0.0},
                             {0.2, 0.2, 0.0},
                             {0.2, 0.6, 0.0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 2D: " + std::to_string(n));
            }
        }

        case 3:
        {
            switch(n)
            {
                case 1:
                    return { {1.0/4.0, 1.0/4.0, 1.0/4.0} };
                case 4:
                    return{ {0.585410196624968, 0.138196601125011, 0.138196601125011},
                            {0.138196601125011, 0.585410196624968, 0.138196601125011},
                            {0.138196601125011, 0.138196601125011, 0.585410196624968},
                            {0.138196601125011, 0.138196601125011, 0.138196601125011} };
                case 5:
                    return{ {1.0/4.0, 1.0/4.0, 1.0/4.0},
                            {1.0/2.0, 1.0/6.0, 1.0/6.0},
                            {1.0/6.0, 1.0/2.0, 1.0/6.0},
                            {1.0/6.0, 1.0/6.0, 1.0/2.0},
                            {1.0/6.0, 1.0/6.0, 1.0/6.0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 3D: " + std::to_string(n));
            }
        }

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

std::vector<double> Mesh::getGaussWeight(unsigned int dimension, unsigned int n) const
{
    switch(dimension)
    {
        case 1:
        {
            switch(n)
            {
                case 1:
                    return {2.0};
                case 2:
                    return {1.0, 1.0};
                case 3:
                    return {5.0/9.0, 8.0/9.0, 5.0/9.0};
                case 4:
                    return {0.347853845137454,0.652145154862546 ,0.652145154862546, 0.347853845137454};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 1D: " + std::to_string(n));
            }
        }

        case 2:
        {
            switch(n)
            {
                case 1:
                    return {1.0};
                case 3:
                    return {1.0/3.0, 1.0/3.0, 1.0/3.0};
                case 4:
                    return {-0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 2D: " + std::to_string(n));
            }
        }

        case 3:
        {
            switch(n)
            {
                case 1:
                    return {1.0};
                case 4:
                    return {0.25, 0.25, 0.25, 0.25};
                case 5:
                    return {-0.8, 0.45, 0.45, 0.45, 0.45};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 3D: " + std::to_string(n));
            }
        }

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

double Mesh::getRefElementSize(unsigned int dimension) const
{
    switch(dimension)
    {
        case 1:
            return 2.0;

        case 2:
            return 0.5;

        case 3:
            return 0.16666666666666666666666666666667;

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

void Mesh::updateBoundingElementSize()
{
    double minValue = std::numeric_limits<double>::infinity();
    double maxValue = 0;
    for (auto it = m_elementsList.begin(); it != m_elementsList.end(); ++it) 
    {
        Element el = *it;
        double L = el.getNaturalMeshSize(false);
        if (L > maxValue) maxValue = L;
        if (L < minValue) minValue = L;
    }
    m_minTargetMeshSize = minValue;
    m_maxTargetMeshSize = maxValue;
}

void Mesh::updateLocalElementSize()
{
    if (m_ImposedTargetMeshSizeBCtags.size() == 0)
    {
        // default mesh refinement
        for (auto it = m_nodesList.begin(); it != m_nodesList.end(); ++it)
        {
            it->updateNaturalMeshSize();
            double L = it->getNaturalMeshSize();
            it->setLocalMeshSize(L);
            /** set the local mesh size to the natural value defined by the local size of the mesh
            (the local mesh size is not imposed yet), it is a first step to be able to run non-uniform simulations */
        }
    }
    else 
    {
        for (auto it = m_nodesList.begin(); it != m_nodesList.end(); ++it)
        {
            it->updateNaturalMeshSize();
            it->setLocalMeshSize(m_maxTargetMeshSize);
            std::size_t inc = 0;
            for (auto it2 = m_ImposedTargetMeshSizeBCtags.begin(); it2 != m_ImposedTargetMeshSizeBCtags.end(); ++it2)
            {
                if (it->m_tag == *it2)
                    it->setLocalMeshSize(m_ImposedTargetMeshSizeValues[inc]);
                inc++;
            }
        }
        smoothingLcLoc();
    }
    //throw std::runtime_error("I want the program to stop here 2! ");
}

std::vector<std::vector<double>> Mesh::getShapeFunctions(unsigned int dimension, unsigned int n) const
{
    std::vector<std::array<double, 3>> gps = getGaussPoints(dimension, n);

    std::vector<std::vector<double>> sfs(gps.size());

    unsigned int counter = 0;
    for(std::array<double, 3> gp : gps)
    {
        std::vector<double> sf(dimension + 1);

        switch(dimension)
        {
            case 1:
                sf[0] = (1 - gp[0])/2;
                sf[1] = (1 + gp[0])/2;
                break;

            case 2:
                sf[0] = 1 - gp[0] - gp[1];
                sf[1] = gp[0];
                sf[2] = gp[1];
                break;

            case 3:
                sf[0] = 1 - gp[0] - gp[1] - gp[2];
                sf[1] = gp[0];
                sf[2] = gp[1];
                sf[3] = gp[2];
                break;

            default:
                throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
        }

        sfs[counter] = std::move(sf);
        counter++;
    }

    return sfs;
}

std::vector<std::vector<double>> Mesh::getGradShapeFunctions(unsigned int dimension) const
{
    std::vector<std::vector<double>> gradsfs(dimension);

    for(unsigned int i = 0 ; i < gradsfs.size() ; ++i)
        gradsfs[i].resize(dimension + 1);

    switch(dimension)
    {
        case 1:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            break;

        case 2:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            gradsfs[0][2] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 1;
            break;

        case 3:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            gradsfs[0][2] = 0;
            gradsfs[0][3] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 1;
            gradsfs[1][3] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 0;
            gradsfs[1][3] = 1;
            break;

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }

    return gradsfs;
}

void Mesh::loadFromFile(const std::string& fileName)
{
    m_nodesList.clear();

    std::map<std::size_t, std::size_t> nodeIndexMap;

    gmsh::initialize();
#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif

    std::ifstream file(fileName);
    if(file.is_open())
        file.close();
    else
        throw std::runtime_error("the input .msh file does not exist!");

    gmsh::open(fileName);

    // Check that the mesh is not 3D
    computeMeshDim();

    if (m_useMeshRefinement) 
    {
        if (m_dim == 2)
        {
            m_refinementFactor = 2. / 3.;
            m_refinementLevels[0] = 2.;
            m_refinementLevels[1] = 3.;
        }
        else
        {
            m_refinementFactor = 1. / 2.;
            m_refinementLevels[0] = 3.;
            m_refinementLevels[1] = 4.;
        }

        for (int k = 2; k <= 5; k++) 
        {
            m_refinementLevels[k] = std::pow(k, double(m_dim));
        }
    }


    if(m_boundingBox.size() != 2*m_dim)
        throw std::runtime_error("bad bounding box format! Format: [xmin, ymin, xmax, ymax]");

    // We retrieve the tags of the physical groups of dimension m_dim and
    // m_dim-1
    std::vector<std::pair<int, int>> physGroupHandlesHD;
    gmsh::model::getPhysicalGroups(physGroupHandlesHD, m_dim);

    std::vector<std::pair<int, int>> physGroupHandlesLD;
    gmsh::model::getPhysicalGroups(physGroupHandlesLD, m_dim - 1);

    //The file should contain physical group for the boundary, the free surface and the fluid

    for(auto physGroupLD : physGroupHandlesLD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim - 1, physGroupLD.second, name);

        std::vector<double> coord;
        std::vector<std::size_t> dummyNodesTagsBoundary;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_dim - 1, physGroupLD.second,
                                                    dummyNodesTagsBoundary, coord);

        for(std::size_t i = 0 ; i < dummyNodesTagsBoundary.size() ; ++i)
        {
            Node node(*this);
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                node.m_position[d] = coord[3*i + d];
            }

            //If the nodes is already on the boundary, we do not add it twice
            if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                continue;

            node.m_isOnBoundary = true;
            if (name != "FreeSurface")
            {
                node.m_isBound = true;

                auto posBCinTagNames = std::find(m_tagNames.begin(), m_tagNames.end(), name);
                if (posBCinTagNames == std::end(m_tagNames))
                {
                    m_tagNames.push_back(name);
                    node.m_tag = static_cast<int>(m_tagNames.size() - 1);
                }
                else
                {
                    node.m_tag = static_cast<int>(std::distance(m_tagNames.begin(), posBCinTagNames));
                }
            }
            nodeIndexMap[dummyNodesTagsBoundary[i]] = m_nodesList.size();
            m_nodesList.push_back(std::move(node));
            m_boundaryInitialPos[m_nodesList.size() - 1] = m_nodesList.back().m_position;
        }

    }

    for(auto physGroupHD : physGroupHandlesHD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim, physGroupHD.second, name);

        std::vector<std::size_t> dummyNodesTags;
        std::vector<double> coord;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_dim, physGroupHD.second,
                                                    dummyNodesTags, coord);

        for(std::size_t i = 0 ; i < dummyNodesTags.size() ; ++i)
        {
            Node node(*this);
            for(unsigned short d = 0 ; d < m_dim ; ++d)
                node.m_position[d] = coord[3*i + d];

            //If the nodes is already on the boundary, we do not add it twice
            if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                continue;

            node.m_isBound = false;

            auto posBCinTagNames = std::find(m_tagNames.begin(), m_tagNames.end(), name);
            if(posBCinTagNames == std::end(m_tagNames))
            {
                m_tagNames.push_back(name);
                node.m_tag = m_tagNames.size() - 1;
            }
            else
            {
                node.m_tag = static_cast<uint16_t>(std::distance(m_tagNames.begin(), posBCinTagNames));
            }

            nodeIndexMap[dummyNodesTags[i]] = m_nodesList.size();
            m_nodesList.push_back(std::move(node));
        }
    }

//#ifndef NDEBUG
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        for(std::size_t n2 = 0 ; n2 < m_nodesList.size() ; ++n2)
        {
            if(n != n2)
            {
                if(m_nodesList[n] == m_nodesList[n2])
                {
                    std::cerr << "Duplicate node found: " << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cerr << m_nodesList[n].m_position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cerr << std::endl;
                    throw std::runtime_error("duplicates nodes!");
                }
            }
        }
    }
    std::cout << "List of founded BC names: " << std::endl;
    for (auto& aString : m_tagNames)
    {
        std::cout << aString << std::endl;
    }
//#endif

    if(m_nodesList.empty())
        throw std::runtime_error("no nodes loaded! Did you add the physical groups in the .geo file?");

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;

    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags);
    //std::cout << "elementType = " << elementTypes[0] << "\n";
    std::size_t nb_types = elementTypes.size();
    //std::cout << "elementType = " << elementTypes[0] << "\n";
    m_elementsList.clear();

    if (m_dim == 2)
    {
        std::vector<std::size_t> triangles;
        for (std::size_t i = 0; i < nb_types; i++)
        {
            //std::cout << "elementType = " << elementTypes[i] << "\n";
            if (elementTypes[i] != 2) continue;//not triangles

            triangles = nodeTags[i];
            break;
        }
        std::size_t L = triangles.size() / 3;
        //std::cout << "triangles size = " << triangles.size() << "\n";

        for (std::size_t i = 0; i < L; ++i)
        {
            Element element(*this);
            std::vector<std::size_t> elementNodes = { nodeIndexMap[triangles[3 * i]], nodeIndexMap[triangles[3 * i + 1]], nodeIndexMap[triangles[3 * i + 2]]};
            element.m_nodesIndexes = elementNodes;
            element.m_index = m_elementsList.size();
            element.build(elementNodes, m_nodesList);
            //std::cout << element.getSize() << "\n";
            m_elementsList.push_back(element);
        }
    }
    else
    {
        std::vector<std::size_t> tetrahedrons;
        for (std::size_t i = 0; i < nb_types; i++)
        {
            if (elementTypes[i] != 4) continue;//not tetrahedrons

            tetrahedrons = nodeTags[i];
            break;
        }
        std::size_t L = tetrahedrons.size() / 4;

        for (std::size_t i = 0; i < L; ++i)
        {
            Element element(*this);
            std::vector<std::size_t> elementNodes = { nodeIndexMap[tetrahedrons[3 * i]],nodeIndexMap[tetrahedrons[3 * i + 1]], nodeIndexMap[tetrahedrons[3 * i + 2]], nodeIndexMap[tetrahedrons[3 * i + 3]]};
            element.m_nodesIndexes = elementNodes;
            element.build(elementNodes, m_nodesList);
            element.m_index = m_elementsList.size();
            m_elementsList.push_back(std::move(element));

        }
    }
    std::cout << "m_elementsList size = " << m_elementsList.size() << "\n";


    gmsh::finalize();

    //double meanElementSize = 0.;

    if (m_useMeshRefinement)
    {
        std::cout << "updateLocalMeshSize!\n";
        for (auto it = m_elementsList.begin(); it != m_elementsList.end(); ++it)
        {
            Element el = *it;
            //meanElementSize += el.getSize();
            std::size_t elem_index = el.m_index;
            std::size_t L = el.m_nodesIndexes.size();
            for (std::size_t k = 0; k < L; k++)
            {
                std::size_t node_index = el.m_nodesIndexes[k];
                m_nodesList[node_index].m_elements.push_back(elem_index);
            }
        }
        updateBoundingElementSize();
        updateLocalElementSize();
        //meanElementSize /= m_elementsList.size();
        //std::cout << "meanElementSize = " << meanElementSize << "\n";
    }
    triangulateAlphaShape();
}

void Mesh::remesh(bool verboseOutput)
{
    checkBoundingBox(verboseOutput);
    addNodes(verboseOutput);
    removeNodes(verboseOutput);
    triangulateAlphaShape();
}

bool Mesh::removeNodes(bool verboseOutput) noexcept
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    std::vector<bool> touched(m_nodesList.size(), false);       //Is the node next to a node which should be deleted
    std::vector<bool> toBeDeleted(m_nodesList.size(), false);   //Should the node be deleted

    bool removeNodes = false;

    double limitLength = m_gamma*m_hchar;

    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        //We the node is free, it is not in an element; if a node is already tagged,
        //we can skip it (do not delete node twice)
        if(toBeDeleted[i] || m_nodesList[i].isFree())
            continue;

        double tmeshSize = m_nodesList[i].getLocalMeshSize();

        for(unsigned int j = 0 ; j < m_nodesList[i].m_neighbourNodes.size() ; ++j)
        {
            double d = 0;

            double tmeshSizeNeighbour = m_nodesList[m_nodesList[i].m_neighbourNodes[j]].getLocalMeshSize();

            for(unsigned short k = 0 ; k < m_dim ; ++k)
            {
                d +=(m_nodesList[i].m_position[k] - m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_position[k])
                   *(m_nodesList[i].m_position[k] - m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_position[k]);
            }
            d = std::sqrt(d);

            if (m_useMeshRefinement) // if the mesh is non-uniform, the limit length depend on the local mesh size... -R.Falla
            {
                limitLength = m_gamma2 * (tmeshSize + tmeshSizeNeighbour) / 2.; // I don't know if the threshold is well adapted for 3D simulations -R.Falla
            }

            //Two nodes are too close.
            if(d <= limitLength)
            {
                //If the neighbour nodes is touched, we delete the current nodes
                if(touched[m_nodesList[i].m_neighbourNodes[j]])
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[i].m_isBound || m_nodesList[i].m_isOnFreeSurface)
                        continue;

                    toBeDeleted[i] = true;
                    removeNodes = true;
                }
                //If the neighbour nodes is not touched, we delete the neoghbour nodes
                else
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_isBound ||
                       m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_isOnFreeSurface)
                        continue;

                    touched[i] = true;
                    toBeDeleted[m_nodesList[i].m_neighbourNodes[j]] = true;
                    removeNodes = true;
                }
            }
        }
    }

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this, &toBeDeleted, verboseOutput](const Node& node)
    {
       if(toBeDeleted[&node - &*std::begin(m_nodesList)])
       {
           /*if(verboseOutput)
           {
                std::cout << "Removing node " << "(";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
           }*/ 

           return true;
       }
       else
           return false;
    }), m_nodesList.end());

    return removeNodes;
}

void Mesh::restoreNodesList()
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("the nodes list was not saved before or does not exist!");

    m_nodesList = std::move(m_nodesListSave);

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::saveNodesList()
{
    if(m_nodesList.empty())
        throw std::runtime_error("the nodes list does not exist!");

    m_nodesListSave = m_nodesList;
}

void Mesh::triangulateAlphaShape()
{
    if(!m_useMeshRefinement)
    {
        if (m_dim == 2)
            triangulateAlphaShape2D();
        else
            triangulateAlphaShape3D();
    }
    else 
    {
        if (m_dim == 2)
            TriangulateWeightedAlphaShape2D();
        else
            TriangulateWeightedAlphaShape3D();
    }
}

void Mesh::updateNodesPosition(std::vector<double> deltaPos)
{
    if(deltaPos.size() != m_nodesList.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector");

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] += deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::updateNodesPositionFromSave(std::vector<double> deltaPos)
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("you did not save the nodes list!");
    else if(deltaPos.size() != m_nodesListSave.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector");

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] = m_nodesListSave[n].m_position[d] + deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

/**void Mesh::splitElements(Facet* f, std::vector<Element*>& childElements)
{
    Element* elm1 = f->m_elements[0];
    Element* elm2 = f->m_elements[1];

    std::size_t indexOppositeNode1 = getOutNodeIndex(f, elm1);
    std::size_t indexOppositeNode2 = getOutNodeIndex(f, elm2);

    Node newNode(*this);
    newNode = splitElement(elm1, indexOppositeNode1, childElements);
    newNode = splitElement(elm2, indexOppositeNode2, childElements);
    m_nodesList.push_back(newNode);
}*/

std::size_t Mesh::getOutNodeIndex(Facet* face, Element* element)
{
    // will only give consistant results if *element has *face as a facet!!! 
    unsigned short  i = 0;
    for (i = 0; i < m_dim+1; ++i)
    {
        bool found = true;
        for (unsigned short j = 0; j < m_dim; ++j)
        {
            if (element->m_nodesIndexes[i] == face->m_nodesIndexes[j])
            {
                found = false;
                break;
            }
        }
        if (found == true) break;
    }
    return element->m_nodesIndexes[i];
}
/*void Mesh::getLargestFacet(Element* element, Facet* &face )
{
    double max_size = 0.;
    double dummy_size;
    for (auto it = element->m_facets.begin(); it != element->m_facets.end(); ++it) 
    {
        Facet* f = *it;
        dummy_size = f->getDetJ();
        if (dummy_size > max_size)
        {
            max_size = dummy_size;
            face = f; // update the pointer passed by adress;
        }
    }
}*/

ELEMENT_TYPE Mesh::getElementType(std::vector<std::size_t> nodesIndexes)
{
    std::size_t inc = 0;
    std::size_t L = nodesIndexes.size();
    for (std::size_t i = 0; i < L; i++)
    {
        if (m_nodesList[nodesIndexes[i]].isOnBoundary())
            inc++;
    }

    if (m_dim == 2)
    {  
        if (inc < 2)
            return IN_BULK;
        else if (inc == 2)
            return INSIDE_BOUNDARY;
        else
            return OUTSIDE_BOUNDARY;
    }
    else //m_dim==3
    {
        if (inc < 3)
            return IN_BULK;
        else if (inc == 3)
            return INSIDE_BOUNDARY;
        else
            return OUTSIDE_BOUNDARY;
    }
    return IN_BULK;
}

double Mesh::getCircumScribedRadius(std::vector<size_t> nods)
{
    double x1 = m_nodesList[nods[0]].m_position[0];
    double x2 = m_nodesList[nods[1]].m_position[0];
    double x3 = m_nodesList[nods[2]].m_position[0];
    double y1 = m_nodesList[nods[0]].m_position[1];
    double y2 = m_nodesList[nods[1]].m_position[1];
    double y3 = m_nodesList[nods[2]].m_position[1];
    double x21 = x2 - x1;
    double y21 = y2 - y1;
    double x31 = x3 - x1;
    double y31 = y3 - y1;
    double DET = 0.5 / (x21 * y31 - x31 * y21);
    double R21 = (x21 * x21) + (y21 * y21);
    double R31 = (x31 * x31) + (y31 * y31);
    double xdummy = DET * (R21 * y31 - R31 * y21);
    double ydummy = DET * (x21 * R31 - x31 * R21);

    return std::sqrt(xdummy * xdummy + ydummy * ydummy);
}
