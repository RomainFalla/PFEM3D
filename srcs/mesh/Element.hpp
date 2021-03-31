#pragma once
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <array>
#include <cstdint>
#include <vector>
#include <set>

#include "mesh_defines.h"

class Mesh;
class Node;
class Facet;

/// \enumeration to identify the type of element

/**
 * \class Element
 * \brief Represents an element in the PFEM.
 */
class MESH_API Element
{
    public:
        /// \param mesh a reference to the mesh
        Element(Mesh& mesh);
        Element(const Element& element)             = default;
        Element(Element&& element)                  = default;
        ~Element()                                  = default;

        /// \param nodeElmIndex The index of the node inside the element.
        /// \return A constant reference to the node.
        const Node& getNode(unsigned int nodeElmIndex) const noexcept;

        /// \param nodeElmIndex The index of the node inside the element.
        /// \return The index of the node inside the mesh.
        inline std::size_t getNodeIndex(unsigned int nodeElmIndex) const noexcept;

        /// \param stateIndex The index of the state.
        /// \return A vector containing the considered state for each node of the element.
        std::vector<double> getState(unsigned int stateIndex) const noexcept;

        /// \return The determinant of the Jacobian matrix of the change of variable to the reference space.
        inline double getDetJ() const noexcept;

        /**
         * \brief Get an element of the Jacobian matrix of the change of variable to the reference space:
         * \f{eqnarray*}{
         *    &\text{2D: }\boldsymbol{\mathrm{J}} = \begin{bmatrix}
         *     \dfrac{dx}{d\xi} & \dfrac{dx}{d\eta}\\
         *     \dfrac{dy}{d\xi} & \dfrac{dy}{d\eta}
         *    \end{bmatrix}
         *    &\text{3D: }\boldsymbol{\mathrm{J}} = \begin{bmatrix}
         *     \dfrac{dx}{d\xi} & \dfrac{dx}{d\eta} & \dfrac{dx}{d\chi} \\
         *     \dfrac{dy}{d\xi} & \dfrac{dy}{d\eta} & \dfrac{dy}{d\chi} \\
         *     \dfrac{dz}{d\xi} & \dfrac{dz}{d\eta} & \dfrac{dz}{d\chi}
         *    \end{bmatrix}
         * \f}
         * \param i The row index inside the Jacobian matrix.
         * \param j The column index inside the Jacobian matrix.
         * \return The element (i,j) of the Jacobian matrix.
         */
        inline double getJ(unsigned int i, unsigned int j) const noexcept;

        /**
         * \brief Get an element of the inverse of the Jacobian matrix of the change of variable to the reference space:
         * \f{eqnarray*}{
         *    &\text{2D: }\boldsymbol{\mathrm{J}}^{-1}  = \begin{bmatrix}
         *     \dfrac{d\xi}{dx} & \dfrac{d\xi}{dy}\\
         *     \dfrac{d\eta}{dx} & \dfrac{d\eta}{dy}
         *    \end{bmatrix}
         *    &\text{3D: }\boldsymbol{\mathrm{J}}^{-1}  = \begin{bmatrix}
         *     \dfrac{d\xi}{dx} & \dfrac{d\xi}{dy} & \dfrac{d\xi}{dz} \\
         *     \dfrac{d\eta}{dx} & \dfrac{d\eta}{dy} & \dfrac{d\eta}{dz} \\
         *     \dfrac{d\chi}{dx} & \dfrac{d\chi}{dy} & \dfrac{d\chi}{dz}
         *    \end{bmatrix}
         * \f}
         * \param i The row index inside the inverse of the Jacobian matrix.
         * \param j The column index inside the inverse of the Jacobian matrix.
         * \return The element (i,j) of the inverse of the Jacobian matrix.
         */
        inline double getInvJ(unsigned int i, unsigned int j) const noexcept;

        /// \return The highest separation between 2 nodes of the element.
        inline double getLargestExtension() const noexcept;

        /// \get the private value of the circumscribed radius of that element.
        double getCircumscribedRadius() const noexcept;


        /// \update The highest separation between 2 nodes of the element.
        void updateLargestExtension();

        /// \make the connection between the nodes of the elements
        void build(std::vector<std::size_t> nodeIndices, std::vector<Node>& nodesList);
            
        /// \param gp The position of the gauss point in the reference space.
        /// \return The position in real space.
        std::array<double, 3> getPosFromGP(const std::array<double, 3>& gp) const noexcept;

        /// \return The element incircle/insphere.
        double getRin() const noexcept;

        /// \return The area or volume of the element.
        double getSize() const noexcept;

        /// \return The mean mesh size among tho nodes belonging to the element.
        double getLocalMeshSize() const noexcept;

        /// \return The minimal mesh size among tho nodes belonging to the element.
        double getMinMeshSize() const noexcept;

        /// \return The square root of the area or the cubic root of the mesh size of the element.
        double getNaturalMeshSize(bool valueAtNodes = false);

        /// \update the private value of the circumscribed radius of this element.
        void updateCircumscribedRadius();

        /// \return The type of element depending on the boundary node.
        //ELEMENT_TYPE getType() const noexcept;

        /// \return Is the element near a boundary?
        bool isContact() const noexcept;

        friend inline bool operator==(const Element& a, const Element& b) noexcept;

        Element& operator=(const Element& element)   = default;
        Element& operator=(Element&& element)        = default;

    private:
        Mesh* m_pMesh;                                  /**< A pointer to the mesh from which the facet comes from. */
        double m_largest_extension;                     /**< contains the highest separation between 2 nodes of the element*/
        bool m_forcedRefinement;

        double m_r;                                      /**< circumscrived radius*/

        std::size_t m_index;                            /** Index of the element in the element list of the mesh*/
        std::vector<std::size_t> m_nodesIndexes;        /**< Indexes of the nodes in the nodes list which compose this element. */
        std::vector<std::size_t> m_elemsIndexes;        /**< Indexes of the neighbors elements. */
        //std::set<Facet*> m_facets;                      /**< set of pointer to the faces which have this element. */

        double m_detJ;                                  /**< Determinant of the Jacobian matrix of the element. */
        std::array<std::array<double, 3>, 3> m_J;       /**< Jacobian matrix of the element. */
        std::array<std::array<double, 3>, 3> m_invJ;    /**< Inverse Jacobian matrix of the element. */


        /// Compute the Jacobian matrix of the change of variable to the reference space.
        void computeJ();

        /// Compute the determinant of Jacobian matrix of the change of variable to the reference space.
        void computeDetJ();

        /// Compute the inverse of Jacobian matrix of the change of variable to the reference space.
        void computeInvJ();

        friend class Mesh;
        friend class Node;
};

#include "Element.inl"

#endif // ELEMENT_HPP_INCLUDED
