#pragma once
#ifndef FACET_HPP_INCLUDED
#define FACET_HPP_INCLUDED

#include <array>
#include <cstdint>
#include <vector>
#include <set>

#include "mesh_defines.h"

class Mesh;
class Element;
class Node;

/**
 * \class Facet
 * \brief Represents a facet in the PFEM.
 */
class MESH_API Facet
{
    public:
        /// \param mesh a reference to the mesh
        Facet(Mesh& mesh);
        Facet(Mesh& mesh, std::vector<std::size_t> nodeIndexes);
        Facet(const Facet& facet)      = default;
        Facet(Facet&& facet)           = default;
        ~Facet()                       = default;

        /// \param nodeFacetIndex The index of the node in the facet.
        /// \return The index of the facet's node in the mesh.
        inline std::size_t getNodeIndex(unsigned int nodeFacetIndex) const noexcept;

        /// \return The index of the facet's opposite node in the mesh.
        inline std::size_t getOutNodeIndex() const noexcept;

        /// \param nodeFacetIndex The index of the node in the facet.
        /// \return A constant reference to the facet's node.
        const Node& getNode(unsigned int nodeFacetIndex) const noexcept;

        /// \return A constant reference to the facet's opposite node.
        const Node& getOutNode() const noexcept;


        ///  \return true if this facet is a facet with all its node having "isBound=true" ;
        bool isBound() const noexcept;

        ///  \return true if this facet is a facet with all its node of a boundary (whatever the type of it) ;
        bool isOnBoundary() const noexcept;

        /// \return The length or area of the facet.
        double getSize() const noexcept;

        /// \return The local target mesh size of the facet.
        double getLocalMeshSize();

        /// \param stateIndex The index of the state.
        /// \return A vector containing the considered state for each node of the facet.
        std::vector<double> getState(unsigned int stateIndex) const noexcept;

        /// \return The determinant of the Jacobian matrix of the change of variable to the reference space.
        inline double getDetJ() const noexcept;

        /**
         * \brief Get an element of the Jacobian matrix of the change of variable to the reference space:
         * \f{eqnarray*}{
         *    &\text{2D: }\boldsymbol{\mathrm{J}} = \begin{bmatrix}
         *     \dfrac{dx}{d\xi} \\
         *     \dfrac{dy}{d\xi}
         *    \end{bmatrix}
         *    &\text{3D: }\boldsymbol{\mathrm{J}} = \begin{bmatrix}
         *     \dfrac{dx}{d\xi} & \dfrac{dx}{d\eta} \\
         *     \dfrac{dy}{d\xi} & \dfrac{dy}{d\eta} \\
         *     \dfrac{dz}{d\xi} & \dfrac{dz}{d\eta}
         *    \end{bmatrix}
         * \f}
         * \param i The row index inside the Jacobian matrix.
         * \param j The column index inside the Jacobian matrix.
         * \return The element (i,j) of the Jacobian matrix.
         */
        inline double getJ(unsigned int i, unsigned int j) const noexcept;

        /**
         * \brief Get an element of the inverse of Jacobian matrix of the change of variable to the reference space:
         * \f{eqnarray*}{
         *    &\text{2D: }\boldsymbol{\mathrm{J}}^{-1} = \begin{bmatrix}
         *     \dfrac{d\xi}{dx} & \dfrac{d\xi}{dy}
         *    \end{bmatrix}
         *    &\text{3D: }\boldsymbol{\mathrm{J}}^{-1}  = \begin{bmatrix}
         *     \dfrac{d\xi}{dx} & \dfrac{d\xi}{dy} & \dfrac{d\xi}{dz} \\
         *     \dfrac{d\eta}{dx} & \dfrac{d\eta}{dy} & \dfrac{d\eta}{dz}
         *    \end{bmatrix}
         * \f}
         * \param i The row index inside the inverse of Jacobian matrix.
         * \param j The column index inside the inverse of Jacobian matrix.
         * \return The element (i,j) of the inverse of Jacobian matrix.
         */
        inline double getInvJ(unsigned int i, unsigned int j) const noexcept;

        /// \param gp The position of the gauss point in the reference space.
        /// \return The position in real space.
        std::array<double, 3> getPosFromGP(const std::array<double, 3>& gp) const noexcept;

        friend inline bool operator==(const Facet& a, const Facet& b) noexcept;

        Facet& operator=(const Facet& facet)   = default;
        Facet& operator=(Facet&& facet)        = default;

    private:
        Mesh* m_pMesh;                                  /**< A pointer to the mesh from which the facet comes from. */

        std::vector<std::size_t> m_nodesIndexes;        /**< Indexes of the nodes in the nodes list which compose this facet. */
        //std::size_t m_outNodeIndex;                     /**< Indexes of the node which is "in front of" this facet. */
        //std::size_t m_elementIndex;                     /**< Index of the boundary element. */
        std::vector<std::size_t> m_outNodeIndexes;        /**< Indexes of the node which is "in front of" this facet. */
        std::vector<std::size_t> m_elementIndexes;        /**< Index of the elements that share this facet. */

        double m_detJ;                                  /**< Determinant of the Jacobian matrix of the facet. */
        std::array<std::array<double, 2>, 3> m_J;       /**< Jacobian matrix of the facet. */
        std::array<std::array<double, 3>, 2> m_invJ;    /**< Inverse Jacobian matrix of the facet. */

        /// Compute the Jacobian matrix of the change of variable to the reference space.
        void computeJ();

        /// Compute the determinant of Jacobian matrix of the change of variable to the reference space.
        void computeDetJ();

        /// Compute the inverse of Jacobian matrix of the change of variable to the reference space.
        void computeInvJ();

        friend class Mesh;
};

#include "Facet.inl"

#endif // FACET_HPP_INCLUDED
