#pragma once
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <string>
#include <map>

#include "Node.hpp"
#include "Element.hpp"
#include "Facet.hpp"

#include "mesh_defines.h"

/**
 * \struct MeshCreateInfo
 * \brief Structure containing the required parameters to initialize a Mesh class
 */
struct RefinementInfo
{
    double alphaRatio = 1.8; /**< The required \f$ \alpha \f$  ratio value for the weighted \f$ \alpha \f$-shape algorithm*/
    double gamma2 = 0.7; /**< A node should be deleted if ...*/
    double minTargetMeshSize = 0.; /**< The minimal target mesh sizeinside the mesh*/
    double maxTargetMeshSize = 0.; /**< The maximal target mesh sizeinside the mesh*/
    double maxProgressionFactor = 1.25; /**< Maximal progression of the target mesh size from node to node*/
    std::vector<int> imposedTargetMeshSizeBCtags = {};
    std::vector<double> imposedTargetMeshSizeValues = {};
};

struct MeshCreateInfo
{
    double hchar = 0; /**< The characteristic element size of the mesh*/
    double alpha = 1; /**< The required \f$ \alpha \f$ value for the \f$ \alpha \f$-shape algorithm*/
    double gamma = 0; /**< A node should be deleted if \f$ d(n_1, n_2) < \gamma h_{char} \f$ */
    double omega = 1e16; /**< A node should be added in the center of an element if \f$ A_{elm} > \omega h_{char}^{dim} \f$ */
    bool useMeshRefinement = false; /**< Define wether we use mesh refinement or not*/
    RefinementInfo refInfo = {};

    std::vector<double> boundingBox = {}; /**< Nodes and elements outised bounding box are deleted. Format:
                                               \f$ [x_{min}, y_{min}, (z_{min}, )x_{max}, y_{max}, (z_{max}) ] \f$ */
    std::string mshFile = {}; /**< The path to the .msh file to load */
};

enum ELEMENT_TYPE
{
    IN_BULK = 0,
    INSIDE_BOUNDARY = 1,
    OUTSIDE_BOUNDARY = 2,
};
/**
 * \class Mesh
 * \brief Represents a Lagrangian mesh.
 *
 * A mesh in the PFEM is represented as a collection of nodes, which are linked together through triangular/tetrahedral elements.
 * Based on an initial cloud of nodes, a Delaunay triangulation is performed, followed by an alpha-shape algorithm step, which discard
 * too big or to distorted elements, following:
 * \f[
 *    r_{\text{circumcircle/sphere}}^2 > \alpha \rightarrow \text{element discarded}
 * \f]
 * The class allows notably to retrive the nodes, their states, the elements, the boundary facets, their normals and curvatures.
 */
class MESH_API Mesh
{
    public:
        Mesh()                              = delete;
        /// \param meshInfos a reference to a MeshCreateInfo structure
        Mesh(const MeshCreateInfo& meshInfos);
        //Mesh(const MeshCreateInfo& meshInfos, const RefinementInfo& refInfos);
        Mesh(const Mesh& mesh)              = delete;
        Mesh& operator=(const Mesh& mesh)   = delete;
        Mesh(Mesh&& mesh)                   = delete;
        Mesh& operator=(Mesh&& mesh)        = delete;
        ~Mesh()                             = default;

        /// \brief Display the mesh parameters to console.
        void displayToConsole() const noexcept;

        /// \param nodeIndex The index of the boundary node in the nodes list.
        /// \return The initial position of the node.
        inline std::array<double, 3> getBoundNodeInitPos(std::size_t nodeIndex) const;

        /// \return The mesh dimension.
        inline unsigned short getDim() const noexcept;

        /// \param elm The index of the element.
        /// \return A reference to the element.
        inline const Element& getElement(std::size_t elm) const noexcept;

        /// \return The number of elements in the mesh.
        inline std::size_t getElementsCount() const noexcept;

        /// \param facet The index of the face.
        /// \return A reference to the face.
        inline const Facet& getFacet(std::size_t facet) const noexcept;

        /// \return The number of boundary face in the mesh.
        inline std::size_t getFacetsCount() const noexcept;

        /// \param facetIndex The index of the facet in the facets list.
        /// \return The physical group of that facet.
        inline std::string getFacetType(std::size_t facetIndex) const noexcept;

        /// \param nodeIndex The index of the node in the nodes list.
        /// \return The curvature of the boundary at the node.
        inline double getFreeSurfaceCurvature(std::size_t nodeIndex) const;

        /// \return true in the case we use mesh refinement
        inline bool useMeshRefinement() const noexcept;

        /// \param nodeIndex The index of the boudary or free surface node in the nodes list.
        /// \return The exterior normal of the boundary at the node.
        inline std::array<double, 3> getBoundFSNormal(std::size_t nodeIndex) const;

        /// \param dimension The dimension of the reference element on which you want the weights.
        /// \param n The requested number of Gauss points.
        /// \return The Gauss points (x, y, z).
        std::vector<std::array<double, 3>> getGaussPoints(unsigned int dimension, unsigned int n) const;

        /// \param dimension The dimension of the reference element on which you want the weights.
        /// \param n The requested number of Gauss points.
        /// \return The weight associated to those Gauss points.
        std::vector<double> getGaussWeight(unsigned int dimension, unsigned int n) const;

        /// \return The characteristic size of the mesh.
        inline double getHchar() const noexcept;

        /// \return a constant reference of a node in the mesh
        inline const Node& getNode(std::size_t nodeIndex) const noexcept;
        
        /// \return a non-constant pointer to a node in the mesh
        inline Node* getNodePtr(std::size_t nodeIndex);

        /// \return The number of nodes in the mesh.
        inline std::size_t getNodesCount() const noexcept;

        /// \return The number of nodes in an element.
        inline unsigned short getNodesPerElm() const noexcept;

        /// \return The number of nodes in a facet.
        inline unsigned short getNodesPerFacet() const noexcept;

        /// \param nodeIndex The index of the node in the nodes list.
        /// \return The physical group of that node.
        inline std::string getNodeType(std::size_t nodeIndex) const noexcept;


        /// \param dimension The dimension of the reference element on which you want the weights.
        /// \return The size of the element in the reference coordinate system.
        double getRefElementSize(unsigned int dimension) const;

        /// \update the maximal and minimal target mesh size
        void updateBoundingElementSize();

        /// \update the target target mesh size at each node
        void updateLocalElementSize();

        ///  \make the ling between the facets and the elements
        void linkElementsAndFacets();

        /**
         * \param dimension The dimension of the reference element on which you want the shape functions.
         * \param n The requested number of Gauss points.
         * \return Get the shape functions evaluated at each gauss points in the reference space in the format:
         * [[sf1_gp1, ..., sfn_gp1], [sf1_gp2, ..., sfn_gp2], ...].
         */
        std::vector<std::vector<double>> getShapeFunctions(unsigned int dimension, unsigned int n) const;

        /**
         * \param dimension The dimension of the reference element on which you want the gradient of shape functions.
         * \return Get the gradient of the shape functions in the reference space in the format:
         * [dsf1/dx ... dsfn/dx ; dsf1/dy ... dsfn/dy ; dsf1/dz ... dsfn/dz].
         */
        std::vector<std::vector<double>> getGradShapeFunctions(unsigned int dimension) const;

        /// \return If normals and curvature are actually computed.
        inline bool isNormalCurvComputed() const noexcept;

        /// \brief Perform remeshing on the mesh.
        void remesh(bool verboseOutput);

        /// \brief Restore the current nodes list from the list saved in saveNodesList.
        void restoreNodesList();

        /// \brief Save the current nodes list in another variable.
        void saveNodesList();

        /**
         * \brief Activate or not the computation of normals and curvature (default is true).
         * \param activate should the computation be activated
         */
        inline void setComputeNormalCurvature(bool activate) noexcept;

        /**
         * \brief Set the requested flag to true.
         * \param nodeIndex The index of the node which will be tagged.
         * \param flag The flag (0-7) that will be set.
         */
        inline void setNodeFlag(std::size_t nodeIndex, unsigned short flag) noexcept;

        /**
         * \brief Set if the node is in a Dirichlet BC (i.e. a BC which we impose speed but do not move!).
         * \param nodeIndex The index of the node in the internal nodes list;
         * \param isFixed true if the node is in a Dirichlet BC, false otherwise;
         */
        inline void setNodeIsFixed(std::size_t nodeIndex, bool isFixed) noexcept;

        /**
         * \brief Set the position of a node.
         * \param nodeIndex The index of the node in the nodes list.
         * \param stateIndex The index of the state.
         * \param state The new value of the state.
         */
        inline void setNodeState(std::size_t nodeIndex, unsigned int stateIndex, double state) noexcept;

        /**
         * \brief Set the number of states to be stored at node level.
         * \param statesNumber The number of state per nodes.
         */
        inline void setStatesNumber(unsigned int statesNumber);

        /**
         * \brief Update the nodes position.
         * \param deltaPos The variation of the coordinate.
         */
        void updateNodesPosition(std::vector<double> deltaPos);

        /**
         * \brief Update the nodes position (from the saved nodeS List).
         * \param deltaPos The variation of the coordinate.
         */
        void updateNodesPositionFromSave(std::vector<double> deltaPos);

    private:
        double m_hchar; /**< Characteristic size of an element (same as in .geo file). */
        double m_alpha; /**< Alpha parameter of the alpha-shape algorithm (triangles are discared if  r_circumcircle > alpha*hchar). */
        double m_omega; /**< Control the addition of node if a triangle is too big (a node is added if A_triangle > omege*hchar^2). */
        double m_gamma; /**< Control the deletetion of node if two are too close to each other (a node is deleted if d_nodes < gamma*hchar). */
        std::vector<double> m_boundingBox; /**< Box delimiting the zone of nodes existence (format: [xmin, ymin, xmax, ymax]). */
        bool m_useMeshRefinement; /**< Define wether we use mesh refinement or not*/
        

        bool m_computeNormalCurvature;      /**< Control if mesh update should compute normals and curvatures of free surface. */

        double m_alphaRatio; /**< The required \f$ \alpha \f$  ratio value for the weighted \f$ \alpha \f$-shape algorithm*/
        double m_gamma2; /**< A node should be deleted if ...*/
        double m_minTargetMeshSize; /**< The minimal target mesh size inside the mesh*/
        double m_maxTargetMeshSize; /**< The maximal target mesh size inside the mesh*/
        double m_maxProgressionFactor;/**< Maximal progression of the target mesh size from node to node*/

        double m_refinementLevels[6]; /**Area thresholds levels used for mesh refinement*/
        double m_refinementFactor; /**Area thresholds factor used for mesh refinement*/


        unsigned short m_dim;               /**< The mesh dimension. */

        std::vector<Node> m_nodesList;      /**< List of nodes of the mesh. */
        std::vector<Node> m_nodesListSave;  /**< A copy of the nodes list (usefull for non-linear algorithm). */

        std::vector<int> m_ImposedTargetMeshSizeBCtags; /**vector of vectors to pointers to which the target mesh size is imposed.
                                                                          The number of entries correspond to the number of boundaries */
        std::vector<double> m_ImposedTargetMeshSizeValues;                      /**Target mesh size imposed at some boundaries. 
                                                                          The number of entries correspond to the number of boundaries*/

        std::vector<Element> m_elementsList;    /**< The list of elements. */
        std::vector<Facet> m_facetsList;        /**< The list of boundary facets. */
        std::vector<Facet> m_inBulkFacetList;    /**< The list of non-boundary facets. */

        std::vector<std::string> m_tagNames; /**< The name of the tag of the nodes. */
        std::map<std::size_t, std::array<double, 3>> m_boundaryInitialPos;  /**< Initial position of node at the boundary. */
        std::map<std::size_t, std::array<double, 3>> m_boundFSNormal;   /**< Free surface and boundary normals normals */
        std::map<std::size_t, double> m_freeSurfaceCurvature;               /**< Free surface curvatures */

        /**
         * \brief Add nodes in element whose area is too big (A_tringle > omega*hchar^2.
         * \return true if at least one node was added, false otherwise).
         */
        bool addNodes(bool verboseOutput);

        /// \divide an element by adding several nodes depending on a level parameters = 1,2,3.
        void divideElement(Element& el, unsigned int level) ;


        /// \return The index of the node of *element that is opposite to *face;
        std::size_t getOutNodeIndex(Facet* face, Element* element);

        /// \find a pointer to the largest facet of one element and put it in face
        // getLargestFacet(Element* element, Facet*& face);

        /**
         * \brief Check if a node is outside the bounding box and deletes it if so.
         * \return true if at least one node was deleted, false otherwise.
         */
        bool checkBoundingBox(bool verboseOutput) noexcept;

        /// \brief Compute the mesh dimension from the .msh file.
        void computeMeshDim();

        /// \brief Compute the normal and curvature of each boundary and free surface.
        void computeFSNormalCurvature();

        /// \brief Compute the normal and curvature of each boundary and free surface (2D).
        void computeFSNormalCurvature2D();

        /// \brief Compute the normal and curvature of each boundary and free surface (3D).
        void computeFSNormalCurvature3D();

        /// \brief function that goes from the nodes with the smallest target mesh size to the neighbours, and limit the progression from node to node to the value "m_maxProgressionFactor".
        void smoothingLcLoc();
        /**
         * \brief Load the nodes from a file using gmsh.
         * \param fileName The name of the .msh file.
         */
        void loadFromFile(const std::string& fileName);

        /// \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation and alpha-shape).
        void triangulateAlphaShape();

        /// \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation and alpha-shape) (2D).
        void triangulateAlphaShape2D();

        /// \brief Remesh the nodes in nodesList using CGAL (Delaunay triangulation and alpha-shape) (3D).
        void triangulateAlphaShape3D();

        /// \brief Remesh the nodes in nodesList using CGAL (Weighted Delaunay triangulation and weighted alpha-shape) (2D).
        void TriangulateWeightedAlphaShape2D();

        /// \brief Remesh the nodes in nodesList using CGAL (Weighted Delaunay triangulation and weighted alpha-shape) (3D).
        void TriangulateWeightedAlphaShape3D();

        double getCircumScribedRadius(std::vector<size_t> nods);

        ELEMENT_TYPE getElementType(std::vector<std::size_t> nodesIndexes);

        /**
         * \brief Removes nodes if they are too close from each other
         *       (d_nodes < gamma*hchar).
         * \return true if at least one node was deleted, false otherwise.
         */
        bool removeNodes(bool verboseOutput) noexcept;
};

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
