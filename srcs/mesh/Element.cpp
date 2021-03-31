#include "Element.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Facet.hpp"
#include "Mesh.hpp"

Element::Element(Mesh& mesh):
m_pMesh(&mesh)
{
    m_forcedRefinement = false;
    m_largest_extension = 0;
}

void Element::computeJ()
{
    m_J = {{{0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}}};

    if(m_nodesIndexes.size() == 3)
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);

        m_J[0][0] = x1 - x0;
        m_J[0][1] = x2 - x0;
        m_J[1][0] = y1 - y0;
        m_J[1][1] = y2 - y0;
    }
    else
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);
        const Node& n3 = m_pMesh->getNode(m_nodesIndexes[3]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double x3 = n3.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);
        double y3 = n3.getCoordinate(1);
        double z0 = n0.getCoordinate(2);
        double z1 = n1.getCoordinate(2);
        double z2 = n2.getCoordinate(2);
        double z3 = n3.getCoordinate(2);

        m_J[0][0] = x1 - x0;
        m_J[0][1] = x2 - x0;
        m_J[0][2] = x3 - x0;
        m_J[1][0] = y1 - y0;
        m_J[1][1] = y2 - y0;
        m_J[1][2] = y3 - y0;
        m_J[2][0] = z1 - z0;
        m_J[2][1] = z2 - z0;
        m_J[2][2] = z3 - z0;
    }
}

void Element::computeDetJ()
{
    if(m_nodesIndexes.size() == 3)
    {
        m_detJ = m_J[0][0]*m_J[1][1] - m_J[1][0]*m_J[0][1];
    }
    else
    {
        m_detJ = m_J[0][0]*m_J[1][1]*m_J[2][2]
               + m_J[0][1]*m_J[1][2]*m_J[2][0]
               + m_J[0][2]*m_J[1][0]*m_J[2][1]
               - m_J[2][0]*m_J[1][1]*m_J[0][2]
               - m_J[2][1]*m_J[1][2]*m_J[0][0]
               - m_J[2][2]*m_J[1][0]*m_J[0][1];
    }
}

void Element::computeInvJ()
{
    assert(m_detJ != 0);

    m_invJ = {{{0, 0, 0},
               {0, 0, 0},
               {0, 0, 0}}};

    if(m_nodesIndexes.size() == 3)
    {
        m_invJ[0][0] = m_J[1][1]/m_detJ;

        m_invJ[0][1] = - m_J[0][1]/m_detJ;

        m_invJ[1][0] = - m_J[1][0]/m_detJ;

        m_invJ[1][1] = m_J[0][0]/m_detJ;
    }
    else
    {
        m_invJ[0][0] = (m_J[1][1]*m_J[2][2]
                     - m_J[1][2]*m_J[2][1])/m_detJ;

        m_invJ[0][1] = (m_J[2][1]*m_J[0][2]
                     - m_J[2][2]*m_J[0][1])/m_detJ;

        m_invJ[0][2] = (m_J[0][1]*m_J[1][2]
                     - m_J[0][2]*m_J[1][1])/m_detJ;

        m_invJ[1][0] = (m_J[2][0]*m_J[1][2]
                     - m_J[1][0]*m_J[2][2])/m_detJ;

        m_invJ[1][1] = (m_J[0][0]*m_J[2][2]
                     - m_J[2][0]*m_J[0][2])/m_detJ;

        m_invJ[1][2] = (m_J[1][0]*m_J[0][2]
                     - m_J[0][0]*m_J[1][2])/m_detJ;

        m_invJ[2][0] = (m_J[1][0]*m_J[2][1]
                     - m_J[2][0]*m_J[1][1])/m_detJ;

        m_invJ[2][1] = (m_J[2][0]*m_J[0][1]
                     - m_J[0][0]*m_J[2][1])/m_detJ;

        m_invJ[2][2] = (m_J[0][0]*m_J[1][1]
                     - m_J[1][0]*m_J[0][1])/m_detJ;
    }
}

double Element::getLocalMeshSize() const noexcept
{
    if (m_nodesIndexes.size() == 3)
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);

        return (1. / 3.) * (n0.getLocalMeshSize() + n1.getLocalMeshSize() + n2.getLocalMeshSize());
    }
    else
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);
        const Node& n3 = m_pMesh->getNode(m_nodesIndexes[3]);

        return (1. / 4.) * (n0.getLocalMeshSize() + n1.getLocalMeshSize() + n2.getLocalMeshSize() + n3.getLocalMeshSize());
    }
    return 0;
}

double Element::getMinMeshSize() const noexcept
{
    double result = std::numeric_limits<double>::infinity();
    std::size_t L = m_nodesIndexes.size();
    for (std::size_t i = 0; i < L; i++)
    {
        Node n = m_pMesh->getNode(m_nodesIndexes[i]);
        double temp = n.getLocalMeshSize();
        if (temp < result) result = temp;
    }
    return result;
}

/*ELEMENT_TYPE Element::getType() const noexcept
{
    int counter = 0;
    if (m_nodesIndexes.size() == 3)
    {
        for (int k = 0; k < 3; k++)
        {
            if (m_pMesh->getNode(m_nodesIndexes[k]).isOnBoundary()) counter++;
        }

        if (counter < 2)
            return IN_BULK;
        else 
        {
            if (counter == 2)
                return INSIDE_BOUNDARY;
            else //counter==3
                return OUTSIDE_BOUNDARY;
        }
    }
    else
    {
        for (int k = 0; k < 3; k++)
        {
            if (m_pMesh->getNode(m_nodesIndexes[k]).isOnBoundary()) counter++;
        }

        if (counter < 3)
            return IN_BULK;
        else
        {
            if (counter == 3)
                return INSIDE_BOUNDARY;
            else //counter==4
                return OUTSIDE_BOUNDARY;;
        }
    }
    return IN_BULK; // default value to avoid the compilator complaining
}*/

const Node& Element::getNode(unsigned int nodeIndex) const noexcept
{
    return m_pMesh->getNode(m_nodesIndexes[nodeIndex]);
}

std::array<double, 3> Element::getPosFromGP(const std::array<double, 3>& gp) const noexcept
{
    const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);

    std::array<double, 3> pos;
    pos[0] = m_J[0][0]*gp[0] + m_J[0][1]*gp[1] + m_J[0][2]*gp[2] + n0.getCoordinate(0);
    pos[1] = m_J[1][0]*gp[0] + m_J[1][1]*gp[1] + m_J[1][2]*gp[2] + n0.getCoordinate(1);
    pos[2] = m_J[2][0]*gp[0] + m_J[2][1]*gp[1] + m_J[2][2]*gp[2] + n0.getCoordinate(2);

    return pos;
}

double Element::getSize() const noexcept
{
    return abs(m_detJ*m_pMesh->getRefElementSize(m_pMesh->getDim()));
}

double Element::getNaturalMeshSize(bool valueAtNodes)
{
    if (!valueAtNodes) 
    {
        double fraction = 1. / double(m_pMesh->getDim());
        double result = std::pow(getSize(), fraction);
        return result;
    }
    else
    {
        if (m_nodesIndexes.size() == 3)
        {
            const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
            const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
            const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);

            return (1. / 3.) * (n0.getNaturalMeshSize() + n1.getNaturalMeshSize() + n2.getNaturalMeshSize());
        }
        else
        {
            const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
            const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
            const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);
            const Node& n3 = m_pMesh->getNode(m_nodesIndexes[3]);

            return (1. / 4.) * (n0.getNaturalMeshSize() + n1.getNaturalMeshSize() + n2.getNaturalMeshSize() + n3.getNaturalMeshSize());
        }
        return 0;
    }
}

double Element::getRin() const noexcept
{
    if(m_pMesh->getDim() == 2)
    {
        double A = getSize();
        const Node& n0 = this->getNode(0);
        const Node& n1 = this->getNode(1);
        const Node& n2 = this->getNode(2);

        double a = Node::distance(n0, n1);
        double b = Node::distance(n1, n2);
        double c = Node::distance(n0, n2);
        double s = (a + b + c)/2;

        return A/s;
    }
    else
    {
//        //Plane equation for each tetrahedron face
//        std::array<double, 4> a;
//        std::array<double, 4> b;
//        std::array<double, 4> c;
//        std::array<double, 4> d;
//
//        const Node& n0 = this->getNode(0);
//        const Node& n1 = this->getNode(1);
//        const Node& n2 = this->getNode(2);
//        const Node& n3 = this->getNode(3);
//
//        //Position of the tetrahedron center of gravity
//        std::array<double, 3> centerPos = {
//            (n0.getCoordinate(0) + n1.getCoordinate(0) + n2.getCoordinate(0) + n3.getCoordinate(0))/3,
//            (n0.getCoordinate(1) + n1.getCoordinate(1) + n2.getCoordinate(1) + n3.getCoordinate(1))/3,
//            (n0.getCoordinate(2) + n1.getCoordinate(2) + n2.getCoordinate(2) + n3.getCoordinate(2))/3
//        };
//
//        //Determine for each face the plane equation ax + by + cz + d = 0
//        // Then add a line to U matrix [a_i/S b_i/S c_i/S], S = sqrt(a_i + b_i + c_i)
//        Eigen::Matrix<double, 4, 3> U;
//        Eigen::Vector4d p;
//        for(unsigned int f = 0 ; f < 4 ; ++f)
//        {
//            std::array<double, 3> faceCenter = {0, 0, 0};
//
//            //[x0 y0 z0 ; x1 y1 z1 ; x2 y2 z2]*(a, b, c) = (1, 1, 1)
//            Eigen::Matrix3d A;
//            for(unsigned int n = 0 ; n < 3 ; ++n)
//            {
//                for(unsigned int dd = 0 ; dd < m_pMesh->getDim() ; ++dd)
//                {
//                    A(n, dd) = this->getNode((f + n)%4).getCoordinate(dd);
//                    faceCenter[dd] += this->getNode((f + n)%4).getCoordinate(dd);
//                }
//            }
//            Eigen::Vector3d bb = {1, 1, 1};
//
//            for(unsigned int dd = 0 ; dd < m_pMesh->getDim() ; ++dd)
//                faceCenter[dd] /= 3;
//
//            Eigen::Vector3d planeCoeff = A.colPivHouseholderQr().solve(bb);
//            a[f] = planeCoeff[0];
//            b[f] = planeCoeff[1];
//            c[f] = planeCoeff[2];
//            d[f] = 1;
//
//            std::array<double, 3> cToCVec = {
//                centerPos[0] - faceCenter[0],
//                centerPos[1] - faceCenter[1],
//                centerPos[2] - faceCenter[2]
//            };
//
//            //In order to get rid of abs value in the distance formula and create a
//            //system of linear equations, we ensure the plane normal is oriented towards
//            //the tetrahedron center
//            if(a[f]*cToCVec[0] + b[f]*cToCVec[1] + c[f]*cToCVec[2] < 0)
//            {
//                a[f] *= -1;
//                b[f] *= -1;
//                c[f] *= -1;
//                d[f] *= -1;
//            }
//
//            U(f, 0) = a[f];
//            U(f, 1) = b[f];
//            U(f, 2) = c[f];
//            p(f) = d[f];
//        }
//
//        //Least square approach
//        Eigen::Matrix3d UTU = U.transpose()*U;
//        Eigen::Vector3d UTp = U.transpose()*p;
//
//        Eigen::Vector3d xC = UTU.colPivHouseholderQr().solve(UTp);
//
//        return std::abs(a[0]*xC(0) + b[0]*xC(1) + c[0]*xC(2) + d[0])/std::sqrt(a[0]*a[0] + b[0]*b[0] + c[0]*c[0]);

        return m_pMesh->getHchar();
    }
}

std::vector<double> Element::getState(unsigned int stateIndex) const noexcept
{
    std::vector<double> states(m_nodesIndexes.size());

    for(std::size_t i = 0 ; i < states.size() ; ++i)
        states[i] = m_pMesh->getNode(m_nodesIndexes[i]).getState(stateIndex);

    return states;
}

bool Element::isContact() const noexcept
{
    for(std::size_t i = 0 ; i < m_nodesIndexes.size() ; ++i)
    {
        const Node& node = m_pMesh->getNode(m_nodesIndexes[i]);
        if(node.isBound())
            return true;
    }

    return false;
}

//could use a default value, but use another constructor avoid to check  the if condition each time an element is constructed
void Element::updateLargestExtension()
{
    m_largest_extension = 0.;
    double dist = 0.;
    std::size_t L = m_nodesIndexes.size();
    for (std::size_t i = 0; i < L; i++)
    {
        Node n1 = m_pMesh->getNode(m_nodesIndexes[i]);
        for (std::size_t j = 0; j < L; j++)
        {
            Node n2 = m_pMesh->getNode(m_nodesIndexes[j]);
            dist = Node::distance(n1, n2);
            if (dist > m_largest_extension)
                m_largest_extension = dist;
        }
    }
}
void Element::updateCircumscribedRadius()
{
    if (m_nodesIndexes.size() == 3) 
    {
        Eigen::MatrixXd A(3, 3);
        Eigen::MatrixXd Bx(3,3);
        Eigen::MatrixXd By(3,3);
        Eigen::MatrixXd C(3, 3);

        for (std::size_t i = 0; i < 3; i++) 
        {
            Node n = m_pMesh->getNode(m_nodesIndexes[i]);

            double x = n.getCoordinate(0);
            double y = n.getCoordinate(1);
            double R2 = std::pow(x, 2.) + std::pow(y, 2.);

            A(i, 0) = x; A(i, 1) = y; A(i, 2) = 1;

            Bx(i, 0) = R2;  Bx(i, 1) = y; Bx(i, 2) = 1.;
            By(i, 0) = R2; By(i, 1) = x; By(i, 2) = 1.;

            C(i, 0) = R2; C(i, 1) = x;  C(i, 2) = y;
        }
        double a = A.determinant();
        double bx = -Bx.determinant();
        double by = By.determinant();
        double c = -C.determinant();

        m_r = std::sqrt(std::pow(bx, 2) + std::pow(by, 2) - 4 * a * c) / (2. * fabs(a));
    }
    else //m_nodesIndexes.size() == 4
    {
        Eigen::MatrixXd A(4, 4);
        Eigen::MatrixXd Dx(4, 4);
        Eigen::MatrixXd Dy(4, 4);
        Eigen::MatrixXd Dz(4, 4);
        Eigen::MatrixXd C(4, 4);

        for (std::size_t i = 0; i < 4; i++)
        {
            Node n = m_pMesh->getNode(m_nodesIndexes[i]);

            double x = n.getCoordinate(0);
            double y = n.getCoordinate(1);
            double z = n.getCoordinate(2);

            double R2 = std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.);

            A(i, 0) = x; A(i, 1) = y; A(i, 2) = z; A(i, 3) = 1;

            Dx(i, 0) = R2;  Dx(i, 1) = y; Dx(i, 2) = z; Dx(i, 3) = 1.;
            Dy(i, 0) = R2; Dy(i, 1) = x; Dy(i, 2) = z; Dy(i, 3) = 1.;
            Dz(i, 0) = R2;  Dz(i, 1) = x; Dz(i, 2) = y; Dz(i, 3) = 1.;

            C(i, 0) = R2; C(i, 1) = x; C(i, 2) = y; C(i, 3) = z;
        }

        double a = A.determinant();
        double dx = Dx.determinant();
        double dy = -Dy.determinant();
        double dz = Dz.determinant();
        double c = C.determinant();

        m_r = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2) - 4 * a * c) / (2 * fabs(a));
    }
}


/// \make the connection between the nodes of the elements
void Element::build(std::vector<std::size_t> nodeIndices, std::vector<Node> &nodesList)
{
    computeJ();
    computeDetJ();
    computeInvJ();

    std::size_t L = nodeIndices.size();
    for (std::size_t  i = 0; i < L; i++)
    {
        for (std::size_t  j = 0; j < L; j++)
        {
            if (j == i) continue;

            nodesList[nodeIndices[i]].m_neighbourNodes.push_back(nodeIndices[j]);
        }
    }

    updateCircumscribedRadius();
}

