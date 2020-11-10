#include "Facet.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "Node.hpp"

void Facet::computeJ(const std::vector<Node>& nodesList)
{
    m_J = {{{0, 0},
            {0, 0},
            {0, 0}}};

    if(m_nodesIndexes.size() == 2)
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);

        m_J[0][0] = (x1 - x0)/2; //dx/dxi
        m_J[1][0] = (y1 - y0)/2; //dy/dxi
    }
    else
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double x2 = nodesList[m_nodesIndexes[2]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);
        double y2 = nodesList[m_nodesIndexes[2]].getPosition(1);
        double z0 = nodesList[m_nodesIndexes[0]].getPosition(2);
        double z1 = nodesList[m_nodesIndexes[1]].getPosition(2);
        double z2 = nodesList[m_nodesIndexes[2]].getPosition(2);

        m_J[0][0] = x1 - x0; //dx/dxi
        m_J[1][0] = y1 - y0; //dy/dxi
        m_J[2][0] = z1 - z0; //dz/dxi
        m_J[0][1] = x2 - x0; //dx/deta
        m_J[1][1] = y2 - y0; //dy/deta
        m_J[2][1] = z2 - z0; //dz/deta
    }
}

void Facet::computeDetJ()
{
    if(m_nodesIndexes.size() == 2)
    {
        m_detJ = std::sqrt(m_J[0][0]*m_J[0][0] + m_J[1][0]*m_J[1][0]); // || (dx/dxi; dy/dxi)||
    }
    else
    {
        std::array<double, 3> dGammaVec = {
            m_J[1][0]*m_J[2][1] - m_J[1][1]*m_J[2][0],
            m_J[2][0]*m_J[0][1] - m_J[2][1]*m_J[0][0],
            m_J[0][0]*m_J[1][1] - m_J[1][0]*m_J[0][1]
        };

        // || (dx/dxi; dy/dxi; dz/dxi) x (dx/deta; dy/deta; dz/deta)||
        m_detJ = std::sqrt(dGammaVec[0]*dGammaVec[0] + dGammaVec[1]*dGammaVec[1] + dGammaVec[2]*dGammaVec[2]);
    }
}

void Facet::computeInvJ(const std::vector<Node>& nodesList)
{
    m_invJ = {{{0, 0, 0},
               {0, 0, 0}}};

    if(m_nodesIndexes.size() == 2)
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);

        m_invJ[0][0] = 1/m_J[0][0]; //dxi/dx
        m_invJ[0][1] = m_invJ[0][0]*((x1-x0)/(y1-y0)); //dxi/dy = dxi/dx*dx/dy
    }
    else
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double x2 = nodesList[m_nodesIndexes[2]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);
        double y2 = nodesList[m_nodesIndexes[2]].getPosition(1);
        double z0 = nodesList[m_nodesIndexes[0]].getPosition(2);
        double z1 = nodesList[m_nodesIndexes[1]].getPosition(2);
        double z2 = nodesList[m_nodesIndexes[2]].getPosition(2);

        double detB = m_J[0][0]*m_J[1][1] - m_J[0][1]*m_J[1][0];

        m_invJ[0][0] = m_J[1][1]/detB;  //dxi/dx
        m_invJ[0][1] = -m_J[0][1]/detB; //dxi/dy
        m_invJ[1][0] = -m_J[1][0]/detB; //deta/dx
        m_invJ[1][1] = m_J[0][0]/detB;  //deta/dy

        double den = x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1);

        double a_1 = den/(z0*(y1 - y2) + z1*(y1 - y0) + z2*(y0 - y1));
        double b_1 = den/(z0*(x2 - x1) + z1*(x0 - x2) + z2*(x1 - x0));

        m_invJ[0][2] = m_invJ[0][0]*a_1 + m_J[0][1]*b_1; //dxi/dz = dxi/dx*dx/dz + dxi/dy*dy/dz
        m_invJ[1][2] = m_invJ[1][0]*a_1 + m_J[1][1]*b_1; //deta/dz = deta/dx*dx/dz + deta/dy*dy/dz
    }
}
