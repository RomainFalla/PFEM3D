#include "Element.hpp"

#include <cassert>
#include "Node.hpp"

void Element::computeJ(const std::vector<Node>& nodesList)
{
    m_J = {{{0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}}};

    if(m_nodesIndexes.size() == 3)
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double x2 = nodesList[m_nodesIndexes[2]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);
        double y2 = nodesList[m_nodesIndexes[2]].getPosition(1);

        m_J[0][0] = x1 - x0;
        m_J[0][1] = x2 - x0;
        m_J[1][0] = y1 - y0;
        m_J[1][1] = y2 - y0;
    }
    else
    {
        double x0 = nodesList[m_nodesIndexes[0]].getPosition(0);
        double x1 = nodesList[m_nodesIndexes[1]].getPosition(0);
        double x2 = nodesList[m_nodesIndexes[2]].getPosition(0);
        double x3 = nodesList[m_nodesIndexes[3]].getPosition(0);
        double y0 = nodesList[m_nodesIndexes[0]].getPosition(1);
        double y1 = nodesList[m_nodesIndexes[1]].getPosition(1);
        double y2 = nodesList[m_nodesIndexes[2]].getPosition(1);
        double y3 = nodesList[m_nodesIndexes[3]].getPosition(1);
        double z0 = nodesList[m_nodesIndexes[0]].getPosition(2);
        double z1 = nodesList[m_nodesIndexes[1]].getPosition(2);
        double z2 = nodesList[m_nodesIndexes[2]].getPosition(2);
        double z3 = nodesList[m_nodesIndexes[3]].getPosition(2);

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
