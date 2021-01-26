#include "MatricesBuilder.hpp"

#include "../../mesh/Mesh.hpp"

MatrixBuilder::MatrixBuilder(const Mesh& mesh, unsigned int nGPHD, unsigned int nGPLD):
m_mesh(mesh),
m_nGPHD(nGPHD),
m_nGPLD(nGPLD)
{
    m_sfsHD = m_mesh.getShapeFunctions(m_mesh.getDim(), m_nGPHD);
    m_gradsfsHD = m_mesh.getGradShapeFunctions(m_mesh.getDim());
    m_sfsLD = m_mesh.getShapeFunctions(m_mesh.getDim() - 1, m_nGPLD);
    m_gradsfsLD = m_mesh.getGradShapeFunctions(m_mesh.getDim() - 1);

    m_NHD.resize(m_nGPHD);
    m_NLD.resize(m_nGPLD);
    m_NHDtilde.resize(m_nGPHD);
    m_NLDtilde.resize(m_nGPLD);

    unsigned int counter = 0;
    for(std::vector<double> sf : m_sfsHD)
    {
        switch(m_mesh.getDim())
        {
            case 1:
            {
                m_NHD[counter].resize(1, 2);
                m_NHDtilde[counter].resize(1, 2); m_NHDtilde[counter].setZero();
                m_NHDtilde[counter](0,0) = m_NHD[counter](0, 0) = sf[0];
                m_NHDtilde[counter](0,1) = m_NHD[counter](0, 1) = sf[1];
                break;
            }

            case 2:
            {
                m_NHD[counter].resize(1, 3);
                m_NHDtilde[counter].resize(2, 6); m_NHDtilde[counter].setZero();
                m_NHDtilde[counter](0,0) = m_NHDtilde[counter](1,3) =  m_NHD[counter](0, 0) = sf[0];
                m_NHDtilde[counter](0,1) = m_NHDtilde[counter](1,4) =  m_NHD[counter](0, 1) = sf[1];
                m_NHDtilde[counter](0,2) = m_NHDtilde[counter](1,5) =  m_NHD[counter](0, 2) = sf[2];
                break;
            }

            default:
            {
                m_NHD[counter].resize(1, 4);
                m_NHDtilde[counter].resize(3, 12); m_NHDtilde[counter].setZero();
                m_NHDtilde[counter](0,0) = m_NHDtilde[counter](1,4) = m_NHDtilde[counter](2,8)  = m_NHD[counter](0, 0) = sf[0];
                m_NHDtilde[counter](0,1) = m_NHDtilde[counter](1,5) = m_NHDtilde[counter](2,9)  = m_NHD[counter](0, 1) = sf[1];
                m_NHDtilde[counter](0,2) = m_NHDtilde[counter](1,6) = m_NHDtilde[counter](2,10) = m_NHD[counter](0, 2) = sf[2];
                m_NHDtilde[counter](0,3) = m_NHDtilde[counter](1,7) = m_NHDtilde[counter](2,11) = m_NHD[counter](0, 3) = sf[3];
                break;
            }
        }
        counter++;
    }

    counter = 0;
    for(std::vector<double> sf : m_sfsLD)
    {
        switch(m_mesh.getDim())
        {
            case 1:
                break;

            case 2:
            {
                m_NLD[counter].resize(1, 2);
                m_NLDtilde[counter].resize(2, 4); m_NLDtilde[counter].setZero();
                m_NLDtilde[counter](0, 0) = m_NLDtilde[counter](1, 2) = m_NLD[counter](0, 0) = sf[0];
                m_NLDtilde[counter](0, 1) = m_NLDtilde[counter](1, 3) = m_NLD[counter](0, 1) = sf[1];
                break;
            }

            case 3:
            {
                m_NLD[counter].resize(1, 3);
                m_NLDtilde[counter].resize(3, 9); m_NLDtilde[counter].setZero();
                m_NLDtilde[counter](0, 0) = m_NLDtilde[counter](1, 3) = m_NLDtilde[counter](2, 6) = m_NLD[counter](0, 0) = sf[0];
                m_NLDtilde[counter](0, 1) = m_NLDtilde[counter](1, 4) = m_NLDtilde[counter](2, 7) = m_NLD[counter](0, 1) = sf[1];
                m_NLDtilde[counter](0, 2) = m_NLDtilde[counter](1, 5) = m_NLDtilde[counter](2, 8) = m_NLD[counter](0, 2) = sf[2];
                break;
            }
        }
        counter++;
    }

    m_gaussWeightHD = m_mesh.getGaussWeight(m_mesh.getDim(), m_nGPHD);
    m_gaussPointsHD = m_mesh.getGaussPoints(m_mesh.getDim(), m_nGPHD);
    m_gaussPointsLD = m_mesh.getGaussPoints(m_mesh.getDim() - 1, m_nGPLD);
    m_gaussWeightLD = m_mesh.getGaussWeight(m_mesh.getDim() - 1, m_nGPLD);
}

MatrixBuilder::~MatrixBuilder()
{

}

Eigen::MatrixXd MatrixBuilder::getGradN(const Element& element)
{
    Eigen::MatrixXd gradN;
    //I could use m_gradsfHD though...
    switch(m_mesh.getDim())
    {
        case 1:
            gradN.resize(1, 2); gradN.setZero();
            gradN(0,0) = - element.getInvJ(0, 0);
            gradN(0,1) = element.getInvJ(0, 0);
            break;

        case 2:
            gradN.resize(2, 3); gradN.setZero();
            gradN(0,0) = - element.getInvJ(0, 0) - element.getInvJ(1, 0);
            gradN(0,1) = element.getInvJ(0, 0);
            gradN(0,2) = element.getInvJ(1, 0);

            gradN(1,0) = - element.getInvJ(0, 1) - element.getInvJ(1, 1);
            gradN(1,1) = element.getInvJ(0, 1);
            gradN(1,2) = element.getInvJ(1, 1);
            break;

        default:
            gradN.resize(3, 4); gradN.setZero();
            gradN(0,0) = - element.getInvJ(0, 0) - element.getInvJ(1, 0) - element.getInvJ(2, 0);
            gradN(0,1) = element.getInvJ(0, 0);
            gradN(0,2) = element.getInvJ(1, 0);
            gradN(0,3) = element.getInvJ(2, 0);

            gradN(1,0) = - element.getInvJ(0, 1) - element.getInvJ(1, 1) - element.getInvJ(2, 1);
            gradN(1,1) = element.getInvJ(0, 1);
            gradN(1,2) = element.getInvJ(1, 1);
            gradN(1,3) = element.getInvJ(2, 1);

            gradN(2,0) = - element.getInvJ(0, 2) - element.getInvJ(1, 2) - element.getInvJ(2, 2);
            gradN(2,1) = element.getInvJ(0, 2);
            gradN(2,2) = element.getInvJ(1, 2);
            gradN(2,3) = element.getInvJ(2, 2);
            break;
    }

    return gradN;
}

Eigen::MatrixXd MatrixBuilder::getB(Eigen::MatrixXd gradN)
{
    Eigen::MatrixXd B;
    //With linear shape functions, Be is constant for all gauss points ^^.
    switch(m_mesh.getDim())
    {
        case 1:
            B.resize(1, 2); B.setZero();
            B(0,0) = gradN(0,0);
            B(0,1) = gradN(0,1);
            break;

        case 2:
            B.resize(3, 6); B.setZero();
            B(0,0) = B(2,3) = gradN(0,0);
            B(0,1) = B(2,4) = gradN(0,1);
            B(0,2) = B(2,5) = gradN(0,2);

            B(1,3) = B(2,0) = gradN(1,0);
            B(1,4) = B(2,1) = gradN(1,1);
            B(1,5) = B(2,2) = gradN(1,2);
            break;

        default:
            B.resize(6, 12); B.setZero();
            B(0,0) = B(3,4) = B(4,8) = gradN(0,0);
            B(0,1) = B(3,5) = B(4,9) = gradN(0,1);
            B(0,2) = B(3,6) = B(4,10) = gradN(0,2);
            B(0,3) = B(3,7) = B(4,11) = gradN(0,3);

            B(1,4) = B(3,0) = B(5,8) = gradN(1,0);
            B(1,5) = B(3,1) = B(5,9) = gradN(1,1);
            B(1,6) = B(3,2) = B(5,10) = gradN(1,2);
            B(1,7) = B(3,3) = B(5,11) = gradN(1,3);

            B(2,8) = B(4,0) = B(5,4) = gradN(2,0);
            B(2,9) = B(4,1) = B(5,5) = gradN(2,1);
            B(2,10) = B(4,2) = B(5,6) = gradN(2,2);
            B(2,11) = B(4,3) = B(5,7) = gradN(2,3);
            break;
    }

    return B;
}

Eigen::MatrixXd MatrixBuilder::getM(const Element& element)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();

    Eigen::MatrixXd M(noPerEl, noPerEl); M.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        M += m_Mfunc(element, m_NHD[i])*m_NHD[i].transpose()*m_NHD[i]*m_gaussWeightHD[i];
    }

    M *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return M;
}

Eigen::MatrixXd MatrixBuilder::getMGamma(const Facet& facet)
{
    unsigned int noPerFacet = m_mesh.getNodesPerFacet();

    Eigen::MatrixXd MGamma(noPerFacet, noPerFacet); MGamma.setZero();

    for(unsigned int i = 0 ; i < m_NLD.size() ; ++ i)
    {
        MGamma += m_MGammafunc(facet, m_NLD[i])*m_NLD[i].transpose()*m_NLD[i]*m_gaussWeightLD[i];
    }

    MGamma *= facet.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim() - 1);

    return MGamma;
}

Eigen::VectorXd MatrixBuilder::getQN(const Facet& facet)
{
    unsigned int noPerFacet = m_mesh.getNodesPerFacet();

    Eigen::VectorXd qn(noPerFacet); qn.setZero();
    Eigen::VectorXd n(noPerFacet*m_mesh.getDim());

    for(std::size_t i = 0 ; i < noPerFacet ; ++i)
    {
        std::array<double, 3> nodeNormal = m_mesh.getBoundFSNormal(facet.getNodeIndex(i));

        for(std::size_t d = 0 ; d < m_mesh.getDim() ; ++d)
            n(i*noPerFacet + d) = nodeNormal[d];
    }

    for(unsigned int i = 0 ; i < m_NLD.size() ; ++ i)
    {
        Eigen::VectorXd q = m_QFunc(facet, m_gaussPointsLD[i]);
        double fact = (q.transpose()*m_NLDtilde[i]*n).value();
        qn += fact*m_NLD[i].transpose();
    }

    return qn;
}

Eigen::MatrixXd MatrixBuilder::getK(const Element& element, const Eigen::MatrixXd& B)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();
    unsigned int dim = m_mesh.getDim();

    Eigen::MatrixXd K(noPerEl*dim, noPerEl*dim); K.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        K += m_Kfunc(element, m_NHD[i], B)*B.transpose()*m_ddev*B*m_gaussWeightHD[i];
    }

    K *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return K;
}

Eigen::MatrixXd MatrixBuilder::getD(const Element& element, const Eigen::MatrixXd& B)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();
    unsigned int dim = m_mesh.getDim();

    Eigen::MatrixXd D(noPerEl, noPerEl*dim); D.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        D += m_Dfunc(element, m_NHD[i], B)*m_NHD[i].transpose()*m_m.transpose()*B*m_gaussWeightHD[i];
    }

    D *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return D;
}

Eigen::MatrixXd MatrixBuilder::getL(const Element& element, const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();

    Eigen::MatrixXd L(noPerEl, noPerEl); L.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        L += m_Lfunc(element, m_NHD[i], B)*gradN.transpose()*gradN*m_gaussWeightHD[i];
    }

    L *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return L;
}

Eigen::MatrixXd MatrixBuilder::getC(const Element& element, const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();
    unsigned int dim = m_mesh.getDim();

    Eigen::MatrixXd C(noPerEl, noPerEl*dim); C.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        C += m_Cfunc(element, m_NHD[i], B)*gradN.transpose()*m_NHDtilde[i]*m_gaussWeightHD[i];
    }

    C *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return C;
}

Eigen::VectorXd MatrixBuilder::getF(const Element& element, const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& B)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();
    unsigned int dim = m_mesh.getDim();

    Eigen::VectorXd F(noPerEl*dim); F.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        F += m_Ffunc(element, m_NHD[i], B)*m_NHDtilde[i].transpose()*vec*m_gaussWeightHD[i];
    }

    F *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return F;
}

Eigen::VectorXd MatrixBuilder::getH(const Element& element, const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN)
{
    unsigned int noPerEl = m_mesh.getNodesPerElm();

    Eigen::VectorXd H(noPerEl); H.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        H += m_Hfunc(element, m_NHD[i], B)*gradN.transpose()*vec*m_gaussWeightHD[i];
    }

    H *= element.getDetJ()*m_mesh.getRefElementSize(m_mesh.getDim());

    return H;
}

void MatrixBuilder::setddev(Eigen::MatrixXd ddev)
{
    m_ddev = ddev;
}

void MatrixBuilder::setm(Eigen::VectorXd m)
{
    m_m = m;
}

void MatrixBuilder::setMcomputeFactor(simpleMatFuncElm computeFactor)
{
    m_Mfunc = computeFactor;
}

void MatrixBuilder::setMGammacomputeFactor(simpleMatFuncFacet computeFactor)
{
    m_MGammafunc = computeFactor;
}

void MatrixBuilder::setKcomputeFactor(matFuncElm computeFactor)
{
    m_Kfunc = computeFactor;
}

void MatrixBuilder::setDcomputeFactor(matFuncElm computeFactor)
{
    m_Dfunc = computeFactor;
}

void MatrixBuilder::setLcomputeFactor(matFuncElm computeFactor)
{
    m_Lfunc = computeFactor;
}

void MatrixBuilder::setCcomputeFactor(matFuncElm computeFactor)
{
    m_Cfunc = computeFactor;
}

void MatrixBuilder::setFcomputeFactor(matFuncElm computeFactor)
{
    m_Ffunc = computeFactor;
}

void MatrixBuilder::setHcomputeFactor(matFuncElm computeFactor)
{
    m_Hfunc = computeFactor;
}

void MatrixBuilder::setQFunc(qFuncFacet func)
{
    m_QFunc = func;
}
