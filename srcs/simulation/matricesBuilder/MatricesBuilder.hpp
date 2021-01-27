#pragma once
#ifndef MATRIXBUILDER_HPP_INCLUDED
#define MATRIXBUILDER_HPP_INCLUDED

#include <functional>
#include <Eigen/Dense>

class Element;
class Facet;
class Mesh;

using matFuncElm = std::function<double(const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/)>;
using simpleMatFuncElm = std::function<double(const Element& /** element **/, const Eigen::MatrixXd& /** N **/)>;
using simpleMatFuncFacet = std::function<double(const Facet& /** facet **/, const Eigen::MatrixXd& /** N **/)>;
using qFuncFacet = std::function<Eigen::VectorXd(const Facet& /** facet **/, const std::array<double, 3>& /** gp **/)>;

/**
 * \class MatrixBuilder
 * \brief Class responsible to hold code for building matrices.
 */
class MatrixBuilder
{
	public:
		MatrixBuilder(const Mesh& mesh, unsigned int nGPHD, unsigned int nGPLD);
		~MatrixBuilder();

		Eigen::MatrixXd getB(Eigen::MatrixXd gradN);
		Eigen::MatrixXd getGradN(const Element& element);
		Eigen::MatrixXd getM(const Element& element);
		Eigen::MatrixXd getMGamma(const Facet& facet);
		Eigen::MatrixXd getK(const Element& element, const Eigen::MatrixXd& B);
		Eigen::MatrixXd getD(const Element& element, const Eigen::MatrixXd& B);
		Eigen::MatrixXd getL(const Element& element, const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN);
		Eigen::MatrixXd getC(const Element& element, const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN);
		Eigen::VectorXd getF(const Element& element, const Eigen::VectorXd& vec, const Eigen::MatrixXd& B);
		Eigen::VectorXd getH(const Element& element, const Eigen::VectorXd& vec, const Eigen::MatrixXd& B, const Eigen::MatrixXd& gradN);
		Eigen::VectorXd getQN(const Facet& facet);

		void setddev(Eigen::MatrixXd ddev);
		void setm(Eigen::VectorXd m);
		void setMcomputeFactor(simpleMatFuncElm computeFactor);
		void setMGammacomputeFactor(simpleMatFuncFacet computeFactor);
		void setKcomputeFactor(matFuncElm computeFactor);
		void setDcomputeFactor(matFuncElm computeFactor);
		void setLcomputeFactor(matFuncElm computeFactor);
		void setCcomputeFactor(matFuncElm computeFactor);
		void setFcomputeFactor(matFuncElm computeFactor);
		void setHcomputeFactor(matFuncElm computeFactor);
		void setQFunc(qFuncFacet func);

		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> lump2(const Eigen::MatrixXd& mat)
		{
            if(mat.cols() != mat.rows())
                throw std::runtime_error("only square matrices can be lumped!");

            Eigen::DiagonalMatrix<double, Eigen::Dynamic> lumpedMat;
            lumpedMat.resize(mat.rows());
            lumpedMat.setZero();

            auto& lumpedMatDiag = lumpedMat.diagonal();

            for(auto i = 0; i < mat.cols() ; ++i)
            {
                for(auto j = 0; j < mat.cols() ; ++j)
                    lumpedMatDiag(i) += mat(i,j);
            }

            return lumpedMat;
        }

        static void lump(Eigen::MatrixXd& mat)
		{
            if(mat.cols() != mat.rows())
                throw std::runtime_error("only square matrices can be lumped!");

            for(auto i = 0; i < mat.cols() ; ++i)
            {
                for(auto j = 0; j < mat.cols() ; ++j)
                {
                    if(i != j)
                    {
                        mat(i, i) += mat(i,j);
                        mat(i, j) = 0;
                    }
                }
            }
        }

		static Eigen::MatrixXd diagBlock(const Eigen::MatrixXd& mat, unsigned int nBlocks)
        {
            if(mat.cols() != mat.rows())
                throw std::runtime_error("only square matrices can be diagonalized!");

            Eigen::MatrixXd finalMat(mat.cols()*nBlocks, mat.cols()*nBlocks); finalMat.setZero();

            for(unsigned int i = 0 ; i < nBlocks ; ++i)
                finalMat.block(i*mat.cols(), i*mat.cols(), mat.cols(), mat.cols()) = mat;

            return finalMat;
        }

        static void inverse(Eigen::DiagonalMatrix<double, Eigen::Dynamic>& mat)
        {
            auto& matDiag = mat.diagonal();

            #pragma omp parallel for default(shared)
            for(auto i = 0 ; i < mat.rows() ; ++i)
                matDiag[i] = 1/matDiag[i];
        }

    private:
        const Mesh& m_mesh;

        unsigned int m_nGPHD; /**< Number of Gauss points use for the elements. **/
        std::vector<double> m_gaussWeightHD; /**< Gauss weights for the elements. **/
        std::vector<std::array<double, 3>> m_gaussPointsHD; /**< Gauss points for the elements. **/
        std::vector<std::vector<double>> m_sfsHD; /**< Element shape functions evaluated at each Gauss Points for elements. **/
        std::vector<std::vector<double>> m_gradsfsHD; /**< Element gradient of shape functions for elements (constant). **/

        unsigned int m_nGPLD; /**< Number of Gauss points use for the facets. **/
        std::vector<double> m_gaussWeightLD; /**< Gauss weights for the facets. **/
        std::vector<std::array<double, 3>> m_gaussPointsLD; /**< Gauss points for the facets. **/
        std::vector<std::vector<double>> m_sfsLD; /**< Element shape functions evaluated at each Gauss Points for facets. **/
        std::vector<std::vector<double>> m_gradsfsLD; /**< Element gradient of shape functions for facets (constant). **/

        std::vector<Eigen::MatrixXd> m_NHD;
        std::vector<Eigen::MatrixXd> m_NLD;
        std::vector<Eigen::MatrixXd> m_NHDtilde;
        std::vector<Eigen::MatrixXd> m_NLDtilde;

        Eigen::MatrixXd m_ddev;
        Eigen::VectorXd m_m;

        simpleMatFuncElm m_Mfunc;
        simpleMatFuncFacet m_MGammafunc;
        matFuncElm m_Kfunc;
        matFuncElm m_Dfunc;
        matFuncElm m_Lfunc;
        matFuncElm m_Cfunc;
        matFuncElm m_Ffunc;
        matFuncElm m_Hfunc;
        qFuncFacet m_QFunc;
};

#endif // MATRIXBUILDER_HPP_INCLUDED
