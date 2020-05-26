#define SOLVER_API
#define MESH_API

%include std_string.i
%include std_vector.i

%apply const std::string& {std::string* mshFile};
namespace std
{
  %template(VectorDouble) vector<double>;
  %template(VectorString) vector<string>;
  %template(VectorVectorDouble) vector<vector<double>>;
};

%module pfemSolverw
%{
	/* Includes the header in the wrapper code */
	#include "Solver.hpp"
	#include "SolverIncompressible.hpp"
	#include "SolverCompressible.hpp"
%}

/* Parse the header file to generate wrappers */
%include "../mesh/Mesh.hpp"
%include "Solver.hpp"
%include "SolverIncompressible.hpp"
%include "SolverCompressible.hpp"
