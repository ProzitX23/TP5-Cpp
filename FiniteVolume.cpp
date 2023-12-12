#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh)
{
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
	// Matrix
	this->_mat_flux.resize(this->_msh->Get_triangles().size(),this->_msh->Get_triangles().size());
	// RHS
	this->_BC_RHS.resize(this->_msh->Get_triangles().size());
	this->_BC_RHS.setZero();
	vector<Edge> e=_msh->Get_edges();
	VectorXd l=_msh->Get_edges_length();
	VectorXd a=_msh->Get_triangles_area();
	Matrix<double,Dynamic,2> normale=_msh->Get_edges_normal();
	Matrix<double,Dynamic,2> center=_msh->Get_edges_center();
	Matrix<double,Dynamic,2> centertri=_msh->Get_triangles_center();
	vector<Triplet<double>> triplets;	triplets.clear();
	for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
	{
		int k  = e[i].Get_T1();
		int j  = e[i].Get_T2();
		double longueur = l(i);
		double alphaA; double betaA;
		double aire1=a(k);
		double x=center(i,0);
		double y=center(i,1);
		double Vx=_fct->Velocity_x(x,y,t);
		double Vy=_fct->Velocity_y(x,y,t);
		double nx=normale(i,0);
		double ny=normale(i,1);
		double x1=centertri(k,0);
		double y1=centertri(k,1);
		double vn=Vx*nx+Vy*ny;
		double distancecenter=sqrt(pow(x1-x,2)+pow(y1-y,2));
		double mu=_df->Get_mu();
		if (this->_df->Get_numerical_flux_choice() == "upwind")
		{
			if (vn<0)
			{
				alphaA = 0.;
				betaA = vn;
			}
			else
			{
				alphaA = vn;
				betaA = 0.;
			}
		}
		else if (this->_df->Get_numerical_flux_choice() == "centered")
		{
			alphaA = vn/2;
			betaA = vn/2;
		}
		if (j==-1)
		{
			if (e[i].Get_BC()=="Dirichlet")
			{
				double h=_fct->Dirichlet_Function(x, y, t);
				double alphaD=mu/(distancecenter);
				double betaD=-alphaD;
				triplets.push_back({k,k,((longueur*alphaD)/aire1)*(alphaA-betaA)+(longueur*alphaD)/aire1});
				_BC_RHS(k)+=2*betaA*(longueur/aire1)*h+((longueur*betaD)/aire1)*h;
			}
			if (e[i].Get_BC()=="Neumann")
			{
				double g=_fct->Neumann_Function(x, y, t);
				triplets.push_back({k,k,(longueur/aire1)*(alphaA+betaA)});
				_BC_RHS(k)+=-(mu*longueur*g)/aire1+(betaA*longueur*2*distancecenter*g)/aire1;
			}
		}
		else
		{
			double aire2=a(j);
			double x2=centertri(j,0);
			double y2=centertri(j,1);
			double distancecenter2=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
			double alphaD=mu/(distancecenter2);
			double alpha=alphaA+alphaD;
			double betaD=-alphaD;
			double beta=betaA+betaD;
			triplets.push_back({k,k,(longueur*alpha)/aire1});
			triplets.push_back({j,j,-(longueur*beta)/aire2});
			triplets.push_back({k,j,(longueur*beta)/aire1});
			triplets.push_back({j,k,-(longueur*alpha)/aire2});
		}
	}
	this->_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
	// cout <<"La matrice A est:" << endl;
	// cout << _mat_flux << endl;
	// cout << "Le vecteur B est:" << endl;
	// cout << _BC_RHS << endl;
}


// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1));
	}

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}
	// cout << sourceterm << endl;
	return sourceterm;
}

// Solution exacte au centre des triangles
VectorXd FiniteVolume::Exact_solution(const double t)
{
	VectorXd exactsol(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		exactsol(i) = this->_fct->Exact_solution(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}
	return exactsol;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, std::string st)
{
	double norm = 0;
	for (unsigned int i = 0; i < sol.rows(); i++)
	{
		norm += sol(i)*sol(i)*this->_msh->Get_triangles_area()[i];
	}
	norm = sqrt(norm);

	if (st == "solution")
	{
		cout << "Norme de u = " << norm << endl;
	}

	string name_file = this->_df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();
	assert(((long unsigned int)sol.size() == this->_msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol[i]) << endl;
	}
	solution << endl;

	//solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,this->_df->Get_dt()*fabs(sol[i])/this->_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	if (this->_df->Get_mu() > 1e-10)
	{
		solution << "SCALARS Pe float 1" << endl;
		solution << "LOOKUP_TABLE default" << endl;
		// To avoid strange behaviour (which appear only with Apple)
		// with Paraview when we have very small data (e-35 for example)
		for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
		{
			solution << max(eps,this->_msh->Get_triangles_length()(i)*fabs(sol[i])/this->_df->Get_mu()) << endl;
		}
		solution << endl;
	}

	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
