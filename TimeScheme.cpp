#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_fin_vol(adv),_df(data_file), _sol(adv->Initial_condition()), _t(_df->Get_t0())
{
}

EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
   std::cout << "Build time scheme class." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Euler Explicite
void EulerScheme::Advance()
{
   _fin_vol->Build_flux_mat_and_rhs(this->_t);
   SparseMatrix<double> A=_fin_vol->Get_flux_matrix();
   VectorXd B=_fin_vol->Get_BC_RHS();
   double _dt=_df->Get_dt();
   SparseMatrix<double> Id(A.rows(),A.cols());
   Id.setIdentity();
   vector<Triplet<double>> triplets;	triplets.clear();
   VectorXd S=_fin_vol->Source_term(this->_t);
   this->_sol = (Id-_dt*A)*this->_sol+_dt*(S-B);
   this->_t += _dt;
}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
   double dt = this->_df->Get_dt();
   this->_fin_vol->Build_flux_mat_and_rhs(this->_t+dt);
   SparseMatrix<double> A =this->_fin_vol->Get_flux_matrix();
   VectorXd B=_fin_vol->Get_BC_RHS();
   VectorXd S=_fin_vol->Source_term(this->_t+dt);
   SparseMatrix<double> Id(A.rows(),A.cols());
   Id.setIdentity();
   A=Id+dt*A;
   SparseLU<SparseMatrix<double>> solver;
   B = this->_sol -dt*B + dt*S;
   solver.compute(A);
   this->_sol = solver.solve(B);
   this->_t = this->_t + dt;
}

#define _TIME_SCHEME_CPP
#endif
