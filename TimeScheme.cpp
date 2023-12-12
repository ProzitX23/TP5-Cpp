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

void TimeScheme::Moyenne()
{  
   double Moy;
   double j(_sol.size());
   for (int i(0);i<_sol.size();i++)
   {
      Moy=Moy+_sol(i);
   }
   Moy=Moy/j;
   cout <<"La température moyenne est:"<< Moy << endl;
}

// Euler Explicite
void EulerScheme::Advance()
{
   double _dt=this->_df->Get_dt();
   _fin_vol->Build_flux_mat_and_rhs(this->_t);
   SparseMatrix<double> A=_fin_vol->Get_flux_matrix();
   VectorXd B=_fin_vol->Get_BC_RHS();
   SparseMatrix<double> Id(A.rows(),A.cols());
   Id.setIdentity();
   VectorXd S=_fin_vol->Source_term(this->_t);
   this->_sol = this->_sol -_dt*(A*this->_sol+B)+_dt*S;
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
