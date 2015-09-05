/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "TestADKernel.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

// libmesh includes
#include "libmesh/threads.h"

template<>
InputParameters validParams<TestADKernel>()
{
  InputParameters params = validParams<KernelBase>();
  params.registerBase("Kernel");
  return params;
}

TestADKernel::TestADKernel(const InputParameters & parameters) :
    KernelBase(parameters),
    _u_dot(_var.uDot()),
    _du_dot_du(_var.duDotDu())
{
}

TestADKernel::~TestADKernel()
{
}

ADReal
TestADKernel::computeQpResidual()
{
  return _grad_u[_qp] * _grad_test[_i][_qp];
}

void
TestADKernel::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();


  if (true) //_need_u)
  {
    std::vector<dof_id_type> _dof_indices = _var.dofIndices();
    dof_id_type num_dofs = _dof_indices.size();

    unsigned int nqp = _qrule->n_points();

    unsigned int _var_num = _var.number();

    const NumericVector<Real> & current_solution = *_var.sys().currentSolution();

    _ad_dofs.resize(num_dofs);
    _u.resize(nqp);

    const VariablePhiValue & phi = _assembly.fePhi(_var.feType());
    const VariablePhiGradient & grad_phi = _assembly.feGradPhi(_var.feType());

    if (true) //_need_grad_u)
      _grad_u.resize(nqp);

    /*
    if (false) //_need_ad_second_u)
      _ad_second_u.resize(nqp);
    */

    // Derivatives are offset by the variable number
    size_t ad_offset = _var_num * num_dofs;

    // Hopefully this problem can go away at some point
    if (ad_offset + num_dofs > AD_MAX_DOFS_PER_ELEM)
      mooseError("Current number of dofs per element is greater than AD_MAX_DOFS_PER_ELEM of " << AD_MAX_DOFS_PER_ELEM);

    for (unsigned int qp=0; qp < nqp; qp++)
    {
      _u[qp] = 0;

      if (true) //_need_grad_u)
        _grad_u[qp] = 0;

      /*
      if (_need_ad_second_u)
        _ad_second_u[qp] = 0;
      */
    }

    for (unsigned int i=0; i < num_dofs; i++)
    {
      _ad_dofs[i] = current_solution(_dof_indices[i]);

      // NOTE!  You have to do this AFTER setting the value!

      // NOTE!!! You MUST switch this to the "insert" implementation when using sparse storage!
      _ad_dofs[i].derivatives()[ad_offset + i] = 1.0;
      //_ad_dofs[i].derivatives().insert(ad_offset + i) = 1.0;
    }

    // Now build up the solution at each quadrature point:
    for (unsigned int i=0; i < num_dofs; i++)
    {
      for (unsigned int qp=0; qp < nqp; qp++)
      {
        _u[qp] += _ad_dofs[i] * phi[i][qp];

        if (true) //_need_grad_u)
          _grad_u[qp].add_scaled(grad_phi[i][qp], _ad_dofs[i]); // Note: += does NOT work here!

        /*
        if (_need_ad_second_u)
          _ad_second_u[qp].add_scaled((*_second_phi)[i][qp], _ad_dofs[i]); // Note: += does NOT work here!
        */
      }
    }
  }


  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual().value();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_save_in.size(); i++)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
TestADKernel::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  size_t ad_offset = _var.number() * _test.size();

  for (_i = 0; _i < _test.size(); _i++)
  {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
      ADReal residual = computeQpResidual(); // This will also compute the derivative with respect to all dofs
      for (_j = 0; _j < _phi.size(); _j++)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * residual.derivatives()[ad_offset + _j];
    }
  }

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i=0; i<rows; i++)
      diag(i) = _local_ke(i,i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i=0; i<_diag_save_in.size(); i++)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
TestADKernel::computeOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var.number())
    computeJacobian();
  else
  {
    size_t ad_offset = jvar * _test.size();

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);

    for (_i = 0; _i < _test.size(); _i++)
    {
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
        ADReal residual = computeQpResidual(); // This will also compute the derivative with respect to all dofs

        for (_j = 0; _j < _phi.size(); _j++)
          ke(_i, _j) += _JxW[_qp] * _coord[_qp] * residual.derivatives()[ad_offset + _j];
      }
    }
  }
}

void
TestADKernel::computeOffDiagJacobianScalar(unsigned int jvar)
{
  /*
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < jv.order(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpOffDiagJacobian(jvar);
  */
}
