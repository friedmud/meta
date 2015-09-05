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

#ifndef TESTADKERNEL_H
#define TESTADKERNEL_H

#include "KernelBase.h"


// MetaPhysicL
#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"
#include "metaphysicl/dynamicsparsenumberarray.h"

// The 100 here is for how many DoFs there are per element.
#define AD_MAX_DOFS_PER_ELEM 100

typedef MetaPhysicL::DualNumber<double, MetaPhysicL::NumberArray<AD_MAX_DOFS_PER_ELEM, double> > ADReal;

// NOTE!  if you switch to sparse storage you MUST switch to using "insert" over in the .C file around line 110
//typedef MetaPhysicL::DualNumber<double, MetaPhysicL::DynamicSparseNumberArray<double, dof_id_type> > ADReal;

namespace libMesh
{

template<>
struct CompareTypes<double, ADReal>
{
  typedef ADReal supertype;
};

template<>
struct CompareTypes<ADReal, double>
{
  typedef ADReal supertype;
};

}

typedef VectorValue<ADReal> ADRealVectorValue;
typedef ADRealVectorValue ADRealGradient;

typedef TensorValue<ADReal> ADRealTensorValue;
typedef ADRealTensorValue ADRealTensor;

typedef MooseArray<ADReal>             ADVariableValue;
typedef MooseArray<ADRealGradient>     ADVariableGradient;
typedef MooseArray<ADRealTensor>       ADVariableSecond;

class TestADKernel;

template<>
InputParameters validParams<TestADKernel>();

class TestADKernel :
  public KernelBase
{
public:
  TestADKernel(const InputParameters & parameters);

  virtual ~TestADKernel();

  // See KernelBase base for documentation of these overridden methods
  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);

protected:
  /// Compute this Kernel's contribution to the residual at the current quadrature point
  virtual ADReal computeQpResidual();

  /// Holds the solution at current quadrature points
  ADVariableValue _u;

  /// Holds the solution gradient at the current quadrature points
  ADVariableGradient _grad_u;

  /// Time derivative of u
  VariableValue & _u_dot;

  /// Derivative of u_dot with respect to u
  VariableValue & _du_dot_du;

  std::vector<ADReal> _ad_dofs;
};

#endif /* TESTADKERNEL_H */
