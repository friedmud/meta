#include "MetaApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"
#include "metaphysicl/dynamicsparsenumberarray.h"

// The 100 here is for how many DoFs there are per element.
#define AD_MAX_DOFS_PER_ELEM 100
typedef MetaPhysicL::DualNumber<double, MetaPhysicL::NumberArray<AD_MAX_DOFS_PER_ELEM, double> > ADReal;
//typedef MetaPhysicL::DualNumber<double, MetaPhysicL::DynamicSparseNumberArray<double, dof_id_type> > ADReal;

// Create a performance log
PerfLog Moose::perf_log("Meta");

// Begin the main program.
int main(int argc, char *argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  MetaApp::registerApps();

  // This creates dynamic memory that we're responsible for deleting
  MooseApp * app = AppFactory::createApp("MetaApp", argc, argv);

  // Execute the application
  app->run();

  // Free up the memory we created earlier
  delete app;

  return 0;
}
