#ifndef METAAPP_H
#define METAAPP_H

#include "MooseApp.h"

class MetaApp;

template<>
InputParameters validParams<MetaApp>();

class MetaApp : public MooseApp
{
public:
  MetaApp(InputParameters parameters);
  virtual ~MetaApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* METAAPP_H */
