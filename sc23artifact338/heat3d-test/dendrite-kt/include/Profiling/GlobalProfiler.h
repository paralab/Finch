//
// Created by maksbh on 6/16/21.
//

#pragma once
#if BUILD_WITH_PETSC
#ifdef PETSC_LOGGING_IMPL
#define DendriteEvent  int
#else
#define DendriteEvent extern int
#endif

namespace Profiling{
  DendriteEvent matAssembly, matElementalAssembly, matBC, matGaussPointAssembly;
  DendriteEvent vecAssembly, vecElementalAssembly, vecBC, vecGaussPointAssembly;
}




#ifdef PETSC_LOGGING_IMPL
static void registerGlobalProfiler()
{
  PetscLogEventRegister("Dkt-matAssembly", 0, &Profiling::matAssembly);
  PetscLogEventRegister("Dkt-matElementalAssembly", 0, &Profiling::matElementalAssembly);
  PetscLogEventRegister("Dkt-matBC", 0, &Profiling::matBC);
  PetscLogEventRegister("Dkt-matGaussPointAssembly", 0, &Profiling::matGaussPointAssembly);

  PetscLogEventRegister("Dkt-vecAssembly", 0, &Profiling::vecAssembly);
  PetscLogEventRegister("Dkt-vecElementalAssembly", 0, &Profiling::vecElementalAssembly);
  PetscLogEventRegister("Dkt-vecBC", 0, &Profiling::vecBC);
  PetscLogEventRegister("Dkt-vecGaussPointAssembly", 0, &Profiling::vecGaussPointAssembly);


}
#endif
#undef DendriteEvent
#endif
