/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "math.h"
#include "fix_timeboost.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixTimeboost::FixTimeboost(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix timeboost command");

  global_freq = 1;
  
  tset = force->numeric(FLERR,arg[3]);
  id_fix = new char[strlen(arg[4])+1];
  strcpy(id_fix,arg[4]);
  
  if (tset <= 0) error->all(FLERR,"Illegal fix timeboost command");
}

/* ---------------------------------------------------------------------- */

FixTimeboost::~FixTimeboost()
{
  delete[] id_fix;
}

/* ---------------------------------------------------------------------- */

int FixTimeboost::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTimeboost::init()
{
  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find fix ID for fix timeboost");
  fix_energy = modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

void FixTimeboost::end_of_step()
{
  double dt = update->dt;
  double boltz = force->boltz;
  double eboost = fix_energy->compute_scalar();
  double boost = exp(eboost/(boltz*tset));
  update->atime += (boost-1.0)*dt;
}

