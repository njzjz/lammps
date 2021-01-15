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

#ifdef FIX_CLASS

FixStyle(timeboost,FixTimeboost)

#else

#ifndef LMP_FIX_TIMEBOOST_H
#define LMP_FIX_TIMEBOOST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTimeboost : public Fix {
 public:
  FixTimeboost(class LAMMPS *, int, char **);
  ~FixTimeboost();
  int setmask();
  virtual void init();
  virtual void end_of_step();

 private:
  char *id_fix;
  class Fix *fix_energy;
  double tset;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
