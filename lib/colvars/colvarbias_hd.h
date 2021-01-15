/// -*- c++ -*-

#ifndef COLVARBIAS_HD_H
#define COLVARBIAS_HD_H

#include <vector>
#include <list>
#include <sstream>
#include <fstream>

#include "colvarbias.h"

/// Hyperdynamics bias (implementation of \link colvarbias \endlink)
class colvarbias_hd : public colvarbias {

public:

  /// Constructor
  colvarbias_hd(std::string const &conf, char const *key);

  /// Default constructor
  colvarbias_hd();

  /// Destructor
  virtual ~colvarbias_hd();

  virtual cvm::real update();

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart(std::istream &is) {}

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart(std::ostream &os) {}

  /// Write a label to the trajectory file (comment line)
  //virtual std::ostream & write_traj_label(std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  //virtual std::ostream & write_traj(std::ostream &os);

protected:

  cvm::real max_bias;
};



#endif
