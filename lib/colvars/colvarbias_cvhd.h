/// -*- c++ -*-

#ifndef COLVARBIAS_CVHD_H
#define COLVARBIAS_CVHD_H

#include <vector>
#include <list>
#include <sstream>
#include <fstream>

#include "colvarbias.h"

struct hill {
  cvm::real weight;
  cvm::real width;
  cvm::real value;
  colvarvalue center;
};

/// Hyperdynamics bias (implementation of \link colvarbias \endlink)
class colvarbias_cvhd : public colvarbias {

public:

  /// Constructor
  colvarbias_cvhd(std::string const &conf, char const *key);

  /// Default constructor
  colvarbias_cvhd();

  /// Destructor
  virtual ~colvarbias_cvhd();

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
  cvm::real hill_width;
  cvm::real hill_weight;
  cvm::real bias_temp;
  cvm::real eta_max;
  size_t new_hill_freq;
  int eta_old;
  int hybridWait;
  int hybridStart;

  std::vector<hill> hills;

  bool isDynamic;
  bool isStatic;
  bool isWellTempered;
  bool adaptive_eta;
  bool isHybrid;
  bool isAdiabatic;

  cvm::real calc_hills();
  void calc_hills_force();
};

#endif
