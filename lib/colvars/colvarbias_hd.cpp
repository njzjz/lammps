/// -*- c++ -*-

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>


#include "colvar.h"
#include "colvarbias_hd.h"


colvarbias_hd::colvarbias_hd()
  : colvarbias()
{
}


colvarbias_hd::colvarbias_hd(std::string const &conf, char const *key)
  : colvarbias(conf, key)
{
  get_keyval(conf, "maxBias", max_bias, 1.0);
  
  if (colvars.size() != 1) {
    cvm::fatal_error("Error: only one colvar allowed for hyperdynamics!\n");
  }
  
  if (colvars[0]->value().type() != colvarvalue::type_scalar) {
    cvm::fatal_error("Error: hyperdynamics colvar must be scalar!\n");
  }
}


colvarbias_hd::~colvarbias_hd()
{
  if (cvm::n_hd_biases > 0)
    cvm::n_hd_biases -= 1;
}

cvm::real colvarbias_hd::update()
{
  if (cvm::debug())
    cvm::log("Updating the hyperdynamics bias \""+this->name+"\""+".\n");

  // calculate the biasing energy and forces
  bias_energy = 0.0;
  colvar_forces[0].reset();

  cvm::real xi = colvars[0]->value();
  xi = xi - 2.0*floor(xi/2.0);

  if (xi < 1.0) {
    bias_energy = max_bias*(1.0 - xi);
    colvar_forces[0].real_value += max_bias;
  }

  return bias_energy;
}

