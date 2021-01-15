/// -*- c++ -*-

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>


#include "colvar.h"
#include "colvarbias_cvhd.h"


colvarbias_cvhd::colvarbias_cvhd()
  : colvarbias()
{
}


colvarbias_cvhd::colvarbias_cvhd(std::string const &conf, char const *key)
  : colvarbias(conf, key)
{
  get_keyval(conf, "hillWeight", hill_weight, 0.01);
  get_keyval(conf, "newHillFrequency", new_hill_freq, 1000);
  get_keyval(conf, "hillWidth", hill_width, std::sqrt(2.0 * PI) / 2.0);
  get_keyval(conf, "maxBias", max_bias, 1.0);
  get_keyval(conf, "wellTempered", isWellTempered, false);
  get_keyval(conf, "biasTemperature", bias_temp, -1.0);
  get_keyval(conf, "dynamicBias", isDynamic, false);
  get_keyval(conf, "hybridWaitTime", hybridWait, 1000);
  get_keyval(conf, "hybridBias", isHybrid, false);
  get_keyval(conf, "adaptiveEta", adaptive_eta, false);
  get_keyval(conf, "adiabatic", isAdiabatic, false);

  if (isHybrid) isDynamic = false;
  isStatic = !isDynamic;
  eta_old = 0;
  hybridStart = cvm::step_absolute();

  if (adaptive_eta) {
    eta_max = 0.1;
  } else {
    eta_max = 0.9;
  }

  if (hill_weight == 0.0)
    cvm::log("Warning: hillWeight has been set to zero, "
              "this bias will have no effect.\n");
  if ((bias_temp == -1.0) && isWellTempered) {
    cvm::fatal_error("Error: biasTemperature is not set.\n");
  }
  if (colvars.size() != 1) {
    cvm::fatal_error("Error: only one colvar allowed for CVHD!\n");
  }
  if (colvars[0]->value().type() != colvarvalue::type_scalar) {
    cvm::fatal_error("Error: CVHD colvar must be scalar!\n");
  }
}


colvarbias_cvhd::~colvarbias_cvhd()
{
  if (cvm::n_cvhd_biases > 0)
    cvm::n_cvhd_biases -= 1;
}

/*
 * Calculate forces
 *
 * Note that stuff like waiting times and "moving" of the CV, is
 * handled by the CV codes themselves. This might not be ideal
 * (and not really modular), but it allows to use colvars.traj
 * as a simple way to monitor the "progress" of the system.
 */
cvm::real colvarbias_cvhd::update()
{
  if (cvm::debug())
    cvm::log("Updating the CVHD bias \""+this->name+"\""+".\n");

  // check position of CV
  cvm::real eta = colvars[0]->value();
  if (eta < 0) eta = 0;
  int eta_curr = std::floor(eta/2.0);
  if (eta_curr > eta_old) {
    eta_old = eta_curr;
    hills.clear();
    if (adaptive_eta) eta_max = 0.1;
    if (isHybrid) hybridStart = cvm::step_absolute();
    if (isHybrid) isDynamic = false;
    if (isAdiabatic) bias_energy = -1.0;
    if (isAdiabatic) return bias_energy;
  }
  eta = eta - 2.0*eta_curr;

  if (isHybrid && (cvm::step_absolute() - hybridStart) > hybridWait) isDynamic = true;

  // calculate the biasing energy and forces
  bias_energy = 0.0;
  colvar_forces[0].reset();

  if (isStatic && eta < 1.0) {
    bias_energy = max_bias*(1.0 - eta);
    colvar_forces[0].real_value += max_bias;
  }

  if (isDynamic) {
    // check whether moving eta max must be updated
    cvm::real sigma = 0.5*hill_width*colvars[0]->width;
    if (adaptive_eta && eta > eta_max+4*sigma) eta_max = eta-4*sigma;

    // manage addition of hills
    if ((cvm::step_absolute() % new_hill_freq) == 0 && eta < eta_max) {
      hill newHill;
      if (isWellTempered) {
        cvm::real bias_curr = calc_hills();
        cvm::real expWeight = std::exp(-1.0*bias_curr/(bias_temp*cvm::boltzmann()));
        newHill.weight = expWeight*hill_weight;
      } else {
        newHill.weight = hill_weight;
      }
      newHill.center = colvars[0]->value();
      newHill.width  = hill_width*colvars[0]->width;
      hills.push_back(newHill);
    }

    // calculation of energy
    bias_energy += calc_hills();
    calc_hills_force();
  }

  return bias_energy;
}

/*
 * This is a re-implementation of metadynamics, but much simpler than
 * the one already present in Colvars.
 */

cvm::real colvarbias_cvhd::calc_hills()
{
  cvm::real energy = 0.0;
  for (int i = 0; i < hills.size(); i++) {
    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    colvarvalue const &x       = colvars[0]->value();
    colvarvalue const &center  = hills[i].center;
    cvm::real const half_width = 0.5 * hills[i].width;
    cv_sqdev += (colvars[0]->dist2(x, center)) / (half_width*half_width);

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-05)
      hills[i].value = 0.0;
    } else {
      hills[i].value = std::exp(-0.5*cv_sqdev);
    }
    energy += hills[i].value*hills[i].weight;
  }
  return energy;
}

void colvarbias_cvhd::calc_hills_force()
{
  colvarvalue const &x = colvars[0]->value();
  for (int i = 0; i < hills.size(); i++) {
    if (hills[i].value == 0.0) continue;
    colvarvalue const &center  = hills[i].center;
    cvm::real const half_width = 0.5 * hills[i].width;
    colvar_forces[0].real_value +=
      ( hills[i].weight * hills[i].value * (0.5 / (half_width*half_width)) *
      (colvars[0]->dist2_lgrad(x, center)).real_value );
    }
}
