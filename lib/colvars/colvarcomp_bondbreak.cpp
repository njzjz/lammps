/// -*- c++ -*-

#include <cmath>
#include <vector>

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



template<bool calculate_gradients>
cvm::real colvar::bondbreak::switching_function(cvm::real const &rmin,
                                                cvm::real const &rmax,
                                                cvm::real const &power,
                                                cvm::real const &prefactor,
                                                cvm::atom &A1,
                                                cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::real const dist = diff.norm();

  if (dist > rmax) {
    // cvm::log("Bond broken between " + cvm::to_str(A1.id) + " and " + cvm::to_str(A2.id) + ", length: " + cvm::to_str(dist) + ".\n");
    return cvm::real(1.0);
  }
  if (dist < rmin) return cvm::real(0.0);

  cvm::real const xi = (dist - rmin)/(rmax - rmin);

  // gradients
  if (calculate_gradients) {
    cvm::real const dFdr = prefactor*pow(xi, power-1.0)/(rmax-rmin);
    cvm::rvector const drdx = diff/dist;
    A1.grad += (-1.0)*dFdr*drdx;
    A2.grad +=        dFdr*drdx;
  }

  return xi;
}

colvar::bondbreak::bondbreak(std::string const &conf)
  : distance(conf)
{
  function_type = "bondbreak";
  x.type(colvarvalue::type_scalar);

  // need to specify this explicitly because the distance() constructor
  // has set it to true
  b_inverse_gradients = false;
  oneGroup = false;
  rebuildList = true;
  purgeList = false;
  offset = 0.0;
  pairsum = 0.0;
  laststep = 0;
  jumpstep = 0;
  progress = 0.0;

  get_keyval(conf, "rmin", rmin,
              cvm::real(-1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmax", rmax,
              cvm::real(2.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rcut", rcut,
              cvm::real(-1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "power", power,
              cvm::real(12.0));

  get_keyval(conf, "waitTime", nwait,
              int(100));

  get_keyval(conf, "fastEventTime", fastEventTresh,
              int(-1));

  if (fastEventTresh < 0) fastEvent = false;
  if (fastEventTresh > 0) fastEvent = true;

  if (rcut < 0) rcut = rmax;

}

colvar::bondbreak::bondbreak()
{
  function_type = "bondbreak";
  x.type(colvarvalue::type_scalar);
}

void colvar::bondbreak::build_bondlist()
{
  rebuildList = false;
  cleanList = true;
  cleanstep = cvm::step_absolute();
  cvm::real const cut2 = rcut*rcut;
  cvm::real r2;

  pairList1.clear();
  pairList2.clear();
  pairList3.clear();

  // assume that first atom same == all atoms same (might turn this in to an option)
  if (group1[0].pos.x == group2[0].pos.x &&
      group1[0].pos.y == group2[0].pos.y &&
      group1[0].pos.z == group2[0].pos.z) {
    oneGroup = true;
  }

  // rmin < 0: use average value of BL as ref
  // this average is computed in the cleanList routine
  if (oneGroup) {
    for (size_t i = 0; i < group1.size() - 1; i++) {
      for (size_t j = i + 1; j < group1.size(); j++) {
        r2 = cvm::position_dist2(group1[i].pos, group1[j].pos);
        if (r2 < cut2) {
          pairList1.push_back(i);
          pairList2.push_back(j);
          // cvm::log("Bond found: " + cvm::to_str(group1[i].id) + " and " + cvm::to_str(group1[j].id) + ", length: " + cvm::to_str(sqrt(r2)) + ".\n");
          if (rmin < 0.0) {
            pairList3.push_back(0.0);
          } else {
            pairList3.push_back(rmin);
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i < group1.size(); i++) {
      for (size_t j = 0; j < group2.size(); j++) {
        r2 = cvm::position_dist2(group1[i].pos, group2[j].pos);
        if (r2 < cut2) {
          pairList1.push_back(i);
          pairList2.push_back(j);
          if (rmin < 0.0) {
            pairList3.push_back(0.0);
          } else {
            pairList3.push_back(rmin);
          }
        }
      }
    }
  }
}

void colvar::bondbreak::calc_value()
{
  // bookkeeping: check whether bond list rebuild is needed.
  // also, keep in mind that using '==' for floating point numbers
  // is a bad idea, so add a small tolerance for progress check
  // cases:
  // pairsum ~= 1: start/keep counting how long this takes (and eventually reset lists)
  //   -> there's also a fast event check included here
  // laststep == current step: we already did this step, quit
  // else, we increment laststep in order to have a moving frame of reference when the first case starts
  if (pairsum > 0.9999) {
    if ((cvm::step_absolute() - laststep) == nwait) {
      laststep = cvm::step_absolute();
      if (fastEvent && cvm::step_absolute()-jumpstep < fastEventTresh) {
        rebuildList = false;
        purgeList = true;
        jumpstep = cvm::step_absolute();
      } else {
        rebuildList = true;
        jumpstep = cvm::step_absolute();
        offset++;
      }
    }
  } else if (laststep == cvm::step_absolute()) {
    x.real_value = 2.0*cvm::real(offset) + progress;
    return;
  } else {
    laststep = cvm::step_absolute();
  }

  // construct a new bondlist if requested
  if (rebuildList) build_bondlist();

  // after the new bondlist is constructed, check during nwait steps
  // whether each bond is actually a bond
  // we set the CV as something we haven't had before (no bias there yet + new bias not used later)
  // we set progress = 1 (no forces)
  // pairsum = 0 (do not trigger a bondlist rebuild as checked at above)
  // We move *backwards* through the list b/c stuff is going to be deleted
  // (use int i, not tagint, b/c i should become -1 at some point)
  if (cleanList) {
    if (cvm::step_absolute() - cleanstep < nwait) {
      x.real_value = 2.0*cvm::real(offset) - 0.5;
      progress = 1.0;
      pairsum = 0.0;
      cvm::real cut2 = rcut*rcut;
      cvm::real r2;
      for (int i = pairList1.size()-1; i >= 0; i--) {
        if (oneGroup) {
          r2 = cvm::position_dist2(group1[pairList1[i]].pos, group1[pairList2[i]].pos);
        } else {
          r2 = cvm::position_dist2(group1[pairList1[i]].pos, group2[pairList2[i]].pos);
        }
        if (rmin < 0) {
          pairList3[i] += sqrt(r2)/nwait;
        }
        if (r2 > cut2) {
          pairList1.erase(pairList1.begin() + i);
          pairList2.erase(pairList2.begin() + i);
          pairList3.erase(pairList3.begin() + i);
        } 
      }
      return;
    } else {
      cleanList = false;
    }
  }

  if (purgeList) {
    cvm::real cut2 = rmax*rmax;
    cvm::real r2;
    purgeList = false;
    for (int i = pairList1.size()-1; i >= 0; i--) {
      if (oneGroup) {
        r2 = cvm::position_dist2(group1[pairList1[i]].pos, group1[pairList2[i]].pos);
      } else {
        r2 = cvm::position_dist2(group1[pairList1[i]].pos, group2[pairList2[i]].pos);
      }
      if (r2 > cut2) {
        pairList1.erase(pairList1.begin() + i);
        pairList2.erase(pairList2.begin() + i);
        pairList3.erase(pairList3.begin() + i);
      }
    }
  }

  cvm::real pairterm;
  cvm::real dummy = 0.0;
  pairsum = 0.0;
  x.real_value = 0.0;

  // get all pairwise terms
  for (size_t i = 0; i < pairList1.size(); i++) {
    if (oneGroup) {
      pairterm = switching_function<false> (pairList3[i], rmax, power, dummy, group1[pairList1[i]], group1[pairList2[i]]);
    } else {
      pairterm = switching_function<false> (pairList3[i], rmax, power, dummy, group1[pairList1[i]], group2[pairList2[i]]);
    }
    pairsum += pow(pairterm, power); 
  }

  // finish colvar calculation
  progress = pow(pairsum, 1.0/power);
  if (progress < 1.0) {
    progress = 0.5*(1.0-cos(PI*progress*progress));
  } else {
    progress = 1.0;
  }
  x.real_value = 2.0*cvm::real(offset) + progress;
}


void colvar::bondbreak::calc_gradients()
{
  // need to get full sum of pairwise terms before pairwise forces can be calculated
  calc_value();
  if (progress < 1.0) {
    cvm::real xi_t = pow(pairsum, 1.0/power);
    cvm::real const prefact = PI*xi_t*sin(PI*xi_t*xi_t)
                                *pow(pairsum, 1.0/power-1.0);

    for (size_t i = 0; i < pairList1.size(); i++) {
      if (oneGroup) {
        switching_function<true> (pairList3[i], rmax, power, prefact, group1[pairList1[i]], group1[pairList2[i]]);
      } else {
        switching_function<true> (pairList3[i], rmax, power, prefact, group1[pairList1[i]], group2[pairList2[i]]);
      }
    }
  }
}

void colvar::bondbreak::apply_force(colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force(force.real_value);

  if (!group2.noforce && !oneGroup)
    group2.apply_colvar_force(force.real_value);
}

cvm::real colvar::angleswitch::ang_val(cvm::atom &a1, cvm::atom &a2, cvm::atom &a3, cvm::atom &a4)
{
  // Usual sign convention: r12 = r2 - r1
  r12 = cvm::position_distance(a1.pos, a2.pos);
  r23 = cvm::position_distance(a2.pos, a3.pos);
  r34 = cvm::position_distance(a3.pos, a4.pos);

  cvm::rvector const n1 = cvm::rvector::outer(r12, r23);
  cvm::rvector const n2 = cvm::rvector::outer(r23, r34);

  cvm::real const cos_phi = n1 * n2;
  cvm::real const sin_phi = n1 * r34 * r23.norm();

  cvm::real ang = (180.0/PI) * std::atan2(sin_phi, cos_phi);
  if (ang > 180)  ang -= 360;
  if (ang < -180) ang += 360;

  return ang;
}

void colvar::angleswitch::ang_val_grad(cvm::rvector &F1, cvm::rvector &F2, cvm::rvector &F3, cvm::rvector &F4)
{
  cvm::rvector A = cvm::rvector::outer(r12, r23);
  cvm::real   rA = A.norm();
  cvm::rvector B = cvm::rvector::outer(r23, r34);
  cvm::real   rB = B.norm();
  cvm::rvector C = cvm::rvector::outer(r23, A);
  cvm::real   rC = C.norm();

  cvm::real const cos_phi = (A*B)/(rA*rB);
  cvm::real const sin_phi = (C*B)/(rC*rB);

  cvm::rvector f1, f2, f3;

  rB = 1.0/rB;
  B *= rB;

  if (std::fabs(sin_phi) > 0.1) {
    rA = 1.0/rA;
    A *= rA;
    cvm::rvector const dcosdA = rA*(cos_phi*A-B);
    cvm::rvector const dcosdB = rB*(cos_phi*B-A);
    rA = 1.0;

    cvm::real const K = (1.0/sin_phi) * (180.0/PI);

	f1 = K * cvm::rvector::outer(r23, dcosdA);
	f3 = K * cvm::rvector::outer(dcosdB, r23);
	f2 = K * (cvm::rvector::outer(dcosdA, r12)
		   +  cvm::rvector::outer(r34, dcosdB));
  }
  else {
    rC = 1.0/rC;
    C *= rC;
    cvm::rvector const dsindC = rC*(sin_phi*C-B);
    cvm::rvector const dsindB = rB*(sin_phi*B-C);
    rC = 1.0;

    cvm::real    const K = (-1.0/cos_phi) * (180.0/PI);

    f1.x = K*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
              - r23.x*r23.y*dsindC.y
              - r23.x*r23.z*dsindC.z);
    f1.y = K*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
              - r23.y*r23.z*dsindC.z
              - r23.y*r23.x*dsindC.x);
    f1.z = K*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
              - r23.z*r23.x*dsindC.x
              - r23.z*r23.y*dsindC.y);

    f3 = cvm::rvector::outer(dsindB, r23);
    f3 *= K;

    f2.x = K*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
              +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
              +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
              +dsindB.z*r34.y - dsindB.y*r34.z);
    f2.y = K*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
              +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
              +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
              +dsindB.x*r34.z - dsindB.z*r34.x);
    f2.z = K*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
              +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
              +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
              +dsindB.y*r34.x - dsindB.x*r34.y);
  }

  F1 = -f1;
  F2 = -f2 + f1;
  F3 = -f3 + f2;
  F4 = f3;
}


template<bool calculate_gradients>
cvm::real colvar::angleswitch::switching_function(cvm::real const &ang_min,
                                                cvm::real const &ang_width,
                                                cvm::real const &power,
                                                cvm::real const &prefactor,
                                                cvm::atom &A1,
                                                cvm::atom &A2,
                                                cvm::atom &A3,
                                                cvm::atom &A4)
{
  cvm::real ang = this->ang_val(A1, A2, A3, A4);
  cvm::real angdist = ang - ang_min;
  cvm::real sign = 1.0;

  if (angdist > 180)  angdist -= 360;
  if (angdist < -180) angdist += 360;
  if (angdist < 0) sign = -1.0;

  //cvm::log(cvm::to_str(ang) + " " + cvm::to_str(ang_min) + " " + cvm::to_str(angdist) + "\n");

  angdist *= sign;

  if (angdist > ang_width) return cvm::real(1.0);

  cvm::real const xi = angdist/ang_width;

  // gradients
  if (calculate_gradients) {
    cvm::real const dFdang = prefactor*sign*pow(xi, power-1.0)/ang_width;
    cvm::rvector f1, f2, f3, f4;
    this->ang_val_grad(f1, f2, f3, f4);
    A1.grad += dFdang*f1;
    A2.grad += dFdang*f2;
    A3.grad += dFdang*f3;
    A4.grad += dFdang*f4;
  }

  return xi;
}

colvar::angleswitch::angleswitch(std::string const &conf)
  : cvc(conf)
{
  function_type = "angleswitch";
  x.type(colvarvalue::type_scalar);

  b_periodic = false;
  b_inverse_gradients = false;
  b_Jacobian_derivative = false;
  rebuildList = true;

  offset = 0.0;
  pairsum = 0.0;
  laststep = 0;
  progress = 0.0;

  parse_group(conf, "group1", group1);
  parse_group(conf, "group2", group2);
  parse_group(conf, "group3", group3);
  parse_group(conf, "group4", group4);
  atom_groups.push_back(&group1);
  atom_groups.push_back(&group2);
  atom_groups.push_back(&group3);
  atom_groups.push_back(&group4);

  get_keyval(conf, "offset", ang_offset,
              cvm::real(60));

  get_keyval(conf, "repeat", ang_period,
              cvm::real(120));

  get_keyval(conf, "width", ang_fluct,
              cvm::real(30));

  get_keyval(conf, "power", power,
              cvm::real(12.0));

  get_keyval(conf, "waitTime", nwait,
              int(100));

}

colvar::angleswitch::angleswitch()
{
  function_type = "angleswitch";
  x.type(colvarvalue::type_scalar);
  b_periodic = false;
  b_inverse_gradients = false;
  b_Jacobian_derivative = false;
  rebuildList = true;

  offset = 0.0;
  pairsum = 0.0;
  laststep = 0;
  progress = 0.0;
}

void colvar::angleswitch::build_bondlist()
{
  rebuildList = false;
  cvm::real ang;
  angList.clear();

  for (size_t i = 0; i < group1.size(); i++) {
    ang = this->ang_val(group1[i], group2[i], group3[i], group4[i]);
    ang = round((ang-ang_offset)/ang_period)*ang_period + ang_offset;
    if (ang > 180)  ang -= 360;
    if (ang < -180) ang += 360;
    angList.push_back(ang);
  }
}

void colvar::angleswitch::calc_value()
{
  // bookkeeping: check whether bond list rebuild is needed.
  // also, keep in mind that using '==' for floating point numbers
  // is a bad idea, so add a small tolerance for progress check
  // cases:
  // pairsum ~= 1: start/keep counting how long this takes (and eventually reset lists)
  // else, we increment laststep in order to have a moving frame of reference when the first case starts
  if (pairsum > 0.9999) {
    if ((cvm::step_absolute() - laststep) == nwait) {
      rebuildList = true;
      laststep = cvm::step_absolute();
      offset++;
    }
  } else {
    laststep = cvm::step_absolute();
  }

  // construct a new bondlist if requested
  if (rebuildList) build_bondlist();

  cvm::real pairterm;
  cvm::real dummy = 0.0;
  pairsum = 0.0;
  x.real_value = 0.0;

  // get all pairwise terms
  for (size_t i = 0; i < angList.size(); i++) {
    pairterm = this->switching_function<false> (angList[i], ang_fluct, power, dummy, group1[i], group2[i], group3[i], group4[i]);
    pairsum += pow(pairterm, power); 
  }

  // finish colvar calculation
  progress = pow(pairsum, 1.0/power);
  if (progress < 1.0) {
    progress = 0.5*(1.0-cos(PI*progress*progress));
  } else {
    progress = 1.0;
  }
  x.real_value = 2.0*cvm::real(offset) + progress;
}


void colvar::angleswitch::calc_gradients()
{
  // need to get full sum of pairwise terms before pairwise forces can be calculated
  calc_value();
  if (progress < 1.0) {
    cvm::real xi_t = pow(pairsum, 1.0/power);
    cvm::real const prefact = PI*xi_t*sin(PI*xi_t*xi_t)
                                *pow(pairsum, 1.0/power-1.0);

    for (size_t i = 0; i < angList.size(); i++) {
      this->switching_function<true> (angList[i], ang_fluct, power, prefact, group1[i], group2[i], group3[i], group4[i]);
    }
  }
}

void colvar::angleswitch::apply_force(colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force(force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force(force.real_value);

  if (!group3.noforce)
    group3.apply_colvar_force(force.real_value);

  if (!group4.noforce)
    group4.apply_colvar_force(force.real_value);
}


template<bool calculate_gradients>
cvm::real colvar::bondbreakmulti::switching_function(cvm::real const &rmin,
                                                     cvm::real const &rmax,
                                                     cvm::real const &power,
                                                     cvm::real const &prefactor,
                                                     cvm::atom &A1,
                                                     cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::real const dist = diff.norm();

  if (dist > rmax) return cvm::real(1.0);
  if (dist < rmin) return cvm::real(0.0);

  cvm::real const xi = (dist - rmin)/(rmax - rmin);

  // gradients
  if (calculate_gradients) {
    cvm::real const dFdr = prefactor*pow(xi, power-1.0)/(rmax-rmin);
    cvm::rvector const drdx = diff/dist;
    A1.grad += (-1.0)*dFdr*drdx;
    A2.grad +=        dFdr*drdx;
  }

  return xi;
}

colvar::bondbreakmulti::bondbreakmulti(std::string const &conf)
  : distance(conf)
{
  function_type = "bondbreakmulti";
  x.type(colvarvalue::type_scalar);

  // need to specify this explicitly because the distance() constructor
  // has set it to true
  b_inverse_gradients = false;
  rebuildList = true;
  offset = 0.0;
  pairsum = 0.0;
  laststep = 0;
  progress = 0.0;

  get_keyval(conf, "rmin11", rmin11,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmax11", rmax11,
              cvm::real(2.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rcut11", rcut11,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmin12", rmin12,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmax12", rmax12,
              cvm::real(2.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rcut12", rcut12,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmin22", rmin22,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rmax22", rmax22,
              cvm::real(2.0 * cvm::unit_angstrom()));

  get_keyval(conf, "rcut22", rcut22,
              cvm::real(1.0 * cvm::unit_angstrom()));

  get_keyval(conf, "power", power,
              cvm::real(12.0));

  get_keyval(conf, "waitTime", nwait,
              int(100));

}

colvar::bondbreakmulti::bondbreakmulti()
{
  function_type = "bondbreakmulti";
  x.type(colvarvalue::type_scalar);
}

void colvar::bondbreakmulti::build_bondlist()
{
  rebuildList = false;
  cvm::real const cut11 = rcut11*rcut11;
  cvm::real const cut22 = rcut22*rcut22;
  cvm::real const cut12 = rcut12*rcut12;
  cvm::real r2;

  pairList11_1.clear();
  pairList11_2.clear();
  pairList12_1.clear();
  pairList12_2.clear();
  pairList22_1.clear();
  pairList22_2.clear();

  for (size_t i = 0; i < group1.size() - 1; i++) {
    for (size_t j = i + 1; j < group1.size(); j++) {
      r2 = cvm::position_dist2(group1[i].pos, group1[j].pos);
      if (r2 < cut11) {
        pairList11_1.push_back(i);
        pairList11_2.push_back(j);
      }
    }
  }

  for (size_t i = 0; i < group2.size() - 1; i++) {
    for (size_t j = i + 1; j < group2.size(); j++) {
      r2 = cvm::position_dist2(group2[i].pos, group2[j].pos);
      if (r2 < cut22) {
        pairList22_1.push_back(i);
        pairList22_2.push_back(j);
      }
    }
  }

  for (size_t i = 0; i < group1.size(); i++) {
    for (size_t j = 0; j < group2.size(); j++) {
      r2 = cvm::position_dist2(group1[i].pos, group2[j].pos);
      if (r2 < cut12) {
        pairList12_1.push_back(i);
        pairList12_2.push_back(j);
      }
    }
  }

}

void colvar::bondbreakmulti::calc_value()
{
  // bookkeeping: check whether bond list rebuild is needed.
  // also, keep in mind that using '==' for floating point numbers
  // is a bad idea, so add a small tolerance for progress check
  // cases:
  // pairsum ~= 1: start/keep counting how long this takes (and eventually reset lists)
  // else, we increment laststep in order to have a moving frame of reference when the first case starts
  if (pairsum > 0.9999) {
    if ((cvm::step_absolute() - laststep) == nwait) {
      rebuildList = true;
      laststep = cvm::step_absolute();
      offset++;
    }
  } else {
    laststep = cvm::step_absolute();
  }

  // construct a new bondlist if requested
  if (rebuildList) build_bondlist();

  cvm::real pairterm;
  cvm::real dummy = 0.0;
  pairsum = 0.0;
  x.real_value = 0.0;

  // get all pairwise terms
  for (size_t i = 0; i < pairList11_1.size(); i++) {
    pairterm = switching_function<false> (rmin11, rmax11, power, dummy, group1[pairList11_1[i]], group1[pairList11_2[i]]);
    pairsum += pow(pairterm, power);
  }

  for (size_t i = 0; i < pairList22_1.size(); i++) {
    pairterm = switching_function<false> (rmin22, rmax22, power, dummy, group2[pairList22_1[i]], group2[pairList22_2[i]]);
    pairsum += pow(pairterm, power);
  }

  for (size_t i = 0; i < pairList12_1.size(); i++) {
    pairterm = switching_function<false> (rmin12, rmax12, power, dummy, group1[pairList12_1[i]], group2[pairList12_2[i]]);
    pairsum += pow(pairterm, power);
  }

  // finish colvar calculation
  progress = pow(pairsum, 1.0/power);
  if (progress < 1.0) {
    progress = 0.5*(1.0-cos(PI*progress*progress));
  } else {
    progress = 1.0;
  }
  x.real_value = 2.0*cvm::real(offset) + progress;
}


void colvar::bondbreakmulti::calc_gradients()
{
  // need to get full sum of pairwise terms before pairwise forces can be calculated
  calc_value();
  if (progress < 1.0) {
    if (pairsum < pow(0.001, power)) return;
    cvm::real xi_t = pow(pairsum, 1.0/power);
    cvm::real const prefact = PI*xi_t*sin(PI*xi_t*xi_t)
                                *pow(pairsum, 1.0/power-1.0);

    for (size_t i = 0; i < pairList11_1.size(); i++) {
      switching_function<true> (rmin11, rmax11, power, prefact, group1[pairList11_1[i]], group1[pairList11_2[i]]);
    }

    for (size_t i = 0; i < pairList22_1.size(); i++) {
      switching_function<true> (rmin22, rmax22, power, prefact, group2[pairList22_1[i]], group2[pairList22_2[i]]);
    }

    for (size_t i = 0; i < pairList12_1.size(); i++) {
      switching_function<true> (rmin12, rmax12, power, prefact, group1[pairList12_1[i]], group2[pairList12_2[i]]);
    }
  }
}

void colvar::bondbreakmulti::apply_force(colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force(force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force(force.real_value);
}