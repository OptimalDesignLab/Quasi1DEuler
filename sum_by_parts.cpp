/**
 * \file sum_by_parts.cpp
 * \brief member defintions for SumByParts and derivied classes
 * \author  Jason Hicken <jason.hicken@gmail.com>
 * \version 1.0
 */

#include "./sum_by_parts.hpp"

#include <math.h>

#include <ostream>
#include <iostream>

#include "./inner_prod_vector.hpp"

using std::cout;
using std::cerr;
using std::endl;

// ======================================================================

SumByParts::SumByParts(int nodes, int order) {
  Define(nodes, order);
}

// ======================================================================

void SumByParts::Define(int nodes, int order) {
  if ((order <= 1) || (order > 4)) {
    cerr << "SumByParts::Define(): invalid value for order.";
    throw(-1);
  }
  order_ = order;
  if ( ((order_ == 2) && (nodes < 3)) ||
       ((order_ == 3) && (nodes < 9)) ||
       ((order_ == 4) && (nodes < 13)) ) {
    cerr << "SumByParts::Define(): "
         << "invalid number of nodes for order.";
    throw(-1);
  }
  num_nodes_ = nodes;
  ibeg_.resize(num_nodes_);
  iend_.resize(num_nodes_);

  if (order_ == 2) {
    numh_ = 1;
    hinv_.resize(numh_);
    hinv_(0) = 2.0;
  } else if (order_ == 3) {
    numh_ = 4;
    hinv_.resize(numh_);
    hinv_(0) = 48.0/17.0;
    hinv_(1) = 48.0/59.0;
    hinv_(2) = 48.0/43.0;
    hinv_(3) = 48.0/49.0;
  } else if (order_ == 4) {
    numh_ = 6;
    hinv_.resize(numh_);
    hinv_(0) = 43200.0/13649.0;
    hinv_(1) = 8640.0/12013.0;
    hinv_(2) = 4320.0/2711.0;
    hinv_(3) = 4320.0/5359.0;
    hinv_(4) = 8640.0/7877.0;
    hinv_(5) = 43200.0/43801.0;
  }
}

// ======================================================================

double SumByParts::InnerProductSBP(
    const int & num_var, const InnerProdVector & u,
    const InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) ||
       (v.size() != num_var*num_nodes_) ) {
    cerr << "SumByParts::InnerProductSBP(): "
         << "inconsistent sizes -> u.size(), v.size() = "
         << u.size() << " " << v.size() << endl;;
    throw(-1);
  }
  double prod = 0.0;
  for (int i = numh_; i < num_nodes_-numh_; i++) {
    int iptr = i*num_var;
    for (int n = 0; n < num_var; n++)
      prod += u(iptr+n)*v(iptr+n);
  }
  for (int i = 0; i < numh_; i++) {
    int ptr0 = i*num_var;
    int ptr1 = (num_nodes_-i-1)*num_var;
    for (int n = 0; n < num_var; n++) {
      prod += (u(ptr0+n)*v(ptr0+n))/hinv_(i);
      prod += (u(ptr1+n)*v(ptr1+n))/hinv_(i);
    }
  }
  //prod /= static_cast<double>(num_nodes_-1);
  return prod;
}

// ======================================================================

void SumByParts::HTimesVector(const int & num_var, InnerProdVector & u,
                              InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) ||
       (v.size() != num_var*num_nodes_) ) {
    cerr << "SumByParts::HTimesVector(): "
         << "inconsistent sizes -> u.size(), v.size() = "
         << u.size() << " " << v.size() << endl;;
    throw(-1);
  }
  for (int i = 0; i < numh_; i++) {
    int ptr0 = i*num_var;
    int ptr1 = (num_nodes_-i-1)*num_var;
    for (int n = 0; n < num_var; n++) {
      v(ptr0+n) = u(ptr0+n)/hinv_(i);
      v(ptr1+n) = u(ptr1+n)/hinv_(i);
    }
  }
}

// ======================================================================

void SumByParts::HinvTimesVector(const int & num_var,
                                 InnerProdVector & u,
                                 InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) ||
       (v.size() != num_var*num_nodes_) ) {
    cerr << "SumByParts::HTimesVector(): "
         << "inconsistent sizes -> u.size(), v.size() = "
         << u.size() << " " << v.size() << endl;;
    throw(-1);
  }
  for (int i = 0; i < numh_; i++) {
    int ptr0 = i*num_var;
    int ptr1 = (num_nodes_-i-1)*num_var;
    for (int n = 0; n < num_var; n++) {
      v(ptr0+n) = u(ptr0+n)*hinv_(i);
      v(ptr1+n) = u(ptr1+n)*hinv_(i);
    }
  }
}

// ======================================================================

SBP1stDerivative::SBP1stDerivative(int nodes, int order) : 
    SumByParts(nodes, order) {
  Define(nodes, order);
}

// ======================================================================

void SBP1stDerivative::Define(int nodes, int order) {
  SumByParts::Define(nodes, order);
  if (order_ == 2) {
    ja_.resize(6);
    as_.resize(6);
      
    // first row
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0; as_(ptr) = -1.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  1.0;
    iend_(0) = ptr;
    
    // interior rows
    ptr = ptr + 1;
    for (int i = 1; i < num_nodes_-1; i++)  ibeg_(i) = ptr;
    ja_(ptr) = -1; as_(ptr) = -0.5;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  0.5;
    for (int i = 1; i < num_nodes_-1; i++) iend_(i) = ptr;
    
    // last row
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    ja_(ptr) = -1; as_(ptr) = -1.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  1.0;
    iend_(num_nodes_-1) = ptr;

  } else if (order_ == 3) {
    ja_.resize(32);
    as_.resize(32);
      
    // first row
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0; as_(ptr) = -24.0/17.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  59.0/34.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -4.0/17.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) = -3.0/34.0;
    iend_(0) = ptr;

    // second row
    ptr = ptr + 1;
    ibeg_(1) = ptr;
    ja_(ptr) = -1; as_(ptr) = -0.5;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  0.5;
    iend_(1) = ptr;

    // third row
    ptr = ptr + 1;
    ibeg_(2) = ptr;
    ja_(ptr) = -2; as_(ptr) =  4.0/43.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -59.0/86.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  59.0/86.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -4.0/43.0;
    iend_(2) = ptr;

    // fourth row
    ptr = ptr + 1;
    ibeg_(3) = ptr;
    ja_(ptr) = -3; as_(ptr) =  3.0/98.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -59.0/98.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  32.0/49.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -4.0/49.0;
    iend_(3) = ptr;

    // interior rows
    ptr = ptr + 1;
    for (int i = 4; i < num_nodes_-4; i++) ibeg_(i) = ptr;
    ja_(ptr) = -2; as_(ptr) =  1.0/12.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -2.0/3.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  2.0/3.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -1.0/12.0;
    for (int i = 4; i < num_nodes_-4; i++) iend_(i) = ptr;

    // (N-3)rd row
    ptr = ptr + 1;
    ibeg_(num_nodes_-4) = ptr;
    ja_(ptr) = -2; as_(ptr) =  4.0/49.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -32.0/49.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  59.0/98.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) = -3.0/98.0;
    iend_(num_nodes_-4) = ptr;

    // (N-2)nd row
    ptr = ptr + 1;
    ibeg_(num_nodes_-3) = ptr;
    ja_(ptr) = -2; as_(ptr) =  4.0/43.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -59.0/86.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  59.0/86.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -4.0/43.0;
    iend_(num_nodes_-3) = ptr;

    // (N-1)st row
    ptr = ptr + 1;
    ibeg_(num_nodes_-2) = ptr;
    ja_(ptr) = -1; as_(ptr) = -0.5;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  0.5;
    iend_(num_nodes_-2) = ptr;

    // Nth row
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    ja_(ptr) = -3; as_(ptr) =  3.0/34.0;
    ptr = ptr + 1;
    ja_(ptr) = -2; as_(ptr) =  4.0/17.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -59.0/34.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  24.0/17.0;
    iend_(num_nodes_-1) = ptr;

  } else if (order_ == 4) {
    ja_.resize(80);
    as_.resize(80);

    // first row
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0;
    as_(ptr) = -1.582533518939116418785258993332844897062;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  2.033426786468126253898161347360808173712;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.1417052898146741610733887894481170575600;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) = -0.4501096599735708523162117824920488989702;
    ptr = ptr + 1;
    ja_(ptr) =  4;
    as_(ptr) =  0.1042956382142412661862395105494407610836;
    ptr = ptr + 1;
    ja_(ptr) =  5;
    as_(ptr) =  0.03662604404499391209045870736276191879693;
    iend_(0) = ptr;

    // second row
    ptr = ptr + 1;
    ibeg_(1) = ptr;
    ja_(ptr) = -1;
    as_(ptr) = -0.4620701275035953590186631853846278325646;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.2873679417026202568532985205129449923126;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) =  0.2585974499280928196267362923074433487080;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) = -0.06894808744606961472005221923058251153103;
    ptr = ptr + 1;
    ja_(ptr) =  4;
    as_(ptr) = -0.01494717668104810274131940820517799692506;
    iend_(1) = ptr;

    // third row
    ptr = ptr + 1;
    ibeg_(2) = ptr;
    ja_(ptr) = -2;
    as_(ptr) =  0.07134398748360337973038301686379010397038;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.6366933020423417826592908754928085932593;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.6067199374180168986519150843189505198519;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.02338660408468356531858175098561718651857;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) = -0.01798401877459493040442547470431484404443;
    iend_(2) = ptr;
      
    // fourth row
    ptr = ptr + 1;
    ibeg_(3) = ptr;
    ja_(ptr) = -3;
    as_(ptr) =  0.1146397975178068401430112823144985150596;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) = -0.2898424301162697370942324201800071793273;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.3069262456316931913128086944558079603132;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.5203848121857539166740071338174418292578;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.05169127637022742348368508279860701098408;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) =  0.01343534241462959507370778130248180630715;
    iend_(3) = ptr;

    // fifth row
    ptr = ptr + 1;
    ibeg_(4) = ptr;
    ja_(ptr) = -4;
    as_(ptr) = -0.03614399304268576976452921364705641609825;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) =  0.1051508663818248421520867474440761344449;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.01609777419666805778308369351834662756172;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.7080721616106272031118456849378369336023;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.7692160858661111736140494493705980473867;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.1645296432652024882569506157166433921544;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) =  0.01828107147391138758410562396851593246160;
    iend_(4) = ptr;

    // sixth row
    ptr = ptr + 1;
    ibeg_(5) = ptr;
    ja_(ptr) = -5;
    as_(ptr) = -0.01141318406360863692889821914555232596651;
    ptr = ptr + 1;
    ja_(ptr) = -4;
    as_(ptr) =  0.02049729840293952857599941220163960606616;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) =  0.01113095018331244864875173213474522093204;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.06324365883611076515355091406993789453750;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.6916640154753724474963890679085181638850;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.7397091390607520376247117645715851236273;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.1479418278121504075249423529143170247255;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) =  0.01643798086801671194721581699047966941394;
    iend_(5) = ptr;

    // interior rows
    ptr = ptr + 1;
    for (int i = 6; i < num_nodes_-6; i++) ibeg_(i) = ptr;
    ja_(ptr) = -3;
    as_(ptr) = -1.0/60.0;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  3.0/20.0;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -3.0/4.0;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  3.0/4.0;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -3.0/20.0;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) =  1.0/60.0;
    for (int i = 6; i < num_nodes_-6; i++) iend_(i) = ptr;

    // _(N-5)th row
    ptr = ptr + 1;
    ibeg_(num_nodes_-6) = ptr;
    ja_(ptr) =  5;
    as_(ptr) =  0.01141318406360863692889821914555232596651;
    ptr = ptr + 1;
    ja_(ptr) =  4;
    as_(ptr) = -0.02049729840293952857599941220163960606616;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) = -0.01113095018331244864875173213474522093204;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.06324365883611076515355091406993789453750;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.6916640154753724474963890679085181638850;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.7397091390607520376247117645715851236273;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.1479418278121504075249423529143170247255;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) = -0.01643798086801671194721581699047966941394;
    iend_(num_nodes_-6) = ptr;
      
    // _(N-4)th row
    ptr = ptr + 1;
    ibeg_(num_nodes_-5) = ptr;
    ja_(ptr) =  4;
    as_(ptr) =  0.03614399304268576976452921364705641609825;
    ptr = ptr + 1;
    ja_(ptr) =  3;
    as_(ptr) = -0.1051508663818248421520867474440761344449;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) = -0.01609777419666805778308369351834662756172;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.7080721616106272031118456849378369336023;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.7692160858661111736140494493705980473867;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.1645296432652024882569506157166433921544;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) = -0.01828107147391138758410562396851593246160;
    iend_(num_nodes_-5) = ptr;

    // _(N-3)rd row
    ptr = ptr + 1;
    ibeg_(num_nodes_-4) = ptr;
    ja_(ptr) =  3;
    as_(ptr) = -0.1146397975178068401430112823144985150596;
    ptr = ptr + 1;
    ja_(ptr) =  2;
    as_(ptr) =  0.2898424301162697370942324201800071793273;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.3069262456316931913128086944558079603132;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.5203848121857539166740071338174418292578;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.05169127637022742348368508279860701098408;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) = -0.01343534241462959507370778130248180630715;
    iend_(num_nodes_-4) = ptr;
      
    // _(N-2)nd row
    ptr = ptr + 1;
    ibeg_(num_nodes_-3) = ptr;
    ja_(ptr) =  2;
    as_(ptr) = -0.07134398748360337973038301686379010397038;
    ptr = ptr + 1;
    ja_(ptr) =  1;
    as_(ptr) =  0.6366933020423417826592908754928085932593;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.6067199374180168986519150843189505198519;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.02338660408468356531858175098561718651857;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) =  0.01798401877459493040442547470431484404443;
    iend_(num_nodes_-3) = ptr;

    // _(N-1)st row
    ptr = ptr + 1;
    ibeg_(num_nodes_-2) = ptr;
    ja_(ptr) =  1;
    as_(ptr) =  0.4620701275035953590186631853846278325646;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -0.2873679417026202568532985205129449923126;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) = -0.2585974499280928196267362923074433487080;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) =  0.06894808744606961472005221923058251153103;
    ptr = ptr + 1;
    ja_(ptr) = -4;
    as_(ptr) =  0.01494717668104810274131940820517799692506;
    iend_(num_nodes_-2) = ptr;
            
    // Nth row
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    ja_(ptr) =  0;
    as_(ptr) =  1.582533518939116418785258993332844897062;
    ptr = ptr + 1;
    ja_(ptr) = -1;
    as_(ptr) = -2.033426786468126253898161347360808173712;
    ptr = ptr + 1;
    ja_(ptr) = -2;
    as_(ptr) =  0.1417052898146741610733887894481170575600;
    ptr = ptr + 1;
    ja_(ptr) = -3;
    as_(ptr) =  0.4501096599735708523162117824920488989702;
    ptr = ptr + 1;
    ja_(ptr) = -4;
    as_(ptr) = -0.1042956382142412661862395105494407610836;
    ptr = ptr + 1;
    ja_(ptr) = -5;
    as_(ptr) = -0.03662604404499391209045870736276191879693;
    iend_(num_nodes_-1) = ptr;
  }
}

// ======================================================================

void SBP1stDerivative::Apply(const int & num_var, 
                             const InnerProdVector & u,
                             InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) || 
       (v.size() != num_var*num_nodes_) ) {
    cerr << "SBP1stDerivative::Apply(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  // apply the SBP operator
  v = 0.0;
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var; n++)
        v(iptr + n) += coeff*u(jptr + n);
    }
  }
}

// ======================================================================

void SBP1stDerivative::ApplyTranspose(const int & num_var,
                                      const InnerProdVector & u,
                                      InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) || 
       (v.size() != num_var*num_nodes_) ) {
    cerr << "SBP1stDerivative::ApplyTranspose(): "
         << "inconsistent sizes.";
    throw(-1);
  }
  // apply the transposed SBP operator D^{T} = Q^{T} H^{-1}
  v = 0.0;
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var; n++)
        v(jptr + n) += coeff*u(iptr + n);
    }
  }
}

// ======================================================================

SBPDissipation::SBPDissipation(int nodes, int norm_order, int order) : 
    SumByParts(nodes, norm_order) {
  Define(nodes, norm_order, order);
}

// ======================================================================

void SBPDissipation::Define(int nodes, int norm_order, int order) {
  SumByParts::Define(nodes, norm_order);
  order_ = order; // in case norm_order is different from order
  if ( ((order_ == 1) && (nodes < 1)) ||
       ((order_ == 2) && (nodes < 3)) ||
       ((order_ == 3) && (nodes < 5)) ) {
    cerr << "SBPDissipation::Define(): "
         << "invalid number of nodes for order.";
    throw(-1);
  }

  if (order_ == 1) {    
    ja_.resize(2);
    as_.resize(2);

    // interior, half-nodes
    int ptr = 0;
    for (int i = 0; i < num_nodes_-1; i++) ibeg_(i) = ptr;
    ja_(ptr) =  0; as_(ptr) = -1.0;
    ptr = ptr + 1;
    for (int i = 0; i < num_nodes_-1; i++) iend_(i) = ptr;
    ja_(ptr) =  1; as_(ptr) =  1.0;
    
    // last node; by using iend < ibeg, nothing is added to this node
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    iend_(num_nodes_-1) = ptr - 1;

  } else if (order_ == 2) {
    ja_.resize(11);
    as_.resize(11);

    // first row
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0; as_(ptr) =  2.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) = -5.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) =  4.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) = -1.0;
    iend_(0) = ptr;

    // interior rows
    ptr = ptr + 1;
    for (int i = 1; i < num_nodes_-1; i++) ibeg_(i) = ptr;
    ja_(ptr) = -1; as_(ptr) =  1.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) = -2.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  1.0;
    for (int i = 1; i < num_nodes_-1; i++) iend_(i) = ptr;

    // last row
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    ja_(ptr) = -3; as_(ptr) = -1.0;
    ptr = ptr + 1;
    ja_(ptr) = -2; as_(ptr) =  4.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -5.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  2.0;
    iend_(num_nodes_-1) = ptr;
    
  } else if (order_ == 3) {
    ja_.resize(14);
    as_.resize(14);

    // first half node
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0; as_(ptr) = -1.0/18.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  7.0/36.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -1.0/4.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) =  5.0/36.0;
    ptr = ptr + 1;
    ja_(ptr) =  4; as_(ptr) = -1.0/36.0;
    iend_(0) = ptr;

    // interior, half-nodes
    ptr = ptr + 1;
    for (int i = 1; i < num_nodes_-1; i++) ibeg_(i) = ptr;
    ja_(ptr) = -1; as_(ptr) = -1.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  3.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) = -3.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) =  1.0;
    for (int i = 1; i < num_nodes_-1; i++) iend_(i) = ptr;

    // last half node
    ptr = ptr + 1;
    ibeg_(num_nodes_-2) = ptr;
    ja_(ptr) = -3; as_(ptr) = -1.0/36.0;
    ptr = ptr + 1;
    ja_(ptr) = -2; as_(ptr) =  5.0/36.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -1.0/4.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  7.0/36.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) = -1.0/18.0;
    iend_(num_nodes_-2) = ptr;

    // last node; by using iend < ibeg, nothing is added to this node
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    iend_(num_nodes_-1) = ptr - 1;

  } else if (order_ == 4) {
    ja_.resize(25);
    as_.resize(25);

    // first node
    int ptr = 0;
    ibeg_(0) = ptr;
    ja_(ptr) =  0; as_(ptr) =  1.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) =  6.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  4; as_(ptr) =  1.0;
    iend_(0) = ptr;

    // second node
    ptr = ptr + 1;
    ibeg_(1) = ptr;
    ja_(ptr) = -1; as_(ptr) =  1.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  6.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  3; as_(ptr) =  1.0;
    iend_(1) = ptr;

    // interior nodes
    ptr = ptr + 1;
    for (int i = 2; i < num_nodes_-2; i++) ibeg_(i) = ptr;
    ja_(ptr) = -2; as_(ptr) =  1.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  6.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  2; as_(ptr) =  1.0;
    for (int i = 2; i < num_nodes_-2; i++) iend_(i) = ptr;

    // second last node
    ptr = ptr + 1;
    ibeg_(num_nodes_-2) = ptr;
    ja_(ptr) = -3; as_(ptr) =  1.0;
    ptr = ptr + 1      ;
    ja_(ptr) = -2; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) =  6.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  1; as_(ptr) =  1.0;
    iend_(num_nodes_-2) = ptr;
      
    // last node
    ptr = ptr + 1;
    ibeg_(num_nodes_-1) = ptr;
    ja_(ptr) = -4; as_(ptr) =  1.0;
    ptr = ptr + 1      ;
    ja_(ptr) = -3; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) = -2; as_(ptr) =  6.0;
    ptr = ptr + 1;
    ja_(ptr) = -1; as_(ptr) = -4.0;
    ptr = ptr + 1;
    ja_(ptr) =  0; as_(ptr) =  1.0;
    iend_(num_nodes_-1) = ptr;
  }
}

// ======================================================================

void SBPDissipation::Apply(const int & num_var, 
                           const InnerProdVector & u,
                           const InnerProdVector & b,
                           InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) || 
       (v.size() != num_var*num_nodes_) ||
       (b.size() != num_nodes_) ) {
    cerr << "SBPDissipation::Apply(): "
         << "inconsistent sizes.";
    throw(-1);
  }

  InnerProdVector dq(num_var*num_nodes_, 0.0);
  InnerProdVector tmp(num_var, 0.0);

  // set dq = diag(b) * D_p * u and v = 0.0
  int bit = 0;
  if (ibeg_(num_nodes_-1) > iend_(num_nodes_-1)) bit = 1;
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    tmp = 0.0;
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var; n++)
        tmp(n) += coeff*u(jptr+n);
    }
    
    int i_plus_bit = std::min(i + bit, num_nodes_-1);
    double diag = 0.5*( b(i) + b(i_plus_bit) );
    for (int n = 0; n < num_var; n++) {
      dq(iptr+n) = diag * tmp(n);
      v(iptr+n) = 0.0;
    }
  }
  // apply D_p^T to the values in dq  
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    for (int n = 0; n < num_var; n++)
      tmp(n) = dq(iptr+n);
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var;  n++) 
        v(jptr+n) += coeff*tmp(n);
    }
  }
  // apply the inverse diagonal norm
  for (int j = 0; j < numh_; j++) {
    int ptr0 = j*num_var;
    int ptr1 = (num_nodes_-j-1)*num_var;
    for (int n = 0; n < num_var; n++) {
      v(ptr0+n) *= hinv_(j);
      v(ptr1+n) *= hinv_(j);
    }
  }
}

// ======================================================================

void SBPDissipation::ApplyTranspose(const int & num_var, 
                           const InnerProdVector & u,
                           const InnerProdVector & b,
                           InnerProdVector & v) const {
  // check for consistent sizes
  if ( (u.size() != num_var*num_nodes_) || 
       (v.size() != num_var*num_nodes_) ||
       (b.size() != num_nodes_) ) {
    cerr << "SBPDissipation::Apply(): "
         << "inconsistent sizes.";
    throw(-1);
  }

  InnerProdVector dq(num_var*num_nodes_, 0.0);
  InnerProdVector tmp(num_var, 0.0);

  // apply the inverse diagonal norm
  v = u;
  for (int j = 0; j < numh_; j++) {
    int ptr0 = j*num_var;
    int ptr1 = (num_nodes_-j-1)*num_var;
    for (int n = 0; n < num_var; n++) {
      v(ptr0+n) *= hinv_(j);
      v(ptr1+n) *= hinv_(j);
    }
  }

  // set dq = diag(b) * D_p * H^{-1} u
  int bit = 0;
  if (ibeg_(num_nodes_-1) > iend_(num_nodes_-1)) bit = 1;
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    tmp = 0.0;
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var; n++)
        tmp(n) += coeff*v(jptr+n);
    }
    
    int i_plus_bit = std::min(i + bit, num_nodes_-1);
    double diag = 0.5*( b(i) + b(i_plus_bit) );
    for (int n = 0; n < num_var; n++) {
      dq(iptr+n) = diag * tmp(n);
    }
  }
  v = 0.0;

  // apply D_p^T to the values in dq  
  for (int i = 0; i < num_nodes_; i++) {
    int iptr = i*num_var;
    for (int n = 0; n < num_var; n++)
      tmp(n) = dq(iptr+n);
    for (int j = ibeg_(i); j <= iend_(i); j++) {
      int jptr = (i + ja_(j))*num_var;
      double coeff = as_(j);
      for (int n = 0; n < num_var;  n++) 
        v(jptr+n) += coeff*tmp(n);
    }
  }
}
