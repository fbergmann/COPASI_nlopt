// Copyright (C) 2019 by Pedro Mendes, Rector and Visitors of the
// University of Virginia, University of Heidelberg, and University
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2017 - 2018 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and University of
// of Connecticut School of Medicine.
// All rights reserved.

// Copyright (C) 2010 - 2016 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef CLAYOUTENGINE_H
#define CLAYOUTENGINE_H

#include <vector>
#include <sstream>

#include "copasi/odepack++/CLSODA.h"

#include <nlopt.hpp>

class CAbstractLayoutInterface;

class CLayoutEngine
{
public:
  struct Data
  {
    size_t dim;
    CLayoutEngine * pMethod;
  };

  CLayoutEngine(CAbstractLayoutInterface * l, bool so);

  double step();

  void calcForces(std::vector<double> & state, std::vector<double> & forces);
  void calcRHS(std::vector<double> & state, double* rhs);

  CAbstractLayoutInterface * getLayoutInterface()
  {return mpLayout;};

  // /tell the engine that a variable was changed outside
  //void changeVar(int index, double value);

  /**
   * request the current operation to stop
   */
  void requestStop();

  double stepNlopt();

  void setupNlopt(const std::string & algorithmName, const std::string & jsonOptions);
  
  void freeNlopt();

  
protected:
  double stepIntegration();

  double mInitialPot;

  int mStopAfterNthImprovement;
  
  int mDefaultStopAfterNthImprovement;
  
  nlopt::opt* mOpt;
  nlopt::opt* mLocal;


   static double nlopt_objective_function(unsigned n, const double * x, double * grad, void * data);



  CAbstractLayoutInterface * mpLayout;

  unsigned int mIndex;

  std::vector<double> mVariables;
  std::vector<double> mRhs;

  //for RK4
  std::vector<double> mRhsA;
  std::vector<double> mRhsB;
  std::vector<double> mRhsC;
  std::vector<double> mVar2;

  bool mSecondOrder;
//  std::vector<double> mVelocities;
//  std::vector<double> mForce;

  CLSODA mLSODA;

  /**
   * LSODA C_FLOAT64 work area
   */
  std::vector< C_FLOAT64 > mDWork;

  /**
   * LSODA C_INT work area
   */
  std::vector< C_INT > mIWork;

  /**
   * Stream to capture LSODA error messages
   */
  std::ostringstream mErrorMsg;

  /**
   *  LSODA state.
   */
  C_INT mLsodaStatus;

  /**
   * mData.dim is the dimension of the ODE system.
   * mData.pMethod contains CLsodaMethod * this to be used in the static method EvalF
   */
  Data mData;

  /**
   *  Current time.
   */
  C_FLOAT64 mTime;

  static void EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

  /**
   *  This evaluates the derivatives
   */
  void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

  bool mStopRequested;

};

#endif // CLAYOUTENGINE_H
