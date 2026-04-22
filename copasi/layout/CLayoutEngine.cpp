// Copyright (C) 2019 - 2025 by Pedro Mendes, Rector and Visitors of the
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include "copasi/copasi.h"
#include "CLayoutEngine.h"
#include "CAbstractLayoutInterface.h"

#include <copasi/utilities/json.hpp>


CLayoutEngine::CLayoutEngine(CAbstractLayoutInterface * l, bool /* so */)
  : mpLayout(l)
  , mSecondOrder(false)
  , mStopRequested(false)
{
  if (!mpLayout) return;

  mVariables = mpLayout->getInitialValues();

  if (mSecondOrder) mVariables.resize(mVariables.size() * 2, 0);

  mRhs.resize(mVariables.size(), 0);

  //for RK4
  mRhsA.resize(mVariables.size(), 0);
  mRhsB.resize(mVariables.size(), 0);
  mRhsC.resize(mVariables.size(), 0);
  mVar2.resize(mVariables.size(), 0);

  //init LSODA
  mLsodaStatus = 1;
  //mState = 1;
  //mJType = 2;
  mErrorMsg.str("");

  mData.pMethod = this;
  mData.dim = mVariables.size();
  //mAtol = mpModel->initializeAtolVector(*mpAbsoluteTolerance, *mpReducedModel);
  //mY = mpState->beginIndependent();

  //mRtol = *mpRelativeTolerance;

  mTime = 0;

  mDWork.resize(22 + mData.dim * std::max<size_t>(16, mData.dim + 9) + 3 * 0);
  mDWork[4] = mDWork[5] = mDWork[6] = mDWork[7] = mDWork[8] = mDWork[9] = 0.0;
  mDWork[7] = 0; //minstep
  mIWork.resize(20 + mData.dim);
  mIWork[4] = mIWork[6] = mIWork[9] = 0;

  mIWork[5] = 100; //*mpMaxInternalSteps;
  mIWork[7] = 12;
  mIWork[8] = 5;
  mLSODA.setOstream(mErrorMsg);
}

void CLayoutEngine::calcRHS(std::vector<double> & state, double* rhs)
{
  std::vector<double> forces; forces.resize(mpLayout->getNumVariables());
  calcForces(state, forces);

  size_t i, imax = mpLayout->getNumVariables();

  for (i = 0; i < imax; ++i)
    {
      if (mStopRequested)
        break;

      if (mSecondOrder)
        {
          rhs[i + imax] = (forces[i] * 0.04 - state[i + imax] * 0.05) / mpLayout->getMassVector()[i];
          rhs[i] = state[i + imax];
        }
      else
        {
          //first order
          rhs[i] = forces[i] / mpLayout->getMassVector()[i];
        }
    }
}

void CLayoutEngine::calcForces(std::vector<double> & state, std::vector<double> & forces)
{
  if (!mpLayout) return;

  //do numerical derivations...
  mpLayout->setState(state);

  //assymmetric algorithm
  /*double pot0 = mpLayout->getPotential();

  double store, newpot;
  unsigned int i, imax = mpLayout->getNumVariables();
  for (i=0; i<imax; ++i)
  {
    store = mVariables[i];
    mVariables[i]+=1.0;
    mpLayout->setState(mVariables);
    newpot = mpLayout->getPotential();
    mVariables[i] = store;
    mForce[i] = pot0-newpot;
  }*/

  //symmetric algorithm
  double pot0;
  double store, newpot;
  size_t i, imax = mpLayout->getNumVariables();

  for (i = 0; i < imax; ++i)
    {
      if (mStopRequested)
        break;

      store = state[i];
      //std::cout << "var " << store;
      state[i] -= 0.5;
      mpLayout->setState(state);
      pot0 = mpLayout->getPotential();
      //std::cout << "pot0 " << pot0 << std::endl;
      state[i] += 1.0;
      mpLayout->setState(state);
      newpot = mpLayout->getPotential();
      state[i] = store;
      forces[i] = pot0 - newpot;
      //std::cout << "force " << i << " = " << forces[i] << std::endl;
    }
}

// include msvc debug header to debug print 
//#include <windows.h>
//
//void debugPrint(const std::string& message)
//{
//  OutputDebugStringA(message.c_str());
//}

//#define debugPrint(msg) {std::cout << msg;}
#define debugPrint(msg) {}


#include <nlopt.hpp>
double CLayoutEngine::nlopt_objective_function(unsigned n, const double * x, double * grad, void * data)
{
  
  // Cast the data pointer back
  CLayoutEngine * pThis = static_cast< CLayoutEngine * >(data);
  
  if (pThis->mStopRequested || !pThis->mpLayout)
    {
      throw nlopt::forced_stop();      
    }


  auto vars = std::vector< double >(x, x + n);
  auto backup = pThis->mVariables;
  

  pThis->mVariables = vars;

  if (grad)
  {
    std::vector<double> forces; forces.resize(n);
    pThis->calcForces(vars, forces);
    for (size_t i = 0; i < n; ++i)
    {
        grad[i] = -forces[i];
    }
  }

  
  pThis->mpLayout->setState(pThis->mVariables);
  double value = pThis->mpLayout->getPotential();

  if (value < pThis->mInitialPot)
  {
      if (pThis->mStopAfterNthImprovement > 0)
      {
          pThis->mStopAfterNthImprovement--;
      }

      if (pThis->mStopAfterNthImprovement == 0)
      {
          pThis->mStopRequested = true;
      }
  }
  
  return value;
}

double CLayoutEngine::stepNlopt(const std::string& algorithmName, const std::string & jsonConfig)
{
  if (!mpLayout) return -1.0;
  mStopRequested = false;
  size_t i, imax = mVariables.size();
  
// parse json config
  nlohmann::json config = nlohmann::json::parse(jsonConfig, nullptr, false, true);

  int algorithm = nlopt_algorithm_from_string(algorithmName.c_str());
  // if json contains algorithm, parse from name
  if (config.contains("algorithm")) {
      if (config["algorithm"].is_string())
        algorithm = nlopt_algorithm_from_string(config["algorithm"].get< std::string >().c_str());
      else
        algorithm = config["algorithm"].get<int>();
  }

  int maxeval;
  if (config.contains("maxeval")) {
    maxeval = config["maxeval"].get<int>();
  }
  else {
    maxeval = 20000;
  }

  if (config.contains("stopAfterNthImprovement"))
    {
      mStopAfterNthImprovement = config["stopAfterNthImprovement"].get< int >();
    }
  else
  {
      mStopAfterNthImprovement = -1;
  }

  // store state
  std::vector<double> initialState(mVariables);

  // evaluate initial potential
  mpLayout->setState(mVariables);
  mInitialPot = mpLayout->getPotential();
  debugPrint("initial potential: " + std::to_string(mInitialPot) + "\n");
  // setup bounds
  double lowerBound = 0;
  double upperBound = 10000;
  if (config.contains("lowerBound") && config["lowerBound"].get< double >() != -1.0)
    {
    lowerBound = config["lowerBound"].get<double>();
  }
  if (config.contains("upperBound") && config["lowerBound"].get< double >() != -1.0)
    {
    upperBound = config["upperBound"].get<double>();
  }
  std::vector<double> lowerBounds(imax, lowerBound);
  std::vector<double> upperBounds(imax, upperBound);

  double multiplier = 1; 
  if (config.contains("multiplier")) {
    multiplier = config["multiplier"].get<double>();
  }
  if (multiplier != 0)
    {
      for (i = 0; i < imax; ++i)
        {
          lowerBounds[i] = mVariables[i] - mVariables[i] * multiplier;
          upperBounds[i] = mVariables[i] + mVariables[i] * multiplier;
        }
    }

  // optimize 
  std::vector<double> x(mVariables);

  nlopt::opt opt((nlopt::algorithm)algorithm, (unsigned)imax);
  debugPrint(std::string("algorithm: ") + opt.get_algorithm_name());
  //opt.set_xtol_rel(1e-4);
  //opt.set_ftol_rel(1e-4);
  opt.set_maxeval(maxeval);
  opt.set_lower_bounds(lowerBounds);
  opt.set_upper_bounds(upperBounds);

  if (config.contains("initialStep") && config["initialStep"].get<double>() != 0)
    opt.set_initial_step(config["initialStep"].get< double >());



  opt.set_min_objective(nlopt_objective_function, this);

  double minf = -1.0;
  nlopt::result result = (nlopt::result) - 1;
  try
  {
      result = opt.optimize(x, minf);
    }
  catch (std::exception & e)
    {
      debugPrint("NLopt failed: " + std::string(e.what()) + "\n");      
    }
  
  

  if (result < 0 || minf > mInitialPot) {
      debugPrint("NLopt optimization failed or did not improve potential. Result: " + std::to_string(result) + " minf: " + std::to_string(minf) + "\n");
    mVariables = initialState;
      return mInitialPot;
  }

  debugPrint("NLopt optimization succeeded. Result: " + std::to_string(result) + " minf: " + std::to_string(minf) + "\n");

  for (i = 0; i < imax; ++i)
    {
      mVariables[i] = x[i];
    }

  return minf;

}

double CLayoutEngine::step()
{
  if (!mpLayout) return -1.0;

  mStopRequested = false;

  size_t i, imax = mVariables.size();

  //Euler

  mpLayout->setState(mVariables);
  double pot = mpLayout->getPotential();

  calcRHS(mVariables, &mRhs[0]);

  double dt = 3;

  for (i = 0; i < imax; ++i)
    {
      mVariables[i] += mRhs[i] * dt;
    }

  double newpot;

  for (;;)
    {
      mpLayout->setState(mVariables);
      newpot = mpLayout->getPotential();

      if (mStopRequested) break;

      if (newpot < pot) break;

      if (dt < 1e-3) break;

      //step back
      dt *= 0.5;

      for (i = 0; i < imax; ++i)
        {
          mVariables[i] -= mRhs[i] * dt;
        }
    }

  //std::cout << dt << "   " << newpot << std::endl;
  return newpot;
}

double CLayoutEngine::stepIntegration()
{
  if (!mpLayout)
    return -1.0;

  const double dt = 0.2;

  size_t i, imax = mVariables.size();

  //Euler
  /*calcRHS(mVariables, &mRhs[0]);

  for (i = 0; i < imax; ++i)
    {
      mVariables[i] += mRhs[i] * 0.05;
    }*/

  //RK4
  calcRHS(mVariables, &mRhs[0]);

  for (i = 0; i < imax; ++i)
    {
      mVar2[i] = mVariables[i] + mRhs[i] * 0.5 * dt;
    }

  calcRHS(mVar2, &mRhsA[0]);

  for (i = 0; i < imax; ++i)
    {
      mVar2[i] = mVariables[i] + mRhsA[i] * 0.5 * dt;
    }

  calcRHS(mVar2, &mRhsB[0]);

  for (i = 0; i < imax; ++i)
    {
      mVar2[i] = mVariables[i] + mRhsB[i] * dt;
    }

  calcRHS(mVar2, &mRhsC[0]);

  for (i = 0; i < imax; ++i)
    {
      mVariables[i] += dt / 6 * (mRhs[i] + 2 * mRhsA[i] + 2 * mRhsB[i] + mRhsC[i]);
    }

  mpLayout->setState(mVariables);
  return mpLayout->getPotential();
}

void CLayoutEngine::EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot)
{static_cast<Data *>((void *) n)->pMethod->evalF(t, y, ydot);}

void CLayoutEngine::evalF(const C_FLOAT64 * /* t */, const C_FLOAT64 * /* y */, C_FLOAT64 * ydot)
{
  //std::cout << *t << "  " << y << " " << &mVariables[0] << std::endl;

// assert(y == &mVariables[0]);

  calcRHS(mVariables, ydot);

  return;
}

void CLayoutEngine::requestStop()
{
  mStopRequested = true;
}
