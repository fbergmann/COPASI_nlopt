#include "copasi/copasi.h"

#ifdef COPASI_USE_NLOPT


#include "COptMethod.h"
#include "COptProblem.h"
#include "COptItem.h"
#include "COptTask.h"

#include "COptMethodNLopt.h"

#include "copasi/math/CMathContainer.h"
#include "copasi/core/CDataObjectReference.h"
#include "copasi/randomGenerator/CRandom.h"

#include <nlopt.hpp>

const CEnumAnnotation< std::string, COptMethodNLopt::NLoptMethodType > COptMethodNLopt::NLOptMethods(
  {
    // global
    "DIRECT (Dividing RECTangles)",
    "DIRECTL (DIRECT with Local Biasing)",
    "CRS (Controlled Random Search)",
    "MLSL (Multi-Level Single Linkage)",
    "MLSL_LDS (Multi-Level Single Linkage with Low Discrepancy Sampling)",
    "StoGO (Stochastic Global Optimization)",
    "AGS (Adaptive Global Search)",
    "ISRES (Improved Stochastic Ranking Evolution Strategy)",
    "ESCH (evolutionary algorithm)",
    // local derivative-free
    "COBYLA (Constrained Optimization BY Linear Approximations)",
    "BOBYQA",
    "NEWUOA+",
    "PRAXIS (PRincipal AXIS)",
    "Nelder-Mead Simplex",
    "Sbplx (based on Subplex)",
    // local gradient-based
    "MMA (Method of Moving Asymptotes)",
    "SLSQP (Sequential Least Squares Programming)",
    "Low-storage BFGS",
    "Truncated Newton",
    "Truncated Newton Restart",
    "Truncated Newton Preconditioned",
    // hybrid
    "AUGLAG (Augmented Lagrangian algorithm)",
  });

COptMethodNLopt::COptMethodNLopt(const CDataContainer * pParent,
                             const CTaskEnum::Method & methodType,
                             const CTaskEnum::Task & taskType)
  : COptMethod(pParent, methodType, taskType, false)
  , mIterations(100000)
  , mCurrentIteration(0)
  , mIndividual()
  , mValue(std::numeric_limits< C_FLOAT64 >::quiet_NaN())
  , mVariableSize(0)
{
  
  initObjects();
}

COptMethodNLopt::COptMethodNLopt(const COptMethodNLopt & src,
                             const CDataContainer * pParent)
  : COptMethod(src, pParent, false)
  , mIterations(src.mIterations)
  , mCurrentIteration(src.mCurrentIteration)
  , mIndividual(src.mIndividual)
  , mValue(src.mValue)
  , mVariableSize(src.mVariableSize)
{initObjects();}

/**
 * Destructor
 */
COptMethodNLopt::~COptMethodNLopt()
{
  //*** added similar to coptga
  cleanup();
}

void COptMethodNLopt::initObjects()
{
  mpNloptMethod = assertParameter("NLopt Method", CCopasiParameter::Type::STRING, NLOptMethods[COptMethodNLopt::NLoptMethodType::NELDER_MEAD]);
  getParameter("NLopt Method")->setValidValues(NLOptMethods);

  addObjectReference("Current Iteration", mCurrentIteration, CDataObject::ValueInt);

  assertParameter("Number of Iterations", CCopasiParameter::Type::UINT, (unsigned C_INT32) 2000);
  
  mpRtolObjective = assertParameter("Relative Tolerance Objective", CCopasiParameter::Type::DOUBLE, HUGE_VAL);
  mpRtolParameters = assertParameter("Relative Tolerance Parameters", CCopasiParameter::Type::DOUBLE, HUGE_VAL);  

  mpAtolObjective = assertParameter("Absolute Tolerance Objective", CCopasiParameter::Type::DOUBLE, HUGE_VAL);
  mpAtolParameters = assertParameter("Absolute Tolerance Parameters", CCopasiParameter::Type::DOUBLE, HUGE_VAL);  

  assertParameter("Initial Step Size", CCopasiParameter::Type::DOUBLE, HUGE_VAL);

  assertParameter("Population Size", CCopasiParameter::Type::UINT, (unsigned C_INT32) 0);

  mpLocalMethod = assertParameter("Local Optimization Algorithm", CCopasiParameter::Type::STRING, std::string("None"));

  std::vector< std::pair< std::string, std::string > > ValidValues = {
    {  "None", "None"},
    // local derivative-free
  { "COBYLA (Constrained Optimization BY Linear Approximations)", "COBYLA (Constrained Optimization BY Linear Approximations)"},
  { "BOBYQA","BOBYQA"},
  {  "NEWUOA+", "NEWUOA+"},
  { "PRAXIS (PRincipal AXIS)", "PRAXIS (PRincipal AXIS)"},
  { "Nelder-Mead Simplex", "Nelder-Mead Simplex"},
  { "Sbplx (based on Subplex)", "Sbplx (based on Subplex)"},
  // local gradient-based
  { "MMA (Method of Moving Asymptotes)", "MMA (Method of Moving Asymptotes)"},
  { "SLSQP (Sequential Least Squares Programming)", "SLSQP (Sequential Least Squares Programming)"},
  { "Low-storage BFGS", "Low-storage BFGS"},
  { "Truncated Newton", "Truncated Newton"},
  { "Truncated Newton Restart", "Truncated Newton Restart"},
  { "Truncated Newton Preconditioned", "Truncated Newton Preconditioned"},
  };


  getParameter("Local Optimization Algorithm")->setValidValues(ValidValues);

  assertParameter("Local Max FnEvals", CCopasiParameter::Type::INT, -1);

  assertParameter("Local Relative Tolerance Objective", CCopasiParameter::Type::DOUBLE, HUGE_VAL);
  assertParameter("Local Relative Tolerance Parameters", CCopasiParameter::Type::DOUBLE, HUGE_VAL);

  assertParameter("Local Absolute Tolerance Objective", CCopasiParameter::Type::DOUBLE, HUGE_VAL);
  assertParameter("Local Absolute Tolerance Parameters", CCopasiParameter::Type::DOUBLE, HUGE_VAL);
  
  mpMonitorOnImprovement = assertParameter("Monitor only on improvement", CCopasiParameter::Type::BOOL, true);

}

/**
 * Optimizer Function
 * Returns: true if properly initialized
 * should return a boolean
 */
bool COptMethodNLopt::initialize()
{
  cleanup();

  if (!COptMethod::initialize()) return false;

  mIterations = getValue< unsigned C_INT32 >("Number of Iterations");

  mVariableSize = mProblemContext.active()->getOptItemList(true).size();
  mIndividual.resize(mVariableSize);

  mPopulationSize = getValue< unsigned C_INT32 >("Population Size");

  return true;
}

/**
 * Objective function wrapper for NLopt
 * This static function is called by NLopt and evaluates the objective function
 */
/*static*/ double COptMethodNLopt::nlopt_objective_function(unsigned n, const double * x, double * grad, void * data)
{
  // Cast the data pointer back to COptMethodNLopt
  COptMethodNLopt *method = static_cast<COptMethodNLopt *>(data);
  
  if (!method->proceed())
  {
    if (!method->mDidThrow)
    {
      method->mDidThrow = true;
      throw nlopt::forced_stop();
    }
    return HUGE_VAL;
  }
  
  
  // Get the optimization problem context
  COptProblem *problem =  method->mProblemContext.active();
  const std::vector<COptItem *> &optItems = problem->getOptItemList(true);
  
  // Set the parameter values from x
  bool validPoint = true;
  for (unsigned i = 0; i < n; i++)
    {
      double value = x[i];
      method->mIndividual[i] = value;
      validPoint &= optItems[i]->setItemValue(value, COptItem::CheckPolicyFlag::None);
    }
  
  // If point is invalid, return a large penalty value
  if (!validPoint)
    {
      return std::numeric_limits<double>::max();
    }

  // Evaluate the objective function
  C_FLOAT64 value = method->evaluate(COptMethod::EvaluationPolicy::Constraints);


  // if gradient is present, compute it (not implemented here)
  if (grad != nullptr)
    {
      // Calculate the gradient
      for (int i = 0; i < n && method->proceed(); i++)
        {
          COptItem & OptItem = *optItems[i];

          if (x[i] != 0.0)
            {
              C_FLOAT64 X = x[i] * 1.001;
              OptItem.setItemValue(X, COptItem::CheckPolicyFlag::None);
              grad[i] = (method->evaluate(EvaluationPolicyFlag::All) - value) / (x[i] * 0.001);
            }
          else
            {
              // why use 1e-7? shouldn't this be epsilon, or something like that?
              C_FLOAT64 X = 1e-7;
              OptItem.setItemValue(X, COptItem::CheckPolicyFlag::None);
              grad[i] = (method->evaluate(EvaluationPolicyFlag::All) - value) / 1e-7;

              if (method->mLogVerbosity > 2)
                {
                  std::ostringstream auxStream;
                  auxStream << "Calculating gradient for zero valued parameter " << i << ", using 1e-7, results in " << grad[i] << ".";
                  method->mMethodLog.enterLogEntry(COptLogEntry(auxStream.str()));
                }
            }

          double val = x[i];
          OptItem.setItemValue(val, COptItem::CheckPolicyFlag::None);
        }
    }

  method->mValue = value;

  if (!method->mpMonitorOnImprovement)
  method->mpParentTask->output(COutputInterface::MONITORING);
  
  if (value < method->getBestValue())
    {
      method->setSolution(value, method->mIndividual, true);
      if (method->mpMonitorOnImprovement)
      method->mpParentTask->output(COutputInterface::MONITORING);
    }

  return value;
}


nlopt::algorithm COptMethodNLopt::methodToAlgorithm(const std::string & method, nlopt::algorithm defaultAlgorithm)
{
  nlopt::algorithm alg = defaultAlgorithm;
  if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::NELDER_MEAD])
    alg = nlopt::LN_NELDERMEAD;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::PRAXIS])
    alg = nlopt::LN_PRAXIS;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::SBPLX])
    alg = nlopt::LN_SBPLX;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::ISRES])
    alg = nlopt::GN_ISRES;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::LBFGS])
    alg = nlopt::LD_LBFGS;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::TRUNCATED_NEWTON])
    alg = nlopt::LD_TNEWTON;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::TRUNCATED_NEWTON_RESTART])
    alg = nlopt::LD_TNEWTON_RESTART;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::TRUNCATED_NEWTON_PRECOND])
    alg = nlopt::LD_TNEWTON_PRECOND;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::COBYLA])
    alg = nlopt::LN_COBYLA;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::BOBYQA])
    alg = nlopt::LN_BOBYQA;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::NEWUOA])
    alg = nlopt::LN_NEWUOA;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::MMA])
    alg = nlopt::LD_MMA;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::SLSQP])
    alg = nlopt::LD_SLSQP;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::DIRECT])
    alg = nlopt::GN_DIRECT;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::DIRECTL])
    alg = nlopt::GN_DIRECT_L;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::CRS])
    alg = nlopt::GN_CRS2_LM;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::StoGO])
    alg = nlopt::GD_STOGO;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::AGS])
    alg = nlopt::GN_AGS;

  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::MLSL])
    alg = nlopt::GN_MLSL;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::MLSL_LDS])
    alg = nlopt::GN_MLSL_LDS;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::AUGLAG])
    alg = nlopt::LN_AUGLAG;
  else if (method == NLOptMethods[COptMethodNLopt::NLoptMethodType::ESCH])
    alg = nlopt::GN_ESCH;


  return alg;
}

bool needsDerivatives(nlopt::algorithm alg)
{
  switch (alg)
    {
      case nlopt::LD_LBFGS:
      case nlopt::LD_TNEWTON:
      case nlopt::LD_TNEWTON_RESTART:
      case nlopt::LD_TNEWTON_PRECOND:
      case nlopt::LD_MMA:
      case nlopt::LD_SLSQP:
        return true;
      default:
        return false;
    }
}

/**
 * Optimizer Function
 * Returns: nothing
 * should return a boolean
 */
//C_INT32 CRandomSearch::optimize()
bool COptMethodNLopt::optimise()
{
  unsigned C_INT32 j;

  // current value is the initial guess
  bool pointInParameterDomain = true;
  
  mDidThrow = false;

  if (!initialize()) return false;

  mMethodLog.enterLogEntry(
      COptLogEntry(
        "Algorithm started.",
        "For more information about this method see: http://copasi.org/Support/User_Manual/Methods/Optimization_Methods/NLopt/"
      )
    );

  // check initial guess
  for (j = 0; j < mVariableSize; j++)
    {
      C_FLOAT64 & mut = mIndividual[j];
      COptItem & OptItem = *mProblemContext.active()->getOptItemList(true)[j];

      mut = OptItem.getStartValue();
      pointInParameterDomain &= OptItem.setItemValue(mut, COptItem::CheckPolicyFlag::All);
      pointInParameterDomain &= mut == OptItem.getStartValue();
    }

  if (!pointInParameterDomain && (mLogVerbosity > 0))
    mMethodLog.enterLogEntry(COptLogEntry("Initial point outside parameter domain."));
  
  // evaluate intitial guess and set as best solution
  mValue = evaluate(EvaluationPolicy::Constraints);
  setSolution(mValue, mIndividual, true);

  // Set up NLopt optimizer
  try
    {
      auto localAlg = methodToAlgorithm(*mpLocalMethod, nlopt::NUM_ALGORITHMS);
      auto alg = methodToAlgorithm(*mpNloptMethod, nlopt::LN_NELDERMEAD);
      bool localNeedsDerivatives = needsDerivatives(localAlg);
      bool needToSetLocal = localAlg != nlopt::NUM_ALGORITHMS;
      if (localNeedsDerivatives && localAlg != nlopt::NUM_ALGORITHMS)
        {
          switch (alg)
            {
            case nlopt::GN_MLSL:
              alg = nlopt::GD_MLSL;
              needToSetLocal = true;
              break;
            case nlopt::GN_MLSL_LDS:
              alg = nlopt::GD_MLSL_LDS;
              needToSetLocal = true; 
              break;
            case nlopt::AUGLAG:
              alg = nlopt::LD_AUGLAG;
              needToSetLocal = true;
              break;
            default:
              break;
            }
        }

      nlopt::opt opt(alg, (int)mVariableSize);

      mMethodLog.enterLogEntry(COptLogEntry(nlopt::algorithm_name(alg)));

      // Set lower and upper bounds
      std::vector<double> lower_bounds(mVariableSize);
      std::vector<double> upper_bounds(mVariableSize);
      std::vector<double> initial_values(mVariableSize);
      
      const std::vector<COptItem *> &OptItemList = mProblemContext.active()->getOptItemList(true);
      
      for (j = 0; j < mVariableSize; j++)
        {
          COptItem &OptItem = *OptItemList[j];
          lower_bounds[j] = *OptItem.getLowerBoundValue();
          upper_bounds[j] = *OptItem.getUpperBoundValue();
          initial_values[j] = mIndividual[j];
        }
      
      opt.set_lower_bounds(lower_bounds);
      opt.set_upper_bounds(upper_bounds);
      
      // Set the objective function
      opt.set_min_objective(nlopt_objective_function, this);
      
      double initialStepSize = getValue< double >("Initial Step Size");
      if (std::isfinite(initialStepSize))
        opt.set_initial_step(initialStepSize);

      
      if (needToSetLocal)
        {
          auto local = nlopt::opt(localAlg, (int)mVariableSize);
          local.set_lower_bounds(lower_bounds);
          local.set_upper_bounds(upper_bounds);

          double value = getValue< double >("Local Relative Tolerance Objective");
          if (std::isfinite(value))
          local.set_ftol_rel(value);
          value = getValue< double >("Local Relative Tolerance Parameters");
          if (std::isfinite(value))
          local.set_xtol_rel(value);

          value = getValue< double >("Local Absolute Tolerance Objective");
          if (std::isfinite(value))
          local.set_ftol_abs(value);
          value = getValue< double >("Local Absolute Tolerance Parameters");
          if (std::isfinite(value))
          local.set_xtol_abs(value);

          int nValue =getValue< int >("Local Max FnEvals");
          if (nValue >=0)
          local.set_maxeval(nValue);

          // Set the objective function
          local.set_min_objective(nlopt_objective_function, this);
      
          opt.set_local_optimizer(local);
          mMethodLog.enterLogEntry(COptLogEntry("Using local optimizer: " + std::string(nlopt::algorithm_name(localAlg))));
        }


      // Set stopping criteria
      opt.set_maxeval(mIterations);

      // set population size if applicable
      if (mPopulationSize > 0)
        opt.set_population(mPopulationSize);
      
      // set relative tolerance on function value
      if (std::isfinite(*mpRtolObjective))
      opt.set_ftol_rel(*mpRtolObjective);
      if (std::isfinite(*mpAtolObjective))
        opt.set_ftol_abs(*mpAtolObjective);
      
      // set relative tolerance on optimization parameters
      if (std::isfinite(*mpRtolParameters))
      opt.set_xtol_rel(*mpRtolParameters);
      if (std::isfinite(*mpAtolParameters))
        opt.set_xtol_abs(*mpAtolParameters);
      
      // Run the optimization
      double minf = mValue;
      nlopt::result result = opt.optimize(initial_values, minf);
      
//      
//      try {
//        std::stringstream init_step;
//        init_step << "[";
//        for (double v : opt.get_initial_step())
//          init_step << v << ",";
//        init_step << "]";
//        
//        
//        mMethodLog.enterLogEntry(COptLogEntry("Used initial step: " + init_step.str()));
//        
//      } catch (...) {
//      }
      
      // Update the solution with NLopt's result
      for (j = 0; j < mVariableSize; j++)
        {
          mIndividual[j] = initial_values[j];
        }
      
      mValue = minf;
      setSolution(mValue, mIndividual, true);
      
      // Get the number of evaluations performed
      mCurrentIteration = opt.get_numevals();
      
      // Log the result
      if (mLogVerbosity > 0)
        {
          std::string resultMsg;
          switch (result)
            {
              case nlopt::SUCCESS:
                resultMsg = "Generic success";
                break;
              case nlopt::STOPVAL_REACHED:
                resultMsg = "Stop value reached";
                break;
              case nlopt::FTOL_REACHED:
                resultMsg = "Function tolerance reached";
                break;
              case nlopt::XTOL_REACHED:
                resultMsg = "Parameter tolerance reached";
                break;
              case nlopt::MAXEVAL_REACHED:
                resultMsg = "Maximum evaluations reached";
                break;
              case nlopt::MAXTIME_REACHED:
                resultMsg = "Maximum time reached";
                break;
              default:
                resultMsg = "Optimization stopped";
                break;
            }
          
          mMethodLog.enterLogEntry(
            COptLogEntry("Algorithm finished.",
                         resultMsg + ". Terminated after " + std::to_string(mCurrentIteration) + " evaluations."));
        }
      
      mpParentTask->output(COutputInterface::MONITORING);
      
      return true;
    }
  catch (std::exception &e)
    {
      mMethodLog.enterLogEntry(COptLogEntry("NLopt error", e.what()));      
      return false;
    }
}

unsigned C_INT32 COptMethodNLopt::getMaxLogVerbosity() const
{
  return 0;
}

C_FLOAT64 COptMethodNLopt::getCurrentValue() const
{
  return mValue;
}

const CVector< C_FLOAT64 > * COptMethodNLopt::getBestParameters() const
{
  return &mIndividual;
}

const CVector< C_FLOAT64 > * COptMethodNLopt::getCurrentParameters() const
{
  return &mIndividual;
}


#endif // COPASI_USE_NLOPT
