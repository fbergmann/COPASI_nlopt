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

// Copyright (C) 2008 - 2009 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2002 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 * COptMethodNLopt class
 */

 #ifndef COPASI_COptMethodNLopt
 #define COPASI_COptMethodNLopt
 
 #include "copasi/optimization/COptMethod.h"
 #include "copasi/core/CVector.h"
 
#ifdef COPASI_USE_NLOPT

#  include <nlopt.hpp>

 class COptMethodNLopt : public COptMethod
 {
   // Operations
 private:
   /**
    * Default Constructor
    */
   COptMethodNLopt();
 
   /**
    * Initialize arrays and pointer.
    * @return bool success
    */
   bool initialize() override;
 
   /**
    * Initialize contained objects.
    */
   void initObjects();
 
   /**
    * Find the best individual at this generation
    * @return size_t fittest
    */
   size_t fittest();
 
   /**
    * number of iterations
    */
   unsigned C_INT32 mIterations;
 
   /**
    * The current iteration
    */
   unsigned C_INT32 mCurrentIteration;
 
   /**
    * a vector of the number of individuals.
    */
   CVector<C_FLOAT64> mIndividual;
 
   /**
    * array of values of objective function f/ individuals
    */
   C_FLOAT64 mValue; 
 
   /**
    * *** Perhaps this is actually not needed ****number of parameters
    */
   size_t mVariableSize;

   double * mpRtolObjective;
   double * mpRtolParameters;

   double * mpAtolObjective;
   double * mpAtolParameters;

   /**
    * population size
    */
   unsigned C_INT32 mPopulationSize;

   std::string * mpNloptMethod;

   std::string * mpLocalMethod;

   bool mDidThrow {false};
   
   bool *mpMonitorOnImprovement;
   
  enum NLoptMethodType
   {
    // global
     DIRECT = 0,
     DIRECTL,
     CRS,
     MLSL,
     MLSL_LDS,
     StoGO,
     AGS,
     ISRES,
     ESCH,
     // local derivative free
     COBYLA,
     BOBYQA,
     NEWUOA,
     PRAXIS,
     NELDER_MEAD,
     SBPLX,
     // local gradient based
     MMA,
     SLSQP,
     LBFGS,
     TRUNCATED_NEWTON,
     TRUNCATED_NEWTON_RESTART,
     TRUNCATED_NEWTON_PRECOND,
     // hybrid
     AUGLAG,
     

     __SIZE
   };

   static const CEnumAnnotation< std::string, NLoptMethodType > NLOptMethods;


 public:
   /**
    * Specific constructor
    * @param const CDataContainer * pParent
    * @param const CTaskEnum::Method & methodType (default: NLopt)
    * @param const CTaskEnum::Task & taskType (default: optimization)
    */
   COptMethodNLopt(const CDataContainer * pParent,
                 const CTaskEnum::Method & methodType = CTaskEnum::Method::NLopt,
                 const CTaskEnum::Task & taskType = CTaskEnum::Task::optimization);
 
   /**
    * Copy Constructor
    * @param const COptMethodNLopt & src
    * @param const CDataContainer * pParent (default: NULL)
    */
   COptMethodNLopt(const COptMethodNLopt & src,
                 const CDataContainer * pParent);
 
   /**
    * Destructor
    */
   virtual ~COptMethodNLopt();
 
   /**
    * Execute the optimization algorithm calling simulation routine
    * when needed. It is noted that this procedure can give feedback
    * of its progress by the callback function set with SetCallback.
    */
   bool optimise() override;
 
   /**
    * Returns the maximum verbosity at which the method can log.
    */
   unsigned C_INT32 getMaxLogVerbosity() const override;
 
   C_FLOAT64 getCurrentValue() const override;
 
   const CVector< C_FLOAT64 > * getBestParameters() const override;
 
   const CVector< C_FLOAT64 > * getCurrentParameters() const override;

   private:
   /**
    * The objective function for nlopt
    */
    static double nlopt_objective_function(unsigned n, const double * x, double * grad, void * data);

    static nlopt::algorithm methodToAlgorithm(const std::string & method, nlopt::algorithm defaultAlgorithm);
 };
 #endif // COPASI_USE_NLOPT

 #endif  // COPASI_COptMethodNLopt
 
