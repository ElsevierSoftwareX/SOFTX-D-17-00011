/* compile: [g,i]cc compuGUTv01.c -lsundials_kinsol \
                                  -lsundials_cvode \
                                  -lsundials_nvecserial \
                                  -lm
 * -----------------------------------------------------------------------------
 * compuGUT: In silico platform for the design and simulation of intestinal
 * fermentation experiments
 *
 * Revision: 0.2
 * Date: 2015/02/06
 * -----------------------------------------------------------------------------
 * Contributors: Arun S. Moorthy and Hermann J. Eberl
 * Affiliation:  Biophysics Interdepartmental Group, Department of Mathematics
 *               and Statistics, University of Guelph
 *
 * Contact:  [amoorthy, heberl]@uoguelph.ca
 *           http://compugut.sourceforge.net/
 * 
 * Past Contributors: Jesse Knight (2013)
 *                    Kathleen Songin (2014)
 *                    Richard Yam (2014)
 *
 * File Description:
 * The following is the single-file source code behind the compuGUT operation.
 * It follows a non-linear structure, presenting users with many options in
 * terms of problem design, and solution scheme. To summarize:
 *
 *
 * 1)Collect simulation conditions:
 *   (a) Operation Instructions: This includes colon sizing, and bacterial
 *       representation (between 1 and 10 representatives within a functional
 *       group). Additional inputs include the level of accuracy for
 *       integration scheme (grid resolution, k), biological and physical
 *       parameter variances, and print frequency options.
 *   (b) Input Meal/Feeding Strategy (i) continuous feeding of known
 *       concentrations of fiber (~prebiotics), and biomass (~probiotics).
 *       (ii) customized feeding of known concentrations of fiber (~prebiotics)
 *       and biomass (~probiotics).
 *   (c) Initial Conditions: (i) read in solution file (*_r*.txt) from a
 *       previous simulation, using linear interpolation to ensure
 *       compatibility with solution strategy. (ii) read in intial conditions
 *       as average values from file `user_colon_initializer.txt'.
 *
 * 2) Collect simulation type:
 *   (a) Continuous colon: Based on the model of Moorthy et. al (2015),
 *       describing the colon as a continuous pipe/space. Results in a system
 *       of 20 - 100 PDEs, depending on the biological representation
 *       simulated. Unlike results that are obtainable using current (2015) in
 * 	     vivo and in vitro techniques. Files produced using this simulation
 *       model will be prefixed with `ccm'.
 *   (b) Three-stage reactor system: Based on the model of
 *       Munoz-Tamayo et. al (2010), modelling the colon as a system of
 *       sequential reactors. Results in a system of 60 - 300 ODEs, depending
 *       on the biological representation simulated. Most comparible to in
 *       vitro simulations currently employed in literature. Files produced
 * 	 using this simulation model will be prefixed with `tsr'.
 *   (c) Gradostat reactor system: System is modeled as a series of `m' equal
 *       volume reactors. Results in a system of (20-100)*m ODEs, depending on
 *       the biological representation simulated. Files produced using this
 *       simulation model will be prefixed with `gsm'.
 *   (d) All of three methods can be simulated simultaneously, best used for
 *       methods comparisons.
 *
 * 3) Parameter Loading:
 *   (a) ADM1 based Reaction parameters, Munoz-Tamayo et al based Transport
 *       parameters. Loaded and slotted into 4 matrices:
 * 	 (i) Peterson Matrix of Yield Coefficients (called as `scope'_PMat)
 *       (ii) Kinetic Rates Matrix (called as `scope'_KRMat)
 *       (iii) Half-Saturation Constant Matrix (called as `scope'_KHSMat)
 *       (iv) Transportation Matrix (called as `scope'_TMat)
 *
 * 4) Initialization:
 *    Prepare for main simulation loop (define loop variables,
 *    stopping criteria, etc.).
 *
 * 5) Integration:
 *    For simulation type (a), we employ KINSol 2.7.0 from the Sundials Suite
 *    (see below) for solving the systems of non-linear equations that arise
 *    through numerical integration. Currently, integration of the system of
 *    PDEs is accomplished using a central scheme for balance laws as presented
 *    in Liotta et. al. (2001), however, we are looking to employ other
 *    integration schemes in the near future. All systems of ODEs are solved
 *    using cvODE 2.7.0, and all non-linear systems of equations are solved
 *    using KINSol 2.7.0, also from the Sundials Suite (see below).
 *
 * 6) Printing:
 *    As noted, the prefix of the output textfile will indicate which model had
 *    created the solution (i.e. tsr_r1.txt from three stage model, ccm_r1.txt
 *    from continuous colon model). All output files will have the same
 *    structure, An initial header with the time stamp of the output, and a
 *    text delimited matrix of size n x m, where n is the number of discrete
 *    grid points (3 for tsm, 51*pow(2,k) for ccm) and m is the number of state
 *    variables (60:300 for tsm, 20:100 for ccm)
 *
 * Notes:
 * - Global variables are named with the following convention: g_varname_unit.
 *   If the variable has no obvious units, for example, g_nsv, which is a count
 *   of the number of state variables, will not have a unit suffix.
 * - Global variables are copied into modules and given the following naming
 *   convention: `scope'_varname_unit.
 * ---------------------------------------------------------------------------
 *
 * External Libraries Used:
 * Sundials Suite, as produced by:
 * Center for Applied Scientific Computing, Lawrence
 * Livermore National Laboratory
 * http://computation.llnl.gov/casc/sundials/main.html
 *
 * Specific Libraries:
 *  - Sundials CVODE v2.7.0
 *  - Sundials KINSOL v2.7.0
 *
 * External Functions Used:
 * To quickly allocate 2D arrays, the following blog post was used:
 * http://pleasemakeanote.blogspot.ca/2008/06/2d-arrays-in-c-using-malloc.html
 * -----------------------------------------------------------------------------
 *
 * Complimentary Tools:
 * Documentation for these tools is not guaranteed to be complete, rather, are
 * current works in progress. However, they should be of benefit to first-time
 * users.
 *
 * user_prompt.r - R program that wraps compuGUT tool, allowing for users to
 *                 conveniently edit input files, execute compuGUT, and create
 *                 basic plots. Developed August 2014 by Kathleen Songin and
 *                 Richard Yam during an Undergraduate Summer Research Term.
 *
 * batch_op.sh  - bash script that runs compuGUT binary in serial, for use with
 *                larger scale simulation experiments (with stochastic
 *                parameters)
 *
 * quick-compile.sh - bash script to quickly compile compuGUT source code with
 *                    sundials libraries
 *
 * Other tools described and/or used in our below cited papers can be acquired
 * through request.
 *
 * ----------------------------------------------------------------------------- 	
 *
 * Mathematical/Model References:
 * Munoz-Tamayo et. al. Mathematical modelling of carbohydrate degradation by
 * human colonic microbiota. Journal of Theoretical Biology (2010)
 *
 * Moorthy & Eberl. Assessing the influence of reactor design criteria on the
 * performance of model colon fermentation units. Journal of Bioscience and
 * Bioengineering (2014)
 *
 * Moorthy, Brooks, Kalmokoff & Eberl. A spatially continuous model of 
 * carbohydrate digestion and transport processes in the large intestine. 
 * (2015) - pre-print available on request
 *
 * Software Reference:
 * Moorthy & Eberl. compuGUT: An in silico tool for studying flora composition
 * in the large intestine (2015) - pre-print available on request
 *
 * User documentation available at http://compugut.sourceforge.net
 *
 * -arun
 *
 * ---------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// standard headers:
#include <stdlib.h>                  /* standard library */
#include <stdio.h>                   /* input/output */
#include <math.h>                    /* mathematics functions */
#include <string.h>                  /* creation of filenames */
#include <time.h>                    /* assessing computing times */
#include <errno.h>

// sundials headers:
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <kinsol/kinsol.h>           /* prototypes for KINSol fcts., consts. */
#include <kinsol/kinsol_dense.h>     /* prototype for KINDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definitions for math functions */

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// DEFINITIONS:
#define FTOL  RCONST(1e-08)      /* function tolerance (KINSol)*/
#define STOL  RCONST(1e-08)      /* step tolerance (KINSol)*/
#define RTOL  RCONST(1e-10)      /* scalar relative tol. (CVODE) */
#define ATOLs RCONST(1e-10)      /* abs. tol. small initial guess (CVODE) */
#define ATOLb RCONST(1e-2)       /* abs. tol. large initial guess (CVODE) */
#define EPSI  RCONST(1e-16) 	   /* for contois rate equation */
#define PI    RCONST(3.1415926)	 /* value of pi */
#define E     RCOST(2.7182818)	 /* value of e */

// sundials macros:
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

// Reaction/Situation-Specific Definitions:
#define NOTHERS 10      /* Number of non-biomass state variables */
#define MPS     100     /* Maximum problem size */
#define MB      10      /* Maximum number per functional biomass group */
#define MR      52001   /* Maximum grid resolution */
#define MinR 	51      /* Minimum grid resolution */
#define NAP     9       /* Number of ADM processes */

// Colon specific
#define PT 0.14	   /* Perc. of total length colon where proximal->transverse */
#define TD 0.42	   /* Perc. of total length colon where transverse->distal */
#define SpA 0.1    /* Percentage splined */

// Diet specific
#define MEALLENGTH 0.25  /* Length of a meal in hours */
#define MMD	   1704  /* Maximum meal data 168 hours = 7 days */

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// STRUCTURES: -----------------------------------------------------------------
typedef struct{
  /* User Operation Variables are those that define the type of simulation,
   * volume/sizing of the physical problem, and the length of simulation time.
   * Also included are computational variables, grid resolution: the level of
   * spatial discretization used in full colon simulations, and print-frequency:
   * the number of iterations between solution printing.
   *
   * call structure in scope as: `scope'_op_inst
   */

  int 	    simtype;
  double    lcolon_m;
  double    dcolon_m;
  double    lsintestine_m;
  double    dsintestine_m;
  double    otime_d;

  int       grid_resolution;
  int       print_frequency;
} *UOperationVariables;

typedef struct {
  /* Reaction Data Variables are those that are required to compute the yield of
   * the various reaction and transport equations that are presented in the
   * mathematical problem.
   *
   * call structure as: `scope'_reac_data
   */

  double** PMat;
  double** KRMat;
  double** KHSMat;
  double** TMat;
  double* pH_v;
  double dx;
  double dt;
  double LocX;
  double tfac;
  double* interim;
} *ReactionData;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// GLOBAL VARIABLES:	Initialization Value (over-written in ROI function) --
int g_nsv       = MPS;
int g_mbacteria = MB;
int g_nsd       = MB;
int g_nld       = MB;
int g_nhda      = MB;
int g_nhdm      = MB;

double g_fr_lpd	        = 0;
double g_fr_mpd         = 0;
double g_lcolon_m       = 0;
double g_vcolon_l       = 0;
double g_vproximal_l    = 0;
double g_vtransverse_l  = 0;
double g_vdistal_l      = 0;
double g_vsintest_l     = 0;

double g_bvariance      = 0;
double g_pvariance      = 0;

int gLocation           = 0;

int g_num_scheme 	= 1; /* 1-Liota et. al, 2-MOL/upwinding*/
int g_MOL_sch_order 	= 2; /* 1-first order, 2-second order, for MOL */

int g_hydro_scheme 	= 1; /* 1-contois, 2-linear, 3-monod */

int g_nreactors = 50;

int g_MealData = MMD; /* size of meal data in hours */
int g_sdiv = 3;	      /* number of hydrographs - pre-colon black box */

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// FUNCTION DECLARATIONS: ------------------------------------------------------
// Helper Functions:
int PauseMessage(int a);

int Clean1DArray(double* array, int ncol);

int Clean2DArray(double** array, int nrows, int ncol);

double** Make2DDoubleArray(int arraySizeX, int arraySizeY);

int CodeCheckParameters(void* user_param_data);

// Primary Modules:
int ReadOperationInstructions(void* user_data); // ROI

double NormRand(double std_dev, double mean, int seed); //NR

static int PMG(double** PMat, double* DefaultYParameters);

static int KMG(double** KRateMatrix, double** KHSConstant,
  double* KRateParameters_vector);

static int TMG(double** TransportParameters, double* Parameters);

int ReadReactionParameters(void* user_data); //RRP

static int InitialConditions(double** PSolution, int ngrids); //IC

int PrintFullSolution(double** solution, int ngrids, char a[],
  int fileiter, double t);

int PrintCondensedSolution(double** solution, int ngrids, char a[],
  int fileiter, double t, int increment);

int MOL_numScheme(double** PSolution, double* Input, int ngrids,
     void* user_reac_data);

int ThreeStageReactorScheme(double **PSolution, double *Input,
  void* user_reac_data);

int GradostatScheme(double **PSolution, double *Input, int nReactors,
  void* user_reac_data);

int Liotta_numScheme(double** PSolution, double* Input, int ngrids,
  void* user_reac_data);

int LoadMealPlan(double** mealPlan);

static int PreGut(double* input, double** PSol, double** mealPlan,
  double time, double dt);
  
static int PGUT_interpolator(double** mealPlan, double mealTime, 
  double SimTime, double** OutputMeal, int MaxIter, double dt);
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main()
{
  // Main scope indexing variables:
  int indexi,indexj,indexk;
  int aflag;

  // Initialize Program:
  clock_t runtime = clock();
  printf(\
  "\ncompuGUT: In silico platform for simulating intestinal fermentation\n");
  printf("and microbiota response\n");
  printf("Version: 0.1 \n");
  printf("Copyright: 2015, Moorthy & Eberl\n\n");

//------------------------------------------------------------------------------
// Initialize UOV Structure:
	UOperationVariables main_op_inst;
	main_op_inst = NULL;
	main_op_inst = (UOperationVariables)malloc(sizeof *main_op_inst);

	aflag = ReadOperationInstructions(main_op_inst);
//------------------------------------------------------------------------------
// Check that global variables have appropriate values (comment on use):
printf("Problem Size: %d x %d\n",g_nsv,main_op_inst->grid_resolution);
printf("System Flow Rate [L/d]: %f\n",g_fr_lpd);
printf("Convective velocity [m/d]: %f\n",g_fr_mpd);
printf("Colon Volume [L]: %f (%f + %f + %f)\n",\
		g_vcolon_l,g_vproximal_l,g_vtransverse_l,g_vdistal_l);
printf("Small Intestine Volume [L]: %f\n",g_vsintest_l);
//------------------------------------------------------------------------------
// Initialize ReactionData Structure:
	ReactionData main_reac_data;
	main_reac_data = NULL;
	main_reac_data = (ReactionData)malloc(sizeof *main_reac_data);
	main_reac_data->PMat = Make2DDoubleArray(NAP,g_nsv/2);
	main_reac_data->KRMat = Make2DDoubleArray(NAP,g_mbacteria);
	main_reac_data->KHSMat = Make2DDoubleArray(5,g_mbacteria);
	main_reac_data->TMat = Make2DDoubleArray(4,g_nsv/2);
	main_reac_data->pH_v = malloc(3 *sizeof (double));
	main_reac_data->interim = malloc(g_nsv * sizeof (double));

	aflag = ReadReactionParameters(main_reac_data);
        //aflag = CodeCheckParameters(main_reac_data);
//------------------------------------------------------------------------------
// Load meal data;
double** main_meal_data;
main_meal_data = Make2DDoubleArray(MMD,3);
aflag = Clean2DArray(main_meal_data,MMD,3);

aflag = LoadMealPlan(main_meal_data);

//~ //for(indexi=0;indexi<g_MealData;indexi++){
	//~ //for(indexj=0;indexj<3;indexj++){
		//~ //printf("%f ",main_meal_data[indexi][indexj]);
	//~ //}
	//~ //printf("\n");
//~ //}
//~ //PauseMessage(396);
//------------------------------------------------------------------------------

  // Primary Solution variables:
  double dx_m = main_op_inst->lcolon_m/(main_op_inst->grid_resolution - 1);
  double dt_d = dx_m/g_fr_mpd;

  main_reac_data->dx = dx_m;
  main_reac_data->dt = dt_d;

  printf("delx = %f [m] \n",dx_m);
  printf("delt = %f [d] \n",dt_d);

  int 	ngrids = main_op_inst->grid_resolution;
  int 	nreactors = g_nreactors;
  int 	nsv = g_nsv;
  double** tsr_psolution;
  double** ccm_psolution;
  double** gs_psolution;

	  tsr_psolution = Make2DDoubleArray(3, nsv);
	  aflag = Clean2DArray(tsr_psolution,3,nsv);

	  gs_psolution = Make2DDoubleArray(nreactors,nsv);
	  aflag = Clean2DArray(gs_psolution,nreactors,nsv);

	  ccm_psolution = Make2DDoubleArray(ngrids, nsv);
	  aflag = Clean2DArray(ccm_psolution,ngrids,nsv);

  int simType = main_op_inst->simtype;

  // Initial Conditions
  if (simType == 1 || simType == 4){
	  aflag = InitialConditions(ccm_psolution,ngrids);
	}

	if (simType == 2 || simType == 4){
	  aflag = InitialConditions(tsr_psolution,3);
	}

	if (simType == 3 || simType == 4){
		aflag = InitialConditions(gs_psolution,nreactors);
	}

  // iteration variables:
	double time_d = 0;
	double mealADJ_time_d = 0;
	double max_meal_data_time = ceil(g_MealData/24);
	// printf("%f\n",max_meal_data_time);
	// PauseMessage(445);

	double Otime_d = main_op_inst->otime_d;
	int    iter = 0;
	int    MaxIter = Otime_d/dt_d + 1;
	char   a[7]="tsr_r";
	char   b[7]="ccm_r";
	char   c[7]="ccm_i";
	char   d[7]="gsm_r";
	int   fileiter = 0;

	int SolIncrement = (ngrids-1)/(MinR-1);

	// Print Initial Solution
	if (simType == 1 || simType == 4){
	 aflag = PrintFullSolution(ccm_psolution,ngrids,b,fileiter,time_d);
	 if(SolIncrement != 1){
	  aflag = PrintCondensedSolution(ccm_psolution,ngrids,c,fileiter,time_d,
	   SolIncrement);}
	}

	if (simType == 2 || simType == 4){
	  aflag = PrintFullSolution(tsr_psolution,3,a,fileiter,time_d);
	}

	if (simType == 3 || simType == 4){
	 aflag = PrintFullSolution(gs_psolution,nreactors,d,fileiter,time_d);
	}

FILE *f4;
f4 = fopen("TimeValues.txt","w");
fprintf(f4,"%f\n",time_d);

FILE *f5;
f5 = fopen("RecordedTimeValues.txt","w");
fprintf(f5,"%f\n",time_d);


double** InputMeal; 
InputMeal = Make2DDoubleArray(MaxIter,nsv);

aflag = PGUT_interpolator(main_meal_data, \
                          max_meal_data_time, \
                          Otime_d, \
                          InputMeal, \
                          MaxIter, \
                          dt_d);

	double *inputmeal;
	inputmeal = malloc(nsv * sizeof (double));
	aflag = Clean1DArray(inputmeal,nsv);


printf("\nSimulating colon:\n");

//~ // MAIN LOOP
 do{
	
	   //~ mealADJ_time_d = fmod(time_d,max_meal_data_time);
	   //~ aflag = PreGut(inputmeal,pre_gut_solution, \
                  //~ main_meal_data,mealADJ_time_d,dt_d);
                  
   for(indexi=0;indexi<nsv;indexi++){
		 inputmeal[indexi] = InputMeal[iter][indexi];
	 }               
   
   if (simType == 1 || simType == 4){
    switch (g_num_scheme) {
     case 1:
		//~ //    printf("\nUsing central scheme of Liotta (2001) \
     //~ to solve continuous colon model.\n");
     aflag = Liotta_numScheme(ccm_psolution,inputmeal,ngrids,main_reac_data);
     break;
     case 2:
		//~ //     printf("\nUsing Method of Lines with explicity upwinding \
     //~ of order %d to solve continuous colon model.\n",g_MOL_sch_order);
     aflag = MOL_numScheme(ccm_psolution,inputmeal,ngrids,main_reac_data);
     break;
    }
   }

   if (simType == 2 || simType == 4){
   //printf("\nUsing cvODE to solve three stage reactor model.\n");
    aflag = ThreeStageReactorScheme(tsr_psolution,inputmeal,main_reac_data);
   }

   if (simType == 3 || simType ==4){
   //printf("\nUsing cvODE to solve gradostat model.\n");
    aflag = GradostatScheme(gs_psolution,inputmeal,nreactors,main_reac_data);
   }

 iter++;
 time_d = time_d + dt_d;
 fprintf(f4,"%f\n",time_d);




  // Print Solution
  if(iter % main_op_inst->print_frequency == 0 || iter > MaxIter - 2)
  //if(iter > MaxIter - 2)
    {fprintf(f5,"%f\n",time_d);
     fileiter++;
       if (simType == 1 || simType == 4){
	 aflag = PrintFullSolution(ccm_psolution,ngrids,b,fileiter,time_d);
	    if(SolIncrement != 1){
         aflag = PrintCondensedSolution(ccm_psolution,ngrids,c,fileiter,time_d,
		    SolIncrement);}
       }
       if (simType == 2 || simType == 4){
         aflag = PrintFullSolution(tsr_psolution,3,a,fileiter,time_d);
       }
       if (simType == 3 || simType == 4){
         aflag = PrintFullSolution(gs_psolution,nreactors,d,fileiter,time_d);
       }
   }

} while (time_d < Otime_d && iter < MaxIter);
printf("Done!\n");

// End Timer:
  time_t current_time;
  char* string_ctime;
  current_time = time(NULL);
  string_ctime = ctime(&current_time);
  runtime = clock() - runtime;

  fclose(f4);
  fclose(f5);


// Write simulation information to file:
FILE *f2;
f2 = fopen("SimulationSummary.txt","w");
fprintf(f2,"compuGUT_fm Simulation Summary: \n");
fprintf(f2,"Date and time: %s\n\n",string_ctime);
if(simType==1){
  fprintf(f2,"Simulation Type: Continuous Colon Model\n\n");
} else if (simType==2){
  fprintf(f2,"Simulation Type: Three Stage Reactor Model\n\n");
} else if (simType==3){
  fprintf(f2,"Simulation Type: Gradostat Model\n\n");
} else if (simType==4){
  fprintf(f2,"Simulation Type: Comparative Study (all three models)\n\n");
}
fprintf(f2,"Number of State Variables: %d\n", g_nsv);
fprintf(f2,"Number of Sugar Degrading biomass: %d\n",g_nsd);
fprintf(f2,"Number of Lactate Degrading biomass: %d\n",g_nld);
fprintf(f2,"Number of Acetogenic biomass: %d\n",g_nhda);
fprintf(f2,"Number of Methanogenic biomass: %d\n",g_nhdm);
fprintf(f2,"Bioparameter variance: %f\n\n",g_bvariance);
fprintf(f2,"Grid resolution: %d\n", ngrids);
fprintf(f2,"Length of simulated colon [m]: %f\n",g_lcolon_m);
fprintf(f2,"Average diameter of colon [m]: %f\n",main_op_inst->dcolon_m);
fprintf(f2,"Average flowrate [L/d]: %f\n", g_fr_lpd);
fprintf(f2,"\n\nDays simulated: %f\n", Otime_d);
fprintf(f2,"Computing time [s]: %f\n",((float)runtime)/CLOCKS_PER_SEC);
fprintf(f2,"Number of output files created: %d\n", fileiter);
fprintf(f2,"\n\nContact Arun [amoorthy@uoguelph.ca] or Hermann");
fprintf(f2," [heberl@uoguelph.ca] with questions or feedback");
fprintf(f2,"\n\n(c) Moorthy & Eberl, 2015\n\n");
fclose(f2);


//------------------------------------------------------------------------------
// Free Memory:

  for(indexi=0;indexi<MMD;indexi++){
		free(main_meal_data[indexi]);
	}
	free(main_meal_data);

	for(indexi=0;indexi<ngrids;indexi++){
		free(ccm_psolution[indexi]);
	}
	free(ccm_psolution);

	for(indexi=0;indexi<3;indexi++){
		free(tsr_psolution[indexi]);
	}
	free(tsr_psolution);

	for(indexi=0;indexi<nreactors;indexi++){
		free(gs_psolution[indexi]);
	}
	free(gs_psolution);

	for(indexi=0;indexi<MaxIter;indexi++){
		free(InputMeal[indexi]);
	}
	free(InputMeal);
	
	free(inputmeal);

	free(main_op_inst);

	for(indexi=0;indexi<NAP;indexi++){
		free(main_reac_data->PMat[indexi]);
		free(main_reac_data->KRMat[indexi]);
	}
	free(main_reac_data->PMat);
	free(main_reac_data->KRMat);

	for(indexi=0;indexi<5;indexi++){
		free(main_reac_data->KHSMat[indexi]);
	}
	free(main_reac_data->KHSMat);

	for(indexi=0;indexi<4;indexi++){
		free(main_reac_data->TMat[indexi]);
	}
	free(main_reac_data->TMat);
	free(main_reac_data->pH_v);
	free(main_reac_data->interim);
	free(main_reac_data);



return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int PauseMessage(int a)
{
  /* Description: This function is for code checking purposes. If called with
   * line number, can provide a quick way to know where in the source code
   * problems may be occuring.
   *
   *  Date: 2014/07/03

   */

printf("\nMessage %d\nPress enter to continue...\n",a);
getchar();
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int Clean1DArray(double* array, int ncol)
{
  /* Description: This function ensures that all 1d arrays are initialized to
   * zero, as to prevent unpredictable/irradic behaviour after malloc.
   *
   *  Date: 2014/07/03
   */
	int i;
	for(i=0;i<ncol;i++){
		array[i] = 0.0;
	}
	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int Clean2DArray(double** array, int nrows, int ncol)
{
  /* Description: This function ensures that all 2d arrays are initialized to
   * zero, as to prevent unpredictable/irradic behaviour after malloc.
   *
   *  Date: 2014/07/03
   */

	int i;
	int j;

	for(i=0;i<nrows;i++){
		for(j=0;j<ncol;j++){
			array[i][j] = 0.0;
		}
	}

	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double** Make2DDoubleArray(int arraySizeX, int arraySizeY)
{
  /* Description: This function quickly allocates memory for 2d arrays.
   *
   * From:
   * http://pleasemakeanote.blogspot.ca/2008/06/2d-arrays-in-c-using-malloc.html
   *
   *  Date: 2014/07/03
   */

    int i,aflag;
    double** theArray;

    theArray = (double**) malloc(arraySizeX*sizeof(double*));
    for (i = 0; i < arraySizeX; i++){
       theArray[i] = (double*) malloc(arraySizeY*sizeof(double));
    }

    aflag = Clean2DArray(theArray,arraySizeX,arraySizeY);

return theArray;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int CodeCheckParameters(void* user_param_data)
{
  /* Description: This function prints out the values of the various reaction
   * matrices. For error checking.
   *
   * Date: 2014/07/03
   */
	int indexi,indexj,indexk;

	ReactionData CCP_reac_data;
	CCP_reac_data = (ReactionData)user_param_data;

	printf("\nPeterson Matrix: \n\n");
	for(indexi=0;indexi<NAP;indexi++){
		for(indexj=0;indexj<g_nsv/2;indexj++){
			printf("|%1.1f ",CCP_reac_data->PMat[indexi][indexj]);
		}
		printf("\n");
	}
	printf("\nKRate Matrix: \n\n");
	for (indexi=0;indexi<NAP;indexi++){
		for(indexj=0;indexj<g_mbacteria;indexj++){
			printf("%f ",CCP_reac_data->KRMat[indexi][indexj]);
		}
		printf("\n");
	}
	printf("\nKHS Matrix: \n\n");
	for (indexi=0;indexi<5;indexi++){
		for(indexj=0;indexj<g_mbacteria;indexj++){
			printf("%f ",CCP_reac_data->KHSMat[indexi][indexj]);
		}
		printf("\n");
	}
	printf("\nTrans Matrix: \n\n");
	for (indexi=0;indexi<4;indexi++){
		for(indexj=0;indexj<g_nsv/2;indexj++){
			printf("%f ",CCP_reac_data->TMat[indexi][indexj]);
		}
		printf("\n");
	}

	printf("pH values: \n\n");
	printf("Proximal Colon: %f\n",CCP_reac_data->pH_v[0]);
	printf("Transverse Colon: %f\n",CCP_reac_data->pH_v[1]);
	printf("Distal Colon: %f\n",CCP_reac_data->pH_v[2]);

	PauseMessage(7);
	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int ReadOperationInstructions(void* user_data)
{
  /* Description: Reads in user instructions and allocates to the appropriate
   * variables.
   *
   *  Date: 2014/07/03
   */


	int 		indexi;
	int 		max_operation_inputs = 25;
	double 	        *u_operation_inputs;
	const int 	bsz = 300;
	char 		buf[bsz];

	UOperationVariables ROI_op_inst;
	ROI_op_inst = (UOperationVariables)user_data;

  u_operation_inputs = malloc(max_operation_inputs * sizeof (double));

  FILE *f_uop_inputs;
  f_uop_inputs = fopen("InputFiles/user_operation_instr_compuGUT.txt","r");

	  if (f_uop_inputs != NULL){
		indexi = 0;
	    while (fgets(buf, bsz, f_uop_inputs) != NULL) {
	      if (buf[0] != '#' && buf[0]!= '\n') {
		  sscanf(buf,"%lf", &u_operation_inputs[indexi]);
		  indexi++;
	      }
	    }
	  } else {
	    perror("fopen");
	  }
	fclose(f_uop_inputs);

	//~ int npara = indexi;

	//~ for (indexi=0;indexi<npara;indexi++){      //<- error check location
	//~ printf("%d: %lf\n",indexi,u_operation_inputs[indexi]);
        //~ }
	//~ PauseMessage(456);

  int p0 = ROI_op_inst->simtype = (int) u_operation_inputs[0];
  double p1 = g_lcolon_m = ROI_op_inst->lcolon_m = u_operation_inputs[1];
  double p2 = ROI_op_inst->dcolon_m = u_operation_inputs[2];
  double p3 = ROI_op_inst->lsintestine_m = u_operation_inputs[3];
  double p4 = ROI_op_inst->dsintestine_m = u_operation_inputs[4];
  double p5 = g_fr_lpd = u_operation_inputs[5];
  double p6 = ROI_op_inst->otime_d = u_operation_inputs[6];
  int p7 = g_nsd = (int)u_operation_inputs[7];
  int p8 = g_nld = (int)u_operation_inputs[8];
  int p9 = g_nhda = (int)u_operation_inputs[9];
  int p10= g_nhdm = (int)u_operation_inputs[10];

  double p11= g_bvariance = u_operation_inputs[11];
  double p12= g_pvariance = u_operation_inputs[12];

  int p13= (int)u_operation_inputs[13];
  int p14= ROI_op_inst->print_frequency	= (int)u_operation_inputs[14];

  ROI_op_inst->grid_resolution = 50*pow(2,p13) + 1;

free(u_operation_inputs);

// Set additional Global Variables
g_nsv = 2*(NOTHERS + g_nsd + g_nld + g_nhda + g_nhdm);
g_mbacteria = fmax(g_nsd,g_nld);
g_mbacteria = fmax(g_mbacteria,g_nhda);
g_mbacteria = fmax(g_mbacteria,g_nhdm);

g_fr_mpd	= 0.001*p5/(p2*p2*0.25*PI);
g_vcolon_l	= p1*p2*p2*0.25*PI*1000;
g_vproximal_l	= PT*g_vcolon_l;
g_vtransverse_l	= (TD-PT)*g_vcolon_l;
g_vdistal_l	= (1-TD)*g_vcolon_l;
g_vsintest_l	= (p4*p4*0.25)*PI*p3*1000;

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double NormRand(double std_dev, double mean, int seed)
{

int counter = 0;

int random_variable1;
int random_variable2;
double random_percentage1;
double random_percentage2;

srand(rand());
random_variable1 = rand();
random_variable2 = rand();

double r_var1 = random_variable1;
double r_var2 = random_variable2;
double r_max = RAND_MAX;

random_percentage1 = r_var1/r_max;
random_percentage2 = r_var2/r_max;

double u1,u2;
double z0,z1;

u1 = random_percentage1;
u2 = random_percentage2;

z0 = sqrt(-2.0 * log(u1))*cos(2.0*PI*u2);
z1 = sqrt(-2.0 * log(u1))*sin(2.0*PI*u2);

//printf("%f %f\n",z0,z1);

double rfinal = (z0 * std_dev * mean) + mean;

//printf("\nFinal random Number = %f \n", rfinal);


return rfinal;

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int PMG(double** PMat, double* DefaultYParameters)
{
int nsv = g_nsv;
int nSD = g_nsd;
int nLD = g_nld;
int nHDA = g_nhda;
int nHDM = g_nhdm;
int nOther = NOTHERS;
double variance = 0.0; // g_bvariance;
int indexi,indexj;

int seed;
srand(time(NULL));
seed = rand();

// Hydrolysis
PMat[0][0] = DefaultYParameters[0];
PMat[0][nOther-1] = -1;

// Glucose Utilization
PMat[1][0] = -1;
PMat[1][1] = DefaultYParameters[1];
PMat[1][2] = DefaultYParameters[2];
PMat[1][3] = DefaultYParameters[3];
PMat[1][4] = DefaultYParameters[4];
PMat[1][5] = DefaultYParameters[5];
PMat[1][7] = DefaultYParameters[6];
PMat[1][8] = DefaultYParameters[7];

for(indexi=0; indexi< nSD; indexi++){
PMat[1][nOther+indexi] = NormRand(variance,DefaultYParameters[20],seed);
}

//Lactate Utilization
PMat[2][1] = -1;
PMat[2][2] = DefaultYParameters[8];
PMat[2][3] = DefaultYParameters[9];
PMat[2][4] = DefaultYParameters[10];
PMat[2][5] = DefaultYParameters[11];
PMat[2][7] = DefaultYParameters[12];
PMat[2][8] = DefaultYParameters[13];

for(indexi=0;indexi<nLD;indexi++){
PMat[2][nOther+nSD+indexi] = NormRand(variance,DefaultYParameters[21],seed);
}


//Acetogenesis
PMat[3][2] = -1;
PMat[3][3] = DefaultYParameters[14];
PMat[3][7] = DefaultYParameters[15];
PMat[3][8] = DefaultYParameters[16];

for(indexi=0;indexi<nHDA;indexi++){
PMat[3][nOther+nSD+nLD+indexi] = \
  NormRand(variance,DefaultYParameters[22],seed);
}


//Methanogenesis
PMat[4][2] = -1;
PMat[4][6] = DefaultYParameters[17];
PMat[4][7] = DefaultYParameters[18];
PMat[4][8] = DefaultYParameters[19];

for(indexi=0;indexi<nHDM;indexi++){
PMat[4][nOther+nSD+nLD+nHDA+indexi] = \
  NormRand(variance,DefaultYParameters[23],seed);
}


//SD Bacterial Death
for(indexi=0;indexi<nSD;indexi++){
	PMat[5][nOther+indexi] = NormRand(variance,-1,seed);
}
//LD Bacterial Death
for(indexi=0;indexi<nLD;indexi++){
	PMat[6][nOther+nSD+indexi] = NormRand(variance,-1,seed);}
//HDA Bacterial Death
for(indexi=0;indexi<nHDA;indexi++){
	PMat[7][nOther+nSD+nLD+indexi] = NormRand(variance,-1,seed);}
//HDM Bacterial Death
for(indexi=0;indexi<nHDM;indexi++){
	PMat[8][nOther+nSD+nLD+nHDA+indexi] = NormRand(variance,-1,seed);}
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int KMG(double** KRateMatrix, double** KHSConstant, 
  double* KRateParameters_vector)
{
int seed;
srand(time(NULL));
seed = rand();

int nsv = g_nsv;
int nSD = g_nsd;
int nLD = g_nld;
int nHDA = g_nhda;
int nHDM = g_nhdm;
int nOther = NOTHERS;
double variance = g_bvariance;

int indexi,indexj;
for (indexj=0;indexj<nSD;indexj++){
  KRateMatrix[0][indexj] = NormRand(variance,KRateParameters_vector[0],seed);
  KRateMatrix[1][indexj] = NormRand(variance,KRateParameters_vector[1],seed);
  KRateMatrix[5][indexj] = NormRand(variance,KRateParameters_vector[5],seed);

  KHSConstant[0][indexj] = NormRand(variance,KRateParameters_vector[9],seed);
  KHSConstant[1][indexj] = NormRand(variance,KRateParameters_vector[10],seed);
}

for (indexj = 0;indexj<nLD;indexj++){
  KRateMatrix[2][indexj] = NormRand(variance,KRateParameters_vector[2],seed);
  KRateMatrix[6][indexj] = NormRand(variance,KRateParameters_vector[6],seed);

  KHSConstant[2][indexj] = NormRand(variance,KRateParameters_vector[11],seed);
}

for (indexj = 0;indexj<nHDA;indexj++){
  KRateMatrix[3][indexj] = NormRand(variance,KRateParameters_vector[3],seed);
  KRateMatrix[7][indexj] = NormRand(variance,KRateParameters_vector[7],seed);

  KHSConstant[3][indexj] = NormRand(variance,KRateParameters_vector[12],seed);
}

for (indexj = 0;indexj<nHDM;indexj++){
  KRateMatrix[4][indexj] = NormRand(variance,KRateParameters_vector[4],seed);
  KRateMatrix[8][indexj] = NormRand(variance,KRateParameters_vector[8],seed);

  KHSConstant[4][indexj] = NormRand(variance,KRateParameters_vector[13],seed);
}

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int TMG(double** TransportParameters, double* Parameters)
{
int indexi,indexj,indexk;

int nsv = g_nsv;
int nSD = g_nsd;
int nLD = g_nld;
int nHDA = g_nhda;
int nHDM = g_nhdm;
int nOther = NOTHERS;
double variance = g_pvariance;

indexk = 0;
for(indexi=0;indexi<4;indexi++){

  int seed;
  srand(time(NULL));
  seed = rand();

  for(indexj=0;indexj<nOther;indexj++){
  TransportParameters[indexi][indexj] = Parameters[45+indexk];
  indexk=indexk+1;
  }

  for(indexj=nOther;indexj<nOther+nSD;indexj++){
  TransportParameters[indexi][indexj] = \
    NormRand(variance,Parameters[45+indexk],seed);}

  indexk=indexk+1;

  for(indexj=nOther+nSD;indexj<nOther+nSD+nLD;indexj++){
  TransportParameters[indexi][indexj] = \
    NormRand(variance,Parameters[45+indexk],seed);}

  indexk=indexk+1;

  for(indexj=nOther+nSD+nLD;indexj<nOther+nSD+nLD+nHDA;indexj++){
  TransportParameters[indexi][indexj] = \
    NormRand(variance,Parameters[45+indexk],seed);}

  indexk=indexk+1;

  for(indexj=nOther+nSD+nLD+nHDA;indexj<nOther+nSD+nLD+nHDA+nHDM;indexj++){
  TransportParameters[indexi][indexj] = \
    NormRand(variance,Parameters[45+indexk],seed);}

  indexk=indexk+1;
}

return(0);
} // end function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int ReadReactionParameters(void* user_data)
{
// RRP scope variables:
	int indexi,indexj,indexk;

	ReactionData RRP_reac_data;
	RRP_reac_data = (ReactionData)user_data;

//------------------------------------------------------------------------------
// Read in parameters from text file:
	int max_reac_parameters = 200;
	double *parameters;
	const int bsz = 300;
	char buf[bsz];

	parameters =malloc(max_reac_parameters * sizeof (double));

	FILE *fparameters;
	fparameters = fopen("InputFiles/ADM1Parameters_compuGUT.txt","r");

	if (fparameters != NULL){
		indexi = 0;
	    while (fgets(buf, bsz, fparameters) != NULL) {
	      if (buf[0] != '#' && buf[0]!= '\n') {
		  sscanf(buf,"%lf", &parameters[indexi]);
		  indexi++;
	      }
	    }
	} else {
		  printf("\nReaction Parameters\n");
	    perror("fopen");
	}
	fclose(fparameters);

//------------------------------------------------------------------------------	
// Peterson Matrix Generation
  int num_YieldParameters = 24;
  double *PMatVector;
  double **PetersonMatrix_YC;
  PMatVector = malloc(num_YieldParameters * sizeof (double));
  PetersonMatrix_YC = Make2DDoubleArray(NAP,g_nsv/2);

  for (indexi=0;indexi<num_YieldParameters;indexi++){
    PMatVector[indexi] = parameters[7+indexi];
  }

  PMG(PetersonMatrix_YC, PMatVector);

  for (indexi=0;indexi<NAP;indexi++){
    for (indexj=0; indexj<g_nsv/2; indexj++){
     RRP_reac_data->PMat[indexi][indexj] = PetersonMatrix_YC[indexi][indexj];
    }
  }

  time_t current_time;
  char* string_ctime;
  current_time = time(NULL);
  string_ctime = ctime(&current_time);

  FILE *f3;
  f3 = fopen("SimulationParameters.txt","w");
  fprintf(f3,"compuGUT_fm Simulation Parameters: \n");
  fprintf(f3,"Date and time: %s\n\n",string_ctime);
  fprintf(f3,"Peterson Matrix:\n");
  fprintf(f3,"-------------------\n\n");

    for(indexi=0;indexi<NAP;indexi++){
      for(indexj=0;indexj<g_nsv/2;indexj++){
	fprintf(f3,"%f ",PetersonMatrix_YC[indexi][indexj]);
      }
      fprintf(f3,"\n");
    }


  free(PMatVector);
    for(indexi=0;indexi<NAP;indexi++){
      free(PetersonMatrix_YC[indexi]);
    }
    free(PetersonMatrix_YC);
//------------------------------------------------------------------------------	
// Kinetic Rates Matrix Generation
  int num_kinetic_parameters = 14;
  double *KRateParameters_vector;
  double **KRateMatrix;
  double **KHSConstant;
  KRateParameters_vector = malloc(num_kinetic_parameters * sizeof (double));
  KRateMatrix = Make2DDoubleArray(NAP,g_mbacteria);
  KHSConstant = Make2DDoubleArray(NAP,g_mbacteria);

  for(indexi=0;indexi<num_kinetic_parameters;indexi++){
    KRateParameters_vector[indexi] = parameters[31+indexi];
  }

  KMG(KRateMatrix,KHSConstant,KRateParameters_vector);

  for (indexi=0;indexi<NAP;indexi++){
    for (indexj=0;indexj<g_mbacteria;indexj++){
      RRP_reac_data->KRMat[indexi][indexj]=KRateMatrix[indexi][indexj];
    }
  }

  for (indexi=0;indexi<5;indexi++){
    for (indexj=0;indexj<g_mbacteria;indexj++){
      RRP_reac_data->KHSMat[indexi][indexj]= \
                     KHSConstant[indexi][indexj];
    }
  }

  fprintf(f3,"\n\n\n");
  fprintf(f3,"Kinetic Rates Matrix: \n");
  fprintf(f3,"---------------------\n\n");
  for (indexi=0;indexi<NAP;indexi++){
    for (indexj=0;indexj<g_mbacteria;indexj++){
      fprintf(f3,"%f ",KRateMatrix[indexi][indexj]);
    }
  fprintf(f3,"\n");
  }
  fprintf(f3,"\n\n\n");
  fprintf(f3,"Kinetic Half Saturation Matrix: \n");
  fprintf(f3,"---------------------\n\n");
  for (indexi=0;indexi<5;indexi++){
    for (indexj=0;indexj<g_mbacteria;indexj++){
      fprintf(f3,"%f ",KHSConstant[indexi][indexj]);
    }
    fprintf(f3,"\n");
  }


  free(KRateParameters_vector);

  for(indexi=0;indexi<NAP;indexi++){
    free(KRateMatrix[indexi]);
    free(KHSConstant[indexi]);
  }
  free(KRateMatrix);
  free(KHSConstant);
//------------------------------------------------------------------------------	
// Transport Rates Matrix Generation (Donor controlled transport)
  int num_transport_parameters = g_nsv/2;
  int num_colon_locations = 4;// 3 colon location + 1 for reverse transport
  double** TransportParameters;

  TransportParameters = \
  Make2DDoubleArray(num_colon_locations,num_transport_parameters);
  TMG(TransportParameters,parameters);

  for (indexi=0;indexi<num_colon_locations;indexi++){
    for (indexj=0;indexj<num_transport_parameters;indexj++){
      RRP_reac_data->TMat[indexi][indexj]=TransportParameters[indexi][indexj];
    }
  }

  fprintf(f3,"\n\n\n");
  fprintf(f3,"Transport Parameters Matrix: \n");
  fprintf(f3,"---------------------\n\n");
  for(indexi=0;indexi<num_colon_locations;indexi++){
    for(indexj=0;indexj<num_transport_parameters;indexj++){
   fprintf(f3,"%d %d %f \n", indexi,indexj,TransportParameters[indexi][indexj]);
    }
    fprintf(f3,"\n");
  }
  fprintf(f3,"\n\nContact Arun [amoorthy@uoguelph.ca] or Hermann");
  fprintf(f3," [heberl@uoguelph.ca] with questions or feedback");
  fprintf(f3,"\n\n(c) Moorthy & Eberl, 2015\n\n");
  fclose(f3);

  for(indexi=0;indexi<num_colon_locations;indexi++){
    free(TransportParameters[indexi]);
  }
  free(TransportParameters);
//------------------------------------------------------------------------------	
  RRP_reac_data->pH_v[0] = parameters[101];
  RRP_reac_data->pH_v[1] = parameters[102];
  RRP_reac_data->pH_v[2] = parameters[103];


  free(parameters);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int InitialConditions(double** PSolution, int ngrids)
{
	int indexi,indexj,indexk;

	int nsv = g_nsv;
	int nSD = g_nsd;
	int nLD = g_nld;
	int nHDA = g_nhda;
	int nHDM = g_nhdm;
	int nOther = NOTHERS;

	//Load Parameters
	int 	maxNumParameters = 150;
	double *Parameters;;
	const int bsz = 300;
	char *buf;

	Parameters = malloc(maxNumParameters * sizeof (double));
	buf = malloc( bsz * sizeof (char));

	FILE *fparameters;
	fparameters = fopen("InputFiles/user_colon_initializer.txt","r");

		int run_index = 0;
		int num_lines_parameter_file = 9; // volatile code.
		indexi = 0;

		while(indexi<num_lines_parameter_file){
			fgets(buf, bsz, fparameters);

			if (buf[0] != '#' && buf[0]!= '\n') {
			sscanf(buf,"%lf", &Parameters[run_index]);
			run_index++;
			}
			indexi++;
		}

	fclose(fparameters);

	/* for(indexi=0;indexi<5;indexi++){
		//~ printf("%d %f\n",indexi,Parameters[indexi]);
	 }
	PauseMessage(912);*/

	double InitialMucus_gpL = Parameters[0];
	double InitialSDBiomass_gpL = Parameters[1];
	double InitialLDBiomass_gpL = Parameters[2];
	double InitialHDABiomass_gpL = Parameters[3];
	double InitialHDMBiomass_gpL = Parameters[4];
	double TotalBiomass = nSD + nLD + nHDA + nHDM;

	free(Parameters);
	free(buf);

// Assign initial Mucus density along length of colon
for (indexi=0; indexi < ngrids; indexi++){
	PSolution[indexi][19] = InitialMucus_gpL;
}
// Assign initial SD density along length of colon
if (nSD > 0){
	for (indexi=0;indexi<ngrids;indexi++){
	for(indexj = 2*nOther+1; indexj < 2*(nOther+nSD); indexj+=2){
		PSolution[indexi][indexj] = InitialSDBiomass_gpL/nSD;
	}
	}
}
// Assign initial LD density along length of colon
if (nLD > 0){
	for(indexi=0;indexi<ngrids;indexi++){
	for(indexj = 2*(nOther+nSD)+1;
		 indexj < 2*(nOther+nSD+nLD);
		 indexj+=2){
		PSolution[indexi][indexj] = InitialLDBiomass_gpL/nLD;
	}
	}
}
// Assign initial HDA density along length of colon
if (nHDA > 0){
	for(indexi=0;indexi<ngrids;indexi++){
	for(indexj = 2*(nOther+nSD+nLD)+1;
		 indexj < 2*(nOther+nSD+nLD+nHDA);
		 indexj+=2){
		 PSolution[indexi][indexj] = InitialHDABiomass_gpL/nHDA;
	}
	}
}
// Assign initial HDM density along length of colon
if (nHDM > 0){
	for(indexi=0;indexi<ngrids;indexi++){
	for(indexj = 2*(nOther+nSD+nLD+nHDA)+1;
		indexj < 2*(nOther+nSD+nLD+nHDA+nHDM);
		indexj+=2){
		PSolution[indexi][indexj] = InitialHDMBiomass_gpL/nHDM;
	}
	}
}

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int PrintFullSolution(double** solution, int ngrids, char a[], int fileiter,
  double t)
{
	/* Description: prints full solution */

	int indexi,indexj,indexk;
	int nsv = g_nsv;

	char	  filename[64];
	char	  ofilename[64];

	FILE *f1;
	sprintf(filename, "%c%c%c%c%c%d.txt",a[0],a[1],a[2],a[3],a[4],fileiter);
	f1 = fopen(filename, "w");
	fprintf(f1,"#Time: %f\n#\n#\n",t);
	for (indexi=0;indexi<ngrids;indexi++){
		for (indexj=0;indexj<nsv;indexj++){
		fprintf(f1,"%f\t",solution[indexi][indexj]);
		}
		fprintf(f1,"\n");
	}
	fclose(f1);



return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int PrintCondensedSolution(double** solution, int ngrids, char a[],
  int fileiter, double t, int increment)
{
	/* Description: prints condensed solution */

	int indexi,indexj,indexk;
	int nsv = g_nsv;

	char	  filename[64];
	char	  ofilename[64];

	FILE *f1;
  sprintf(filename, "%c%c%c%c%c%d.txt",a[0],a[1],a[2],a[3],a[4],fileiter);
	f1 = fopen(filename, "w");
	fprintf(f1,"#Time: %f\n#\n#\n",t);

	indexk = 0;
	for (indexi=0;indexi<MinR;indexi++){
		for (indexj=0;indexj<nsv;indexj++){
		fprintf(f1,"%f\t",solution[indexk][indexj]);
		}
		fprintf(f1,"\n");
		indexk=indexk+increment;
	}
	fclose(f1);



return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double SplineInterp(double x, double Pp, double Pt, double Pd)
{
	double q,x1,x2,y1,y2,k1,k2,t,a,b;

	if (x >=0 && x < PT - SpA){
		q = Pp;}
	else if (x>=PT-SpA && x <PT+SpA){
		x1 = PT-SpA;
		x2 = PT+SpA;
		y1 = Pp;
		y2 = Pt;
		k1 = 0;
		k2 = 0;

		t = (x-x1)/(x2-x1);
		a = k1*(x2-x1) - (y2 - y1);
		b = -k2*(x2 - x1) + (y2 - y1);

		q = (1-t)*y1 + t*y2 + t*(1-t)*(a*(1-t)+b*t);}

	else if (x > PT+SpA && x < TD-SpA){
		q = Pt;}

	else if (x>= TD-SpA && x <=TD+SpA){
        	x1 = TD-SpA;
        	x2 = TD+SpA;
        	y1 = Pt;
        	y2 = Pd;
        	k1 = 0;
        	k2 = 0;

        	t = (x-x1)/(x2-x1);
        	a = k1*(x2-x1) - (y2 - y1);
        	b = -k2*(x2 - x1) + (y2 - y1);

        	q = (1-t)*y1 + t*y2 + t*(1-t)*(a*(1-t)+b*t);}

	else {q = Pd;}

return q;
} //end function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double fo_kinetics(double mu, double A)
{

double C;
double D;

C = mu*A;
D = fmax(0,C);

return(D);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double fo_reaction_kinetics(double mu, double A, double B)
{

double C;
double D;

C = mu*A*B;
D = fmax(0,C);

return(D);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double contois_kinetics(double mu, double K, double A, double B)
{

double C;
double D;

C = mu*A*B/(K*B + A + EPSI);

D = fmax(0,C);

return(D);
//return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double contois_kinetics_2(double A, double Bt, double Bb)
{
	double C;
	double D;

	C = A*Bt/(Bb+A + EPSI);
	D = fmax(0,C);

return(D);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double monod_kinetics(double mu, double K, double A, double B)
{

double C;
double D;

C = mu*A*B/(K + A + EPSI);

D = fmax(0,C);

return(D);
//return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int reac(double* input, double* output, void* user_reac_data)
{
//~ % 0 - 0,1 	- monomer sugar
//~ % 1 - 2,3 	- lactate
//~ % 2 - 4,5 	- hydrogen
//~ % 3 - 6,7 	- acetate
//~ % 4 - 8,9 	- proprionate
//~ % 5 - 10,11 - butyrate
//~ % 6 - 12,13 - methane
//~ % 7 - 14,15 - carbondioxide
//~ % 8 - 16,17 - water
//~ % 9 - 18,19 - polymer sugar (fiber/mucin)
//~ %
//~ % followed by:
//~ % Groups of Bacteria [odd numbers for mucus, even for lumen]

	int indexi,indexj,indexk;
	int rflag;
	int nsv = g_nsv;
	int pnsv = g_nsv/2;
	int MaxBacteria = g_mbacteria;

	int nSD = g_nsd;
	int nLD = g_nld;
	int nHDA = g_nhda;
	int nHDM = g_nhdm;
	int nOther = NOTHERS;
	double variance = g_bvariance;

	ReactionData reac_reac_data;
	reac_reac_data=(ReactionData)user_reac_data;

	double pH_value;
	double IpH;
	double IpH_power;
	double LocX;

	LocX = reac_reac_data->LocX;

	double* pH_vector = malloc(3 * sizeof(double));
	rflag = Clean1DArray(pH_vector,3);

		for(indexi=0;indexi<3;indexi++){
		pH_vector[indexi] = reac_reac_data->pH_v[indexi];
		}

	pH_value = SplineInterp(LocX,pH_vector[0],pH_vector[1],pH_vector[2]);

	IpH_power = -3.0*(pH_value - 6.69)/((6.69-5.8)*(6.69*5.8)); // from ADM1
	IpH = exp(IpH_power);

	double** PetersonMatrix_reac;
	double** KRateMatrix_reac;
	double** KHSConstant_reac;
	double** TransMatrix_reac;

	PetersonMatrix_reac = Make2DDoubleArray(NAP,pnsv);
	rflag = Clean2DArray(PetersonMatrix_reac,NAP,pnsv);
	KRateMatrix_reac = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(KRateMatrix_reac,NAP,MaxBacteria);
	KHSConstant_reac = Make2DDoubleArray(5,MaxBacteria);
	rflag = Clean2DArray(KHSConstant_reac,5,MaxBacteria);
	TransMatrix_reac = Make2DDoubleArray(4,pnsv);

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  PetersonMatrix_reac[indexi][indexj] = \
		  reac_reac_data->PMat[indexi][indexj];
	  }
    }

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		 KRateMatrix_reac[indexi][indexj] = \
		 reac_reac_data->KRMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<5;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		  KHSConstant_reac[indexi][indexj] = \
		  reac_reac_data->KHSMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<4;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  TransMatrix_reac[indexi][indexj] = \
		  reac_reac_data->TMat[indexi][indexj];
	  }
    }
//------------------------------------------------------------------------------	
// Set up problem---------------------------------------------------------------	

	double* Lin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lin,pnsv);

	double* Min = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Min,pnsv);

	double* Lout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lout,pnsv);

	double* Mout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Mout,pnsv);

	double* ReacL = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacL,pnsv);

	double* ReacM = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacM,pnsv);

	double temp=0.0;


indexj=0;
for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		Lin[indexj] = input[indexi];
	} else if(indexi%2==1) {
		Min[indexj] = input[indexi];
		indexj=indexj+1;
	}
}

	double** rhoL;
	rhoL = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoL,NAP,MaxBacteria);

	double** rhoM;
	rhoM = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoM,NAP,MaxBacteria);

	//~ for(indexi=0;indexi<NAP;indexi++){
	//~ for(indexj=0;indexj<MaxBacteria;indexj++){
		//~ printf("%f ",rhoL[indexi][indexj]);
	//~ }
	//~ printf("\n");
//~ }
//~ printf("\n");
//~ for(indexi=0;indexi<NAP;indexi++){
	//~ for(indexj=0;indexj<MaxBacteria;indexj++){
		//~ printf("%f ",rhoM[indexi][indexj]);
	//~ }
	//~ printf("\n");
//~ }
//~ printf("\n");
//------------------------------------------------------------------------------
// Determine Biochemical reaction Rates-----------------------------------------

/*Hydrolysis (rho1);*/

//~ /*Contois 1 - OLD - Reminder of challenges of contois kinetics and comp.*/
//~ if (nSD > 0) {
  //~ for(indexi=0;indexi<nSD;indexi++){
  //~ rhoL[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
                                     //~ KHSConstant_reac[0][indexi], \
                                     //~ Lin[9], \
                                     //~ Lin[nOther+indexi]);
                                     //~
	//~ rhoM[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
	                                   //~ KHSConstant_reac[0][indexi], \
	                                   //~ Min[9], \
	                                   //~ Min[nOther+indexi]);
  //~ }
//~ }

if (g_hydro_scheme == 2){
  /*First-order disintegration*/
  double cum_mu = 0.0;
  double mean_mu = 0.0;

  for(indexi=0;indexi<nSD;indexi++){
   cum_mu = cum_mu + KRateMatrix_reac[0][indexi];
  }
  mean_mu = cum_mu/nSD;
  rhoL[0][0] = fo_kinetics(mean_mu,Lin[9]);
  rhoM[0][0] = fo_kinetics(mean_mu,Min[9]);

} else if (g_hydro_scheme == 3){
  /*MONOD*/
  double sl_conc=3.0;  // approx. steady state lumen concentration (CSTR)
  double sm_conc=42.0; // approx. steady state mucus concentration (CSTR)

  if (nSD > 0) {
   for(indexi=0;indexi<nSD;indexi++){
    rhoL[0][indexi] = monod_kinetics(KRateMatrix_reac[0][indexi], \
                                     sl_conc*KHSConstant_reac[0][indexi], \
			             Lin[9], \
                                     Lin[nOther+indexi]);

    rhoM[0][indexi] = monod_kinetics(KRateMatrix_reac[0][indexi], \
                                     sm_conc*KHSConstant_reac[0][indexi], \
                                     Min[9], \
                                     Min[nOther+indexi]);
   }
  }
} else {
 /*Contois 2*/
 double LcumGrowth = 0.0;
 double LcumHalfSat = 0.0;
 double McumGrowth = 0.0;
 double McumHalfSat = 0.0;

 if (nSD > 0) {
   for (indexi=0;indexi<nSD;indexi++){
     LcumGrowth = LcumGrowth + KRateMatrix_reac[0][indexi]*Lin[nOther+indexi];
     LcumHalfSat = LcumHalfSat + KHSConstant_reac[0][indexi]*Lin[nOther+indexi];
     McumGrowth = McumGrowth + KRateMatrix_reac[0][indexi]*Min[nOther+indexi];
     McumHalfSat = McumHalfSat + KHSConstant_reac[0][indexi]*Min[nOther+indexi];
   }
 }

 rhoL[0][0] = contois_kinetics_2(Lin[9],LcumGrowth,LcumHalfSat);
 rhoM[0][0] = contois_kinetics_2(Min[9],McumGrowth,McumHalfSat);
}


/*Glucose Utilization (rho2) and Death of Sugar Degrading Bacteria (rho6):*/
if (nSD > 0) {
  for(indexi=0;indexi<nSD;indexi++){
	rhoL[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Lin[0], \
	                                 Lin[nOther+indexi]);

	rhoM[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Min[0], \
	                                 Min[nOther+indexi]);

	rhoL[5][indexi] = KRateMatrix_reac[5][indexi]*Lin[nOther+indexi];
	rhoM[5][indexi] = KRateMatrix_reac[5][indexi]*Min[nOther+indexi];
  }
}

/*Lactate Utilization (rho3) and Death of Lactate Degrading Bacteria (rho7):*/
if (nLD > 0) {
  for (indexi=0;indexi<nLD;indexi++){
	rhoL[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Lin[1], \
	                                 Lin[nOther+nSD+indexi]);

	rhoM[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Min[1], \
	                                 Min[nOther+nSD+indexi]);

	rhoL[6][indexi] = KRateMatrix_reac[6][indexi]*Lin[nOther+nSD+indexi];
	rhoM[6][indexi] = KRateMatrix_reac[6][indexi]*Min[nOther+nSD+indexi];
  }
}

/*Homoacetogenesis (rho4) and Death of acetogenic bacteria (rho8):*/
if (nHDA > 0){
  for (indexi=0;indexi<nHDA;indexi++){
    rhoL[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
                                     KHSConstant_reac[3][indexi], \
                                     Lin[2], \
	                             Lin[nOther+nSD+nLD+indexi]);

    rhoM[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
	                             KHSConstant_reac[3][indexi], \
	                             Min[2], \
	                             Min[nOther+nSD+nLD+indexi]);

    rhoL[7][indexi] = KRateMatrix_reac[7][indexi]*Lin[nOther+nSD+nLD+indexi];
    rhoM[7][indexi] = KRateMatrix_reac[7][indexi]*Min[nOther+nSD+nLD+indexi];
  }
}

/*Homomethanogenesis (rho5) and Death of methanogenic bacteria (rho8):*/
if (nHDM > 0){
  for (indexi=0;indexi<nHDM;indexi++){
    rhoL[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                 KHSConstant_reac[4][indexi], \
	                                 Lin[2], \
	                                 Lin[nOther+nSD+nLD+nHDA+indexi]);

    rhoM[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                 KHSConstant_reac[4][indexi], \
	                                 Min[2], \
	                                 Min[nOther+nSD+nLD+nHDA+indexi]);

  rhoL[8][indexi] = KRateMatrix_reac[8][indexi]*Lin[nOther+nSD+nLD+nHDA+indexi];
  rhoM[8][indexi] = KRateMatrix_reac[8][indexi]*Min[nOther+nSD+nLD+nHDA+indexi];
  }
}

//~ for(indexi=0;indexi<NAP;indexi++){
	//~ for(indexj=0;indexj<MaxBacteria;indexj++){
		//~ printf("%f ",rhoL[indexi][indexj]);
	//~ }
	//~ printf("\n");
//~ }
//~ printf("\n");
//~ for(indexi=0;indexi<NAP;indexi++){
	//~ for(indexj=0;indexj<MaxBacteria;indexj++){
		//~ printf("%f ",rhoM[indexi][indexj]);
	//~ }
	//~ printf("\n");
//~ }
//~
//~ PauseMessage(1766);
//------------------------------------------------------------------------------
// Build Differential Equations-------------------------------------------------
double V = g_vcolon_l;
double Vd = 0.9*V;
double Vm = 0.1*V;
double MuMax = 50; // [g/L]
double MuPro = 500.0; // [g/Ld]
double MuGen = 0.0;
double InterpTransport= 0.0;

MuGen = (1.0 - Min[9]/MuMax)*MuPro;

  /*Non-biomass terms*/
  for(indexi=0;indexi<nOther;indexi++){
    ReacL[indexi]=0;
    ReacM[indexi]=0;

    for(indexj=0;indexj<9;indexj++){
      for(indexk=0;indexk<MaxBacteria;indexk++){
        ReacL[indexi] = ReacL[indexi] +  \
	              PetersonMatrix_reac[indexj][indexi]*rhoL[indexj][indexk];
	ReacM[indexi] = ReacM[indexi] + \
	              PetersonMatrix_reac[indexj][indexi]*rhoM[indexj][indexk];
      }
    }

    if(indexi==0){
        InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = ReacL[indexi] - \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vd;

	Mout[indexi] = ReacM[indexi] + \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vm;

     }else if(indexi==nOther-1){
	Lout[indexi] = ReacL[indexi] + \
	               (Vm/Vd)*TransMatrix_reac[3][indexi]*Min[indexi] ;

	Mout[indexi] = MuGen + ReacM[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];

      }else{
	InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = ReacL[indexi] - \
	               InterpTransport*Lin[indexi];

	Mout[indexi] = ReacM[indexi] + \
	               (Vd/Vm)*InterpTransport*Lin[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];
      }

  }

  /*Biomass Terms*/
  for(indexi=0;indexi<nSD;indexi++){
	ReacL[nOther+indexi] = 0;
	ReacM[nOther+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+indexi] = ReacL[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+indexi] = ReacM[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoM[indexj][indexi];
	}

  InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+indexi], \
  TransMatrix_reac[1][nOther+indexi],TransMatrix_reac[2][nOther+indexi]);

	Lout[nOther+indexi] = ReacL[nOther+indexi] - \
	InterpTransport*Lin[nOther+indexi] + \
	(Vm/Vd)*TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];

	Mout[nOther+indexi] = ReacM[nOther+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+indexi] - \
	TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];
  }

  for(indexi=0;indexi<nLD;indexi++){
	ReacL[nOther+nSD+indexi] = 0;
	ReacM[nOther+nSD+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
          ReacL[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] + \
	  PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoL[indexj][indexi];

	  ReacM[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	  PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoM[indexj][indexi];
	}

InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+nSD+indexi],\
		          TransMatrix_reac[1][nOther+nSD+indexi],\
                  TransMatrix_reac[2][nOther+nSD+indexi]);

	Lout[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] - \
	InterpTransport*Lin[nOther+nSD+indexi] + \
	(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi];

	Mout[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+indexi] - \
	TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi];
  }

  for(indexi=0;indexi<nHDA;indexi++){
    ReacL[nOther+nSD+nLD+indexi] = 0;
    ReacM[nOther+nSD+nLD+indexi] = 0;
    for(indexj=0;indexj<9;indexj++){
      ReacL[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] + \
      PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoL[indexj][indexi];

      ReacM[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
      PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoM[indexj][indexi];
    }
      InterpTransport = SplineInterp(LocX, \
                        TransMatrix_reac[0][nOther+nSD+nLD+indexi], \
	                TransMatrix_reac[1][nOther+nSD+nLD+indexi], \
	                TransMatrix_reac[2][nOther+nSD+nLD+indexi]);

  Lout[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] - \
                              InterpTransport*Lin[nOther+nSD+nLD+indexi] + \
  (Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi];

  Mout[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
  (Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+indexi] - \
  TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi];
  }

  for(indexi=0;indexi<nHDM;indexi++){
    ReacL[nOther+nSD+nLD+nHDA+indexi] = 0;
    ReacM[nOther+nSD+nLD+nHDA+indexi] = 0;
    for(indexj=0;indexj<9;indexj++){
      ReacL[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] +\
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoL[indexj][indexi];

      ReacM[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoM[indexj][indexi];
    }

	InterpTransport = SplineInterp(LocX,\
	TransMatrix_reac[0][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[1][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[2][nOther+nSD+nLD+nHDA+indexi]);

	Lout[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] - \
	InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] + \
	(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi];

	Mout[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] - \
	TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi];
  }

  /*Return Completed Vector*/
  indexj=0;
  for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		output[indexi] = Lout[indexj];
	} else if(indexi%2==1){
		output[indexi] = Mout[indexj];
		indexj=indexj+1;
	}
  }

//------------------------------------------------------------------------------
// Free memory:
	free(pH_vector);

	for(indexi=0;indexi<NAP;indexi++){
		free(PetersonMatrix_reac[indexi]);
		free(KRateMatrix_reac[indexi]);
		free(rhoL[indexi]);
		free(rhoM[indexi]);
	}
	free(PetersonMatrix_reac);
	free(KRateMatrix_reac);
	free(rhoL);
	free(rhoM);

	for(indexi=0;indexi<5;indexi++){
		free(KHSConstant_reac[indexi]);
	}
	free(KHSConstant_reac);

	for(indexi=0;indexi<4;indexi++){
		free(TransMatrix_reac[indexi]);
	}
	free(TransMatrix_reac);

	free(Lin);
	free(Min);
	free(Lout);
	free(Mout);
	free(ReacL);
	free(ReacM);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int reac3sImplicit(double* input, double* pinput, double* output,
  void* user_reac_data, double V1, double V2)
{
//~ % 0 - 0,1 	- monomer sugar
//~ % 1 - 2,3 	- lactate
//~ % 2 - 4,5 	- hydrogen
//~ % 3 - 6,7 	- acetate
//~ % 4 - 8,9 	- proprionate
//~ % 5 - 10,11 - butyrate
//~ % 6 - 12,13 - methane
//~ % 7 - 14,15 - carbondioxide
//~ % 8 - 16,17 - water
//~ % 9 - 18,19 - polymer sugar (fiber/mucin)
//~ %
//~ % followed by:
//~ % Groups of Bacteria [odd numbers for mucus, even for lumen]

	int indexi,indexj,indexk;
	int rflag;
	int nsv = g_nsv;
	int pnsv = g_nsv/2;
	int MaxBacteria = g_mbacteria;

	int nSD = g_nsd;
	int nLD = g_nld;
	int nHDA = g_nhda;
	int nHDM = g_nhdm;
	int nOther = NOTHERS;
	double variance = g_bvariance;

	ReactionData reac_reac_data;
	reac_reac_data=(ReactionData)user_reac_data;

	double pH_value;
	double IpH;
	double IpH_power;
	double LocX;

	LocX = reac_reac_data->LocX;


	double* pH_vector = malloc(3 * sizeof(double));
	rflag = Clean1DArray(pH_vector,3);

		for(indexi=0;indexi<3;indexi++){
		pH_vector[indexi] = reac_reac_data->pH_v[indexi];
		}

	pH_value = SplineInterp(LocX,pH_vector[0],pH_vector[1],pH_vector[2]);

	IpH_power = -3.0*(pH_value - 6.69)/((6.69-5.8)*(6.69*5.8)); // from ADM1
	IpH = exp(IpH_power);

	double** PetersonMatrix_reac;
	double** KRateMatrix_reac;
	double** KHSConstant_reac;
	double** TransMatrix_reac;

	PetersonMatrix_reac = Make2DDoubleArray(NAP,pnsv);
	rflag = Clean2DArray(PetersonMatrix_reac,NAP,pnsv);
	KRateMatrix_reac = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(KRateMatrix_reac,NAP,MaxBacteria);
	KHSConstant_reac = Make2DDoubleArray(5,MaxBacteria);
	rflag = Clean2DArray(KHSConstant_reac,5,MaxBacteria);
	TransMatrix_reac = Make2DDoubleArray(4,pnsv);

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  PetersonMatrix_reac[indexi][indexj] = \
		  reac_reac_data->PMat[indexi][indexj];
	  }
    }

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		 KRateMatrix_reac[indexi][indexj] = \
		 reac_reac_data->KRMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<5;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		  KHSConstant_reac[indexi][indexj] = \
		  reac_reac_data->KHSMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<4;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  TransMatrix_reac[indexi][indexj] = \
		  reac_reac_data->TMat[indexi][indexj];
	  }
    }
//------------------------------------------------------------------------------
// Set up problem---------------------------------------------------------------

	double* Lin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lin,pnsv);

	double* Lpin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lpin,pnsv);

	double* Min = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Min,pnsv);

	double* Mpin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Mpin,pnsv);

	double* Lout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lout,pnsv);

	double* Mout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Mout,pnsv);

	double* ReacL = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacL,pnsv);

	double* ReacM = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacM,pnsv);

	double temp=0.0;


indexj=0;
for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		Lin[indexj] = input[indexi];
		Lpin[indexj] = pinput[indexi];
	} else if(indexi%2==1) {
		Min[indexj] = input[indexi];
		Mpin[indexj] = pinput[indexi];
		indexj=indexj+1;
	}
}

	double** rhoL;
	rhoL = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoL,NAP,MaxBacteria);

	double** rhoM;
	rhoM = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoM,NAP,MaxBacteria);

//------------------------------------------------------------------------------
// Determine Biochemical reaction Rates-----------------------------------------

/*Hydrolysis (rho1);*/

//~ /*Contois 1 - OLD - Reminder of challenges of contois kinetics and comp.*/
//~ if (nSD > 0) {
  //~ for(indexi=0;indexi<nSD;indexi++){
  //~ rhoL[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
                                     //~ KHSConstant_reac[0][indexi], \
                                     //~ Lin[9], \
                                     //~ Lin[nOther+indexi]);
                                       //~
	//~ rhoM[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
	                                   //~ KHSConstant_reac[0][indexi], \
	                                   //~ Min[9], \
	                                   //~ Min[nOther+indexi]);
  //~ }
//~ }

if (g_hydro_scheme == 2){
  /*First-order disintegration*/

  rhoL[0][0] = fo_kinetics(10.62,Lin[9]);
  rhoM[0][0] = fo_kinetics(10.62,Min[9]);

} else if (g_hydro_scheme == 3){
  /*First-order Reaction*/
  if (nSD > 0) {
   for(indexi=0;indexi<nSD;indexi++){
    rhoL[0][indexi] = fo_reaction_kinetics(KRateMatrix_reac[0][indexi], \
                                           Lin[9], \
                                           Lin[nOther+indexi]);

    rhoM[0][indexi] = fo_reaction_kinetics(KRateMatrix_reac[0][indexi], \
                                           Min[9], \
                                           Min[nOther+indexi]);
   }
  }
} else {
 /*Contois NEW*/
 double LcumGrowth = 0.0;
 double LcumHalfSat = 0.0;
 double McumGrowth = 0.0;
 double McumHalfSat = 0.0;

 if (nSD > 0) {
   for (indexi=0;indexi<nSD;indexi++){
     LcumGrowth = LcumGrowth + KRateMatrix_reac[0][indexi]*Lin[nOther+indexi];
     LcumHalfSat = LcumHalfSat + KHSConstant_reac[0][indexi]*Lin[nOther+indexi];
     McumGrowth = McumGrowth + KRateMatrix_reac[0][indexi]*Min[nOther+indexi];
     McumHalfSat = McumHalfSat + KHSConstant_reac[0][indexi]*Min[nOther+indexi];
   }
 }

 rhoL[0][0] = contois_kinetics_2(Lin[9],LcumGrowth,LcumHalfSat);
 rhoM[0][0] = contois_kinetics_2(Min[9],McumGrowth,McumHalfSat);
}

/*Glucose Utilization (rho2) and Death of Sugar Degrading Bacteria (rho6):*/
if (nSD > 0) {
  for(indexi=0;indexi<nSD;indexi++){
	rhoL[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Lin[0], \
	                                 Lin[nOther+indexi]);

	rhoM[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Min[0], \
	                                 Min[nOther+indexi]);

	rhoL[5][indexi] = KRateMatrix_reac[5][indexi]*Lin[nOther+indexi];
	rhoM[5][indexi] = KRateMatrix_reac[5][indexi]*Min[nOther+indexi];
  }
}

/*Lactate Utilization (rho3) and Death of Lactate Degrading Bacteria (rho7):*/
if (nLD > 0) {
  for (indexi=0;indexi<nLD;indexi++){
	rhoL[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Lin[1], \
	                                 Lin[nOther+nSD+indexi]);

	rhoM[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Min[1], \
	                                 Min[nOther+nSD+indexi]);

	rhoL[6][indexi] = KRateMatrix_reac[6][indexi]*Lin[nOther+nSD+indexi];
	rhoM[6][indexi] = KRateMatrix_reac[6][indexi]*Min[nOther+nSD+indexi];
  }
}

/*Homoacetogenesis (rho4) and Death of acetogenic bacteria (rho8):*/
if (nHDA > 0){
  for (indexi=0;indexi<nHDA;indexi++){
	rhoL[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
	                                 KHSConstant_reac[3][indexi], \
	                                 Lin[2], \
	                                 Lin[nOther+nSD+nLD+indexi]);

	rhoM[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
	                                 KHSConstant_reac[3][indexi], \
	                                 Min[2], \
	                                 Min[nOther+nSD+nLD+indexi]);

    rhoL[7][indexi] = KRateMatrix_reac[7][indexi]*Lin[nOther+nSD+nLD+indexi];
    rhoM[7][indexi] = KRateMatrix_reac[7][indexi]*Min[nOther+nSD+nLD+indexi];
  }
}

/*Homomethanogenesis (rho5) and Death of methanogenic bacteria (rho8):*/
if (nHDM > 0){
  for (indexi=0;indexi<nHDM;indexi++){
	rhoL[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                     KHSConstant_reac[4][indexi], \
	                                     Lin[2], \
	                                     Lin[nOther+nSD+nLD+nHDA+indexi]);

	rhoM[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                     KHSConstant_reac[4][indexi], \
	                                     Min[2], \
	                                     Min[nOther+nSD+nLD+nHDA+indexi]);

  rhoL[8][indexi] = KRateMatrix_reac[8][indexi]*Lin[nOther+nSD+nLD+nHDA+indexi];
  rhoM[8][indexi] = KRateMatrix_reac[8][indexi]*Min[nOther+nSD+nLD+nHDA+indexi];
  }
}
//------------------------------------------------------------------------------
// Build Differential Equations-------------------------------------------------
double V = g_vcolon_l;
double Vd = 0.9*V;
double Vm = 0.1*V;
double MuMax = 50.0; // [g/L]
double MuPro = 500.0; // [g/Ld]
double MuGen = 0.0;
double InterpTransport= 0.0;
double fr = g_fr_lpd;
double vratio = V1/V2;
double drate = fr/V2;

MuGen = (1.0 - Min[9]/MuMax)*MuPro;

  /*Non-biomass terms*/
  for(indexi=0;indexi<nOther;indexi++){
	ReacL[indexi]=0;
	ReacM[indexi]=0;

	for(indexj=0;indexj<9;indexj++){
	  for(indexk=0;indexk<MaxBacteria;indexk++){
	  ReacL[indexi] = ReacL[indexi] +  \
	               PetersonMatrix_reac[indexj][indexi]*rhoL[indexj][indexk];
	  ReacM[indexi] = ReacM[indexi] + \
	               PetersonMatrix_reac[indexj][indexi]*rhoM[indexj][indexk];
	  }
	}

	if(indexi==0){
	InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] - \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vd;

	Mout[indexi] = ReacM[indexi] + \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vm;

	}else if(indexi==nOther-1){
	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] + \
	               (Vm/Vd)*TransMatrix_reac[3][indexi]*Min[indexi] ;

	Mout[indexi] = MuGen + ReacM[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];

	}else{
	InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] - \
	               InterpTransport*Lin[indexi];

	Mout[indexi] = ReacM[indexi] + \
	               (Vd/Vm)*InterpTransport*Lin[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];
	}

  }

  /*Biomass Terms*/
  for(indexi=0;indexi<nSD;indexi++){
	ReacL[nOther+indexi] = 0;
	ReacM[nOther+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+indexi] = ReacL[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+indexi] = ReacM[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoM[indexj][indexi];
	}

  InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+indexi], \
  TransMatrix_reac[1][nOther+indexi],TransMatrix_reac[2][nOther+indexi]);

   Lout[nOther+indexi] = drate*(vratio*Lpin[nOther+indexi]-Lin[nOther+indexi])+\
   ReacL[nOther+indexi] - \
   InterpTransport*Lin[nOther+indexi] + \
   (Vm/Vd)*TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];

   Mout[nOther+indexi] = ReacM[nOther+indexi] + \
   (Vd/Vm)*InterpTransport*Lin[nOther+indexi] - \
   TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];
  }

  for(indexi=0;indexi<nLD;indexi++){
	ReacL[nOther+nSD+indexi] = 0;
	ReacM[nOther+nSD+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
	  ReacL[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] + \
	  PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoL[indexj][indexi];

	  ReacM[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	  PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoM[indexj][indexi];
	}

InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+nSD+indexi],\
	          TransMatrix_reac[1][nOther+nSD+indexi],\
                  TransMatrix_reac[2][nOther+nSD+indexi]);

	Lout[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] - \
	InterpTransport*Lin[nOther+nSD+indexi] + \
       (Vm/Vd)*TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi] + \
	drate*(vratio*Lpin[nOther+nSD+indexi] - Lin[nOther+nSD+indexi]);

	Mout[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+indexi] - \
	TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi];
  }

  for(indexi=0;indexi<nHDA;indexi++){
	ReacL[nOther+nSD+nLD+indexi] = 0;
	ReacM[nOther+nSD+nLD+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoM[indexj][indexi];
	}
	InterpTransport = SplineInterp(LocX, \
	                          TransMatrix_reac[0][nOther+nSD+nLD+indexi], \
	                          TransMatrix_reac[1][nOther+nSD+nLD+indexi], \
	                          TransMatrix_reac[2][nOther+nSD+nLD+indexi]);

Lout[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] - \
                              InterpTransport*Lin[nOther+nSD+nLD+indexi] + \
(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi] +\
drate*(vratio*Lpin[nOther+nSD+nLD+indexi] - Lin[nOther+nSD+nLD+indexi]);

Mout[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
(Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+indexi] - \
TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi];
}

  for(indexi=0;indexi<nHDM;indexi++){
	ReacL[nOther+nSD+nLD+nHDA+indexi] = 0;
	ReacM[nOther+nSD+nLD+nHDA+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
	ReacL[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] +\
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoL[indexj][indexi];

       ReacM[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoM[indexj][indexi];
	}

	InterpTransport = SplineInterp(LocX,\
	TransMatrix_reac[0][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[1][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[2][nOther+nSD+nLD+nHDA+indexi]);

	Lout[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] - \
	InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] + \
	(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi] + \
	drate*(vratio*Lpin[nOther+nSD+nLD+nHDA+indexi] - \
	Lin[nOther+nSD+nLD+nHDA+indexi]);

	Mout[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] - \
	TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi];
  }

  /*Return Completed Vector*/
  indexj=0;
  for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		output[indexi] = Lout[indexj];
	} else if(indexi%2==1){
		output[indexi] = Mout[indexj];
		indexj=indexj+1;
	}
  }

//------------------------------------------------------------------------------
// Free memory:
	free(pH_vector);

	for(indexi=0;indexi<NAP;indexi++){
		free(PetersonMatrix_reac[indexi]);
		free(KRateMatrix_reac[indexi]);
		free(rhoL[indexi]);
		free(rhoM[indexi]);
	}
	free(PetersonMatrix_reac);
	free(KRateMatrix_reac);
	free(rhoL);
	free(rhoM);

	for(indexi=0;indexi<5;indexi++){
		free(KHSConstant_reac[indexi]);
	}
	free(KHSConstant_reac);

	for(indexi=0;indexi<4;indexi++){
		free(TransMatrix_reac[indexi]);
	}
	free(TransMatrix_reac);

	free(Lin);
	free(Lpin);
	free(Min);
	free(Mpin);
	free(Lout);
	free(Mout);
	free(ReacL);
	free(ReacM);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int reacnsImplicit(double* input, double* pinput, double* output,
  void* user_reac_data, double V1, double V2)
{
//~ % 0 - 0,1 	- monomer sugar
//~ % 1 - 2,3 	- lactate
//~ % 2 - 4,5 	- hydrogen
//~ % 3 - 6,7 	- acetate
//~ % 4 - 8,9 	- proprionate
//~ % 5 - 10,11 - butyrate
//~ % 6 - 12,13 - methane
//~ % 7 - 14,15 - carbondioxide
//~ % 8 - 16,17 - water
//~ % 9 - 18,19 - polymer sugar (fiber/mucin)
//~ %
//~ % followed by:
//~ % Groups of Bacteria [odd numbers for mucus, even for lumen]

	int indexi,indexj,indexk;
	int rflag;
	int nsv = g_nsv;
	int pnsv = g_nsv/2;
	int MaxBacteria = g_mbacteria;

	int nSD = g_nsd;
	int nLD = g_nld;
	int nHDA = g_nhda;
	int nHDM = g_nhdm;
	int nOther = NOTHERS;
	double variance = g_bvariance;

	ReactionData reac_reac_data;
	reac_reac_data=(ReactionData)user_reac_data;

	double pH_value;
	double IpH;
	double IpH_power;
	double LocX;

	LocX = reac_reac_data->LocX;

	double* pH_vector = malloc(3 * sizeof(double));
	rflag = Clean1DArray(pH_vector,3);

		for(indexi=0;indexi<3;indexi++){
		pH_vector[indexi] = reac_reac_data->pH_v[indexi];
		}

	pH_value = SplineInterp(LocX,pH_vector[0],pH_vector[1],pH_vector[2]);

	IpH_power = -3.0*(pH_value - 6.69)/((6.69-5.8)*(6.69*5.8)); // from ADM1
	IpH = exp(IpH_power);

	double** PetersonMatrix_reac;
	double** KRateMatrix_reac;
	double** KHSConstant_reac;
	double** TransMatrix_reac;

	PetersonMatrix_reac = Make2DDoubleArray(NAP,pnsv);
	rflag = Clean2DArray(PetersonMatrix_reac,NAP,pnsv);
	KRateMatrix_reac = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(KRateMatrix_reac,NAP,MaxBacteria);
	KHSConstant_reac = Make2DDoubleArray(5,MaxBacteria);
	rflag = Clean2DArray(KHSConstant_reac,5,MaxBacteria);
	TransMatrix_reac = Make2DDoubleArray(4,pnsv);

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  PetersonMatrix_reac[indexi][indexj] = \
		  reac_reac_data->PMat[indexi][indexj];
	  }
    }

	for(indexi=0;indexi<NAP;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		 KRateMatrix_reac[indexi][indexj] = \
		 reac_reac_data->KRMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<5;indexi++){
	  for(indexj=0;indexj<MaxBacteria;indexj++){
		  KHSConstant_reac[indexi][indexj] = \
		  reac_reac_data->KHSMat[indexi][indexj];
	  }
    }

    for(indexi=0;indexi<4;indexi++){
	  for(indexj=0;indexj<pnsv;indexj++){
		  TransMatrix_reac[indexi][indexj] = \
		  reac_reac_data->TMat[indexi][indexj];
	  }
    }
//------------------------------------------------------------------------------	
// Set up problem---------------------------------------------------------------	

	double* Lin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lin,pnsv);

	double* Lpin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lpin,pnsv);

	double* Min = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Min,pnsv);

	double* Mpin = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Mpin,pnsv);

	double* Lout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Lout,pnsv);

	double* Mout = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(Mout,pnsv);

	double* ReacL = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacL,pnsv);

	double* ReacM = malloc(pnsv * sizeof (double));
	rflag = Clean1DArray(ReacM,pnsv);

	double temp=0.0;


indexj=0;
for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		Lin[indexj] = input[indexi];
		Lpin[indexj] = pinput[indexi];
	} else if(indexi%2==1) {
		Min[indexj] = input[indexi];
		Mpin[indexj] = pinput[indexi];
		indexj=indexj+1;
	}
}

	double** rhoL;
	rhoL = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoL,NAP,MaxBacteria);

	double** rhoM;
	rhoM = Make2DDoubleArray(NAP,MaxBacteria);
	rflag = Clean2DArray(rhoM,NAP,MaxBacteria);

//------------------------------------------------------------------------------
// Determine Biochemical reaction Rates-----------------------------------------

/*Hydrolysis (rho1);*/

//~ /*Contois 1 - OLD - Reminder of challenges of contois kinetics and comp.*/
//~ if (nSD > 0) {
  //~ for(indexi=0;indexi<nSD;indexi++){
  //~ rhoL[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
                                     //~ KHSConstant_reac[0][indexi], \
                                     //~ Lin[9], \
                                     //~ Lin[nOther+indexi]);
                                       //~
	//~ rhoM[0][indexi] = contois_kinetics(KRateMatrix_reac[0][indexi], \
	                                   //~ KHSConstant_reac[0][indexi], \
	                                   //~ Min[9], \
	                                   //~ Min[nOther+indexi]);
  //~ }
//~ }

if (g_hydro_scheme == 2){
  /*First-order disintegration*/

  rhoL[0][0] = fo_kinetics(10.62,Lin[9]);
  rhoM[0][0] = fo_kinetics(10.62,Min[9]);

} else if (g_hydro_scheme == 3){
  /*First-order Reaction*/
  if (nSD > 0) {
   for(indexi=0;indexi<nSD;indexi++){
    rhoL[0][indexi] = fo_reaction_kinetics(KRateMatrix_reac[0][indexi], \
                                           Lin[9], \
                                           Lin[nOther+indexi]);

    rhoM[0][indexi] = fo_reaction_kinetics(KRateMatrix_reac[0][indexi], \
                                           Min[9], \
                                           Min[nOther+indexi]);
   }
  }
} else {
 /*Contois 2*/
 double LcumGrowth = 0.0;
 double LcumHalfSat = 0.0;
 double McumGrowth = 0.0;
 double McumHalfSat = 0.0;

 if (nSD > 0) {
   for (indexi=0;indexi<nSD;indexi++){
     LcumGrowth = LcumGrowth + KRateMatrix_reac[0][indexi]*Lin[nOther+indexi];
     LcumHalfSat = LcumHalfSat + KHSConstant_reac[0][indexi]*Lin[nOther+indexi];
     McumGrowth = McumGrowth + KRateMatrix_reac[0][indexi]*Min[nOther+indexi];
     McumHalfSat = McumHalfSat + KHSConstant_reac[0][indexi]*Min[nOther+indexi];
   }
 }

 rhoL[0][0] = contois_kinetics_2(Lin[9],LcumGrowth,LcumHalfSat);
 rhoM[0][0] = contois_kinetics_2(Min[9],McumGrowth,McumHalfSat);
}

/*Glucose Utilization (rho2) and Death of Sugar Degrading Bacteria (rho6):*/
if (nSD > 0) {
  for(indexi=0;indexi<nSD;indexi++){
	rhoL[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Lin[0], \
	                                 Lin[nOther+indexi]);

	rhoM[1][indexi] = monod_kinetics(KRateMatrix_reac[1][indexi], \
	                                 KHSConstant_reac[1][indexi], \
	                                 Min[0], \
	                                 Min[nOther+indexi]);

	rhoL[5][indexi] = KRateMatrix_reac[5][indexi]*Lin[nOther+indexi];
	rhoM[5][indexi] = KRateMatrix_reac[5][indexi]*Min[nOther+indexi];
  }
}

/*Lactate Utilization (rho3) and Death of Lactate Degrading Bacteria (rho7):*/
if (nLD > 0) {
  for (indexi=0;indexi<nLD;indexi++){
	rhoL[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Lin[1], \
	                                 Lin[nOther+nSD+indexi]);

	rhoM[2][indexi] = monod_kinetics(KRateMatrix_reac[2][indexi], \
	                                 KHSConstant_reac[2][indexi], \
	                                 Min[1], \
	                                 Min[nOther+nSD+indexi]);

	rhoL[6][indexi] = KRateMatrix_reac[6][indexi]*Lin[nOther+nSD+indexi];
	rhoM[6][indexi] = KRateMatrix_reac[6][indexi]*Min[nOther+nSD+indexi];
  }
}

/*Homoacetogenesis (rho4) and Death of acetogenic bacteria (rho8):*/
if (nHDA > 0){
  for (indexi=0;indexi<nHDA;indexi++){
	rhoL[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
	                                 KHSConstant_reac[3][indexi], \
	                                 Lin[2], \
	                                 Lin[nOther+nSD+nLD+indexi]);

	rhoM[3][indexi] = monod_kinetics(KRateMatrix_reac[3][indexi], \
	                                 KHSConstant_reac[3][indexi], \
	                                 Min[2], \
	                                 Min[nOther+nSD+nLD+indexi]);

       rhoL[7][indexi] = KRateMatrix_reac[7][indexi]*Lin[nOther+nSD+nLD+indexi];
       rhoM[7][indexi] = KRateMatrix_reac[7][indexi]*Min[nOther+nSD+nLD+indexi];
  }
}

/*Homomethanogenesis (rho5) and Death of methanogenic bacteria (rho8):*/
if (nHDM > 0){
  for (indexi=0;indexi<nHDM;indexi++){
	rhoL[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                     KHSConstant_reac[4][indexi], \
	                                     Lin[2], \
	                                     Lin[nOther+nSD+nLD+nHDA+indexi]);

	rhoM[4][indexi] = IpH*monod_kinetics(KRateMatrix_reac[4][indexi], \
	                                     KHSConstant_reac[4][indexi], \
	                                     Min[2], \
	                                     Min[nOther+nSD+nLD+nHDA+indexi]);

  rhoL[8][indexi] = KRateMatrix_reac[8][indexi]*Lin[nOther+nSD+nLD+nHDA+indexi];
  rhoM[8][indexi] = KRateMatrix_reac[8][indexi]*Min[nOther+nSD+nLD+nHDA+indexi];
  }
}
//------------------------------------------------------------------------------
// Build Differential Equations-------------------------------------------------
double V = g_vcolon_l;
double Vd = 0.9*V;
double Vm = 0.1*V;

double MuMax = 50.0; // [g/L]
double MuPro = 500.0; // [g/Ld]
double MuGen = 0.0;

double InterpTransport= 0.0;
double fr = g_fr_lpd;
double vratio = V1/V2;
double drate = fr/V2;

MuGen = (1.0 - Min[9]/MuMax)*MuPro;

  /*Non-biomass terms*/
  for(indexi=0;indexi<nOther;indexi++){
	ReacL[indexi]=0;
	ReacM[indexi]=0;

	for(indexj=0;indexj<9;indexj++){
	  for(indexk=0;indexk<MaxBacteria;indexk++){
	  ReacL[indexi] = ReacL[indexi] +  \
	               PetersonMatrix_reac[indexj][indexi]*rhoL[indexj][indexk];
	  ReacM[indexi] = ReacM[indexi] + \
	               PetersonMatrix_reac[indexj][indexi]*rhoM[indexj][indexk];
	  }
	}

	if(indexi==0){
	InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] - \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vd;

	Mout[indexi] = ReacM[indexi] + \
	               InterpTransport*(Lin[indexi]-Min[indexi])/Vm;

	}else if(indexi==nOther-1){
	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] + \
	               (Vm/Vd)*TransMatrix_reac[3][indexi]*Min[indexi] ;

	Mout[indexi] = MuGen + ReacM[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];

	}else{
	InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][indexi], \
	TransMatrix_reac[1][indexi],TransMatrix_reac[2][indexi]);

	Lout[indexi] = drate*(vratio*Lpin[indexi] - Lin[indexi])+ \
	               ReacL[indexi] - \
	               InterpTransport*Lin[indexi];

	Mout[indexi] = ReacM[indexi] + \
	               (Vd/Vm)*InterpTransport*Lin[indexi] - \
	               TransMatrix_reac[3][indexi]*Min[indexi];
	}

  }

  /*Biomass Terms*/
  for(indexi=0;indexi<nSD;indexi++){
	ReacL[nOther+indexi] = 0;
	ReacM[nOther+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+indexi] = ReacL[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+indexi] = ReacM[nOther+indexi] + \
		PetersonMatrix_reac[indexj][nOther+indexi]*rhoM[indexj][indexi];
	}

  InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+indexi], \
  TransMatrix_reac[1][nOther+indexi],TransMatrix_reac[2][nOther+indexi]);

  Lout[nOther+indexi] = drate*(vratio*Lpin[nOther+indexi]-Lin[nOther+indexi])+\
  ReacL[nOther+indexi] - \
  InterpTransport*Lin[nOther+indexi] + \
  (Vm/Vd)*TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];

	Mout[nOther+indexi] = ReacM[nOther+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+indexi] - \
	TransMatrix_reac[3][nOther+indexi]*Min[nOther+indexi];
  }

  for(indexi=0;indexi<nLD;indexi++){
	ReacL[nOther+nSD+indexi] = 0;
	ReacM[nOther+nSD+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+indexi]*rhoM[indexj][indexi];
	}

InterpTransport = SplineInterp(LocX,TransMatrix_reac[0][nOther+nSD+indexi],\
		          TransMatrix_reac[1][nOther+nSD+indexi],\
                  TransMatrix_reac[2][nOther+nSD+indexi]);

	Lout[nOther+nSD+indexi] = ReacL[nOther+nSD+indexi] - \
	InterpTransport*Lin[nOther+nSD+indexi] + \
       (Vm/Vd)*TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi] + \
	drate*(vratio*Lpin[nOther+nSD+indexi] - Lin[nOther+nSD+indexi]);

	Mout[nOther+nSD+indexi] = ReacM[nOther+nSD+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+indexi] - \
	TransMatrix_reac[3][nOther+nSD+indexi]*Min[nOther+nSD+indexi];
  }

  for(indexi=0;indexi<nHDA;indexi++){
	ReacL[nOther+nSD+nLD+indexi] = 0;
	ReacM[nOther+nSD+nLD+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
		ReacL[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoL[indexj][indexi];

		ReacM[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
	PetersonMatrix_reac[indexj][nOther+nSD+nLD+indexi]*rhoM[indexj][indexi];
	}
	InterpTransport = SplineInterp(LocX, \
	                          TransMatrix_reac[0][nOther+nSD+nLD+indexi], \
	                          TransMatrix_reac[1][nOther+nSD+nLD+indexi], \
	                          TransMatrix_reac[2][nOther+nSD+nLD+indexi]);

Lout[nOther+nSD+nLD+indexi] = ReacL[nOther+nSD+nLD+indexi] - \
                              InterpTransport*Lin[nOther+nSD+nLD+indexi] + \
(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi] +\
drate*(vratio*Lpin[nOther+nSD+nLD+indexi] - Lin[nOther+nSD+nLD+indexi]);

Mout[nOther+nSD+nLD+indexi] = ReacM[nOther+nSD+nLD+indexi] + \
(Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+indexi] - \
TransMatrix_reac[3][nOther+nSD+nLD+indexi]*Min[nOther+nSD+nLD+indexi];
}

  for(indexi=0;indexi<nHDM;indexi++){
	ReacL[nOther+nSD+nLD+nHDA+indexi] = 0;
	ReacM[nOther+nSD+nLD+nHDA+indexi] = 0;
	for(indexj=0;indexj<9;indexj++){
	ReacL[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] +\
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoL[indexj][indexi];

       ReacM[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
   PetersonMatrix_reac[indexj][nOther+nSD+nLD+nHDA+indexi]*rhoM[indexj][indexi];
	}

	InterpTransport = SplineInterp(LocX,\
	TransMatrix_reac[0][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[1][nOther+nSD+nLD+nHDA+indexi],\
	TransMatrix_reac[2][nOther+nSD+nLD+nHDA+indexi]);

	Lout[nOther+nSD+nLD+nHDA+indexi] = ReacL[nOther+nSD+nLD+nHDA+indexi] - \
	InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] + \
	(Vm/Vd)*TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi] + \
	drate*(vratio*Lpin[nOther+nSD+nLD+nHDA+indexi] - \
	Lin[nOther+nSD+nLD+nHDA+indexi]);

	Mout[nOther+nSD+nLD+nHDA+indexi] = ReacM[nOther+nSD+nLD+nHDA+indexi] + \
	(Vd/Vm)*InterpTransport*Lin[nOther+nSD+nLD+nHDA+indexi] - \
	TransMatrix_reac[3][nOther+nSD+nLD+nHDA+indexi]\
	*Min[nOther+nSD+nLD+nHDA+indexi];
  }

  /*Return Completed Vector*/
  indexj=0;
  for(indexi=0;indexi<nsv;indexi++){
	if(indexi%2==0){
		output[indexi] = Lout[indexj];
	} else if(indexi%2==1){
		output[indexi] = Mout[indexj];
		indexj=indexj+1;
	}
  }

//------------------------------------------------------------------------------
// Free memory:
	free(pH_vector);

	for(indexi=0;indexi<NAP;indexi++){
		free(PetersonMatrix_reac[indexi]);
		free(KRateMatrix_reac[indexi]);
		free(rhoL[indexi]);
		free(rhoM[indexi]);
	}
	free(PetersonMatrix_reac);
	free(KRateMatrix_reac);
	free(rhoL);
	free(rhoM);

	for(indexi=0;indexi<5;indexi++){
		free(KHSConstant_reac[indexi]);
	}
	free(KHSConstant_reac);

	for(indexi=0;indexi<4;indexi++){
		free(TransMatrix_reac[indexi]);
	}
	free(TransMatrix_reac);

	free(Lin);
	free(Lpin);
	free(Min);
	free(Mpin);
	free(Lout);
	free(Mout);
	free(ReacL);
	free(ReacM);

return(0);


}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
static int MOL_func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	int indexi;
	int aflag;
	int nsv = g_nsv;
	double* PRin = malloc(nsv * sizeof (double));
	aflag = Clean1DArray(PRin,nsv);
	double* PRout = malloc(nsv * sizeof (double));
	aflag = Clean1DArray(PRout,nsv);
	double* Conv = malloc(nsv * sizeof (double));
	aflag = Clean1DArray(Conv,nsv);

	ReactionData f_reac_data;
	f_reac_data = (ReactionData)user_data;

	for(indexi=0;indexi<nsv;indexi++){
		PRin[indexi] = Ith(y,indexi+1);
		Conv[indexi] = f_reac_data->interim[indexi];
	}

	  //printf("\n%f\n",f_reac_data->delt);
	  //PauseMessage(1571);

		aflag = reac(PRin,PRout,f_reac_data);


	for(indexi=0;indexi<nsv;indexi++){
		Ith(ydot,indexi+1) = Conv[indexi] + PRout[indexi];
	}

	free(PRin);
	free(PRout);
	free(Conv);

  return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
int MOL_ODEsolver(double* IG, double* Sol, void* user_reac_data)
{
	int nsv = g_nsv;
	int indexi=0;
	int flag;

	realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;

  ReactionData cvODE_reac_data;
  cvODE_reac_data = (ReactionData)user_reac_data;
  //~ PauseMessage(1489);

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(nsv);
  abstol = N_VNew_Serial(nsv);

  /* Initial Condition */
  for(indexi=0;indexi<nsv;indexi++){
   Ith(y,indexi+1) = IG[indexi];
}

	/* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
	for(indexi=0;indexi<nsv;indexi++){
  Ith(abstol,indexi+1) = ATOLs; //fmax(ATOLs,IG[indexi]*ATOLs);
  }

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, MOL_func, 0, y);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);

  /* Set pointer to reaction information */
  flag = CVodeSetUserData(cvode_mem,cvODE_reac_data);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, nsv);

  /* Solve ODE */
  t=0.0;
	tout = cvODE_reac_data->dt;
	flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

	/* Return Solution */
  for(indexi=0;indexi<nsv;indexi++){
	  Sol[indexi] = Ith(y,indexi+1);
  }

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double upwind_disc(double Cm, double Bm, double Dm, double Em,
  double gamma, double vl, int order)
 {
 /*-----------------------------------------------------------------------------
 *  STENCIL: cell centered finite volume (original solutions are stored at
 *  cell corners
 *
 *  -------------------------------
 *  |     |     |     |     |     |
 *  |     |	|  A  |     |     |
 *  |     |     |     |	    |     |
 *  -------------------------------
 *  |     |     |     |     |     |
 *  |  D  |  B  |  C  |  E  |     |
 *  |     |     |     |	    |     |
 *  -------------------------------
 *
 * 	FO-Upwind:
 *
 *  conv = vl*gamma*(B-C)
 *
 *  SO-Upwind:
 *  conv = vl*gamma*0.5*(4*B - 3C - D)
 *
 *  TO-Upwind:
 *  conv = vl*gamma*0.1666667*(6*B-3*C-D-2*E)
 *
 *----------------------------------------------------------------------------*/

	double Conv;

	if (order == 1){
	  Conv = vl*gamma*(Bm-Cm);
	}else if (order == 2) {
	  Conv = vl*gamma*0.5*(4*Bm - Dm - 3*Cm);
	}else{
		Conv = vl*gamma*0.166666667*(6*Bm - Dm - 3*Cm - 2*Em);
  }

	Conv = fmax(0,Conv);

return (Conv);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int MOL_numScheme(double** PSolution, double* Input, int ngrids,
     void* user_reac_data)
{
/*--------------------------------------------------------------------------
 * MOL numerical scheme (testing):
 * Breaks down to an edge centred finite difference difference scheme with
 * explicit upwinding of the convection term (see upwind_disc.c):
 *
 *  ------------------A------------
 *  |     |     |     |     |     |
 *  |     |	|     |     |     |
 *  |     |     |     |	    |     |
 *  ------D-----B-----C-----E------
 *
 * and is evaluated using 4 solutions from the previous time step (DBCE) A
 * backwards difference formula (SUNDIALS - cvODE) is used to evaluate the
 * resulting set of ODEs.
 */

  int indexi,indexj,indexk;
	int aflag;
	int nsv = g_nsv;
	double fr_mpd = g_fr_mpd;

	double**	iSolution = Make2DDoubleArray(ngrids,nsv);
	aflag = Clean2DArray(iSolution,ngrids,nsv);

	// Reaction Terms
	ReactionData MOL_nS_reac_data;
	MOL_nS_reac_data = (ReactionData)user_reac_data;

	double* Am = malloc(nsv * sizeof(double)); aflag = Clean1DArray(Am,nsv);
	double* Bm = malloc(nsv * sizeof(double)); aflag = Clean1DArray(Bm,nsv);
	double* Cm = malloc(nsv * sizeof(double)); aflag = Clean1DArray(Cm,nsv);
	double* Dm = malloc(nsv * sizeof(double)); aflag = Clean1DArray(Dm,nsv);
	double* Em = malloc(nsv * sizeof(double)); aflag = Clean1DArray(Em,nsv);

     double* Conv = malloc(nsv * sizeof(double));aflag = Clean1DArray(Conv,nsv);
	double* IS = malloc(nsv * sizeof(double));aflag = Clean1DArray(IS,nsv);

	double vl;
	double delt = MOL_nS_reac_data->dt;
	double delx = MOL_nS_reac_data->dx;

	double gamma = 1/delx;
	double lcolon = g_lcolon_m;


  for (indexi=0;indexi<ngrids;indexi++){
	MOL_nS_reac_data->LocX = indexi*delx/lcolon;

		if(indexi==0){
		  for(indexj=0;indexj<nsv;indexj++){
		    Dm[indexj] = 2*Input[indexj] - PSolution[indexi][indexj];
		    Bm[indexj] = Input[indexj];
		    Cm[indexj] = PSolution[indexi][indexj];
		    Em[indexj] = PSolution[indexi+1][indexj];
		  }
		} else if (indexi==1){
			for(indexj=0;indexj<nsv;indexj++){
				Dm[indexj] = Input[indexj];
				Bm[indexj] = PSolution[indexi-1][indexj];
				Cm[indexj] = PSolution[indexi][indexj];
				Em[indexj] = PSolution[indexi+1][indexj];
			}
		} else if (indexi==ngrids-1){
			for(indexj=0;indexj<nsv;indexj++){
				Dm[indexj] = PSolution[indexi-2][indexj];
				Bm[indexj] = PSolution[indexi-1][indexj];
				Cm[indexj] = PSolution[indexi][indexj];
				Em[indexj] = 2*Cm[indexj]-Bm[indexj];
			}
		} else {
			for(indexj=0;indexj<nsv;indexj++){
				Dm[indexj] = PSolution[indexi-2][indexj];
				Bm[indexj] = PSolution[indexi-1][indexj];
				Cm[indexj] = PSolution[indexi][indexj];
				Em[indexj] = PSolution[indexi+1][indexj];
			}
		}

  //printf("\n");
  for(indexj=0;indexj<nsv;indexj++){
    if(indexj%2==0){vl=fr_mpd;} else {vl=0.0;}
      Conv[indexj] = upwind_disc(Cm[indexj],Bm[indexj],Dm[indexj],Em[indexj],\
	             gamma,vl,g_MOL_sch_order);

      Conv[indexj] = fmax(1e-26,Conv[indexj]);
    if(indexj%2==0){IS[indexj] = Cm[indexj];} else {IS[indexj] = Cm[indexj];}

	MOL_nS_reac_data->interim[indexj] = Conv[indexj];
  }

  aflag = MOL_ODEsolver(IS,Am,MOL_nS_reac_data);

  for(indexj=0;indexj<nsv;indexj++){
    iSolution[indexi][indexj] = Am[indexj];
  }


  } //end for indexi

	for(indexi=0;indexi<ngrids;indexi++){
		for(indexj=0;indexj<nsv;indexj++){
			PSolution[indexi][indexj] = iSolution[indexi][indexj];
		}
	}

//------------------------------------------------------------------------------
// Free function memory:

	for(indexi=0;indexi<ngrids;indexi++){
          free(iSolution[indexi]);
	}
	free(iSolution);

	free(Am);
	free(Bm);
	free(Cm);
	free(Dm);
	free(Em);

	free(Conv);
	free(IS);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
static int TSRS_func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int indexi,indexj,indexk;
  int aflag;
  int nsv = 3*g_nsv;
  double* PRin = malloc(nsv * sizeof (double));
  aflag = Clean1DArray(PRin,nsv);
  double* PRout = malloc(nsv * sizeof (double));
  aflag = Clean1DArray(PRout,nsv);

  double* C = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double* P = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double* O = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double Vc;
  double Vp;

  ReactionData f_reac_data;
  f_reac_data = (ReactionData)user_data;

  for(indexi=0;indexi<nsv;indexi++){
  	PRin[indexi] = Ith(y,indexi+1);
  }

  // convert from 3n to 3 sets of n variables (for reaction terms)
  for(indexi=0;indexi<3;indexi++){
    if(indexi==0){

  	  for(indexj=0;indexj<g_nsv;indexj++){
  			P[indexj] = f_reac_data->interim[indexj];
  			C[indexj] = PRin[indexi*g_nsv + indexj];
  	  }

  	  Vp = g_vproximal_l;
  	  Vc = g_vproximal_l;
  	  f_reac_data->LocX = 0;

  		  aflag = reac3sImplicit(C,P,O,f_reac_data,Vp,Vc);

  		  for(indexj=0;indexj<g_nsv;indexj++){
  				PRout[indexi*g_nsv + indexj] = O[indexj];
  			}
	//---------------
	  } else if (indexi==1) {

			for(indexj=0;indexj<g_nsv;indexj++){
				P[indexj] = PRin[(indexi-1)*g_nsv + indexj];
				C[indexj] = PRin[indexi*g_nsv + indexj];
		  }

		  Vp = g_vproximal_l;
		  Vc = g_vtransverse_l;
		  f_reac_data->LocX = 0.3;

		  aflag = reac3sImplicit(C,P,O,f_reac_data,Vp,Vc);

		  for(indexj=0;indexj<g_nsv;indexj++){
				PRout[indexi*g_nsv + indexj] = O[indexj];
			}
		//-----------------
		} else if (indexi==2) {
			for(indexj=0;indexj<g_nsv;indexj++){
				P[indexj] = PRin[(indexi-1)*g_nsv + indexj];
				C[indexj] = PRin[indexi*g_nsv + indexj];
		  }

		  Vp = g_vtransverse_l;
		  Vc = g_vdistal_l;
		  f_reac_data->LocX = 0.9;

		  aflag = reac3sImplicit(C,P,O,f_reac_data,Vp,Vc);

		  for(indexj=0;indexj<g_nsv;indexj++){
				PRout[indexi*g_nsv + indexj] = O[indexj];
			}


		}
	//-----------------

	}

	for(indexi=0;indexi<nsv;indexi++){
		Ith(ydot,indexi+1) = PRout[indexi];
	}

	free(PRin);
	free(PRout);
	free(C);
	free(P);
	free(O);

  return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
int TSRS_ODEsolver(double* IG, double* Sol, void* user_reac_data)
{
	int nsv = 3*g_nsv;
	int indexi=0;
	int flag;

	realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;

  ReactionData cvODE_reac_data;
  cvODE_reac_data = (ReactionData)user_reac_data;

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(nsv);
  abstol = N_VNew_Serial(nsv);

  /* Initial Condition */
  for(indexi=0;indexi<nsv;indexi++){
   Ith(y,indexi+1) = IG[indexi];
}

	/* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
	for(indexi=0;indexi<nsv;indexi++){
  Ith(abstol,indexi+1) = ATOLs; //fmax(ATOLs,IG[indexi]*ATOLs);
  }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, TSRS_func, 0, y);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);

  /* Set pointer to reaction information */
  flag = CVodeSetUserData(cvode_mem,cvODE_reac_data);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, nsv);

  /* Solve ODE */
  t=0.0;
	tout = cvODE_reac_data->dt;
	flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

  /* Return Solution */
  for(indexi=0;indexi<nsv;indexi++){
	  Sol[indexi] = Ith(y,indexi+1);
  }

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int ThreeStageReactorScheme(double **PSolution, double *Input,
  void* user_reac_data)
{
	/* Reformulate the system from 3 reactors with n variables to a complete 
	 * system of 3n variables (ie, 3 x 28 should be a vector of length 84
	 */

	 int indexi,indexj,indexk;
	 int aflag;
	 int nEQ = 3*g_nsv;

	 double* F = malloc(nEQ * sizeof (double));
	 aflag = Clean1DArray(F,nEQ);
	 double* G = malloc(nEQ * sizeof (double));
	 aflag = Clean1DArray(G,nEQ);

	 ReactionData TSRS_reac_data;
	 TSRS_reac_data = (ReactionData)user_reac_data;

	 indexk=0;
	 for(indexi=0;indexi<3;indexi++){
		 for(indexj=0;indexj<g_nsv;indexj++){
			 F[indexk] = PSolution[indexi][indexj];
			 indexk++;
		 }
	 }

	 for(indexj=0;indexj<g_nsv;indexj++){
		 TSRS_reac_data->interim[indexj] = Input[indexj];
	 }

	 aflag = TSRS_ODEsolver(F,G,TSRS_reac_data);

	 indexk=0;
	 for(indexi=0;indexi<3;indexi++){
		 for(indexj=0;indexj<g_nsv;indexj++){
			 PSolution[indexi][indexj] = G[indexk];
			 indexk++;
		 }
	 }

	 free(F);
	 free(G);
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
static int GS_func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int indexi,indexj,indexk;
  int aflag;
  int nsv = g_nreactors*g_nsv;
  double* PRin = malloc(nsv * sizeof (double));
  aflag = Clean1DArray(PRin,nsv);
  double* PRout = malloc(nsv * sizeof (double));
  aflag = Clean1DArray(PRout,nsv);

  double* C = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double* P = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double* O = malloc(g_nsv * sizeof (double)); aflag = Clean1DArray(C,g_nsv);
  double Vc;
  double Vp;

	ReactionData f_reac_data;
	f_reac_data = (ReactionData)user_data;

	for(indexi=0;indexi<nsv;indexi++){
		PRin[indexi] = Ith(y,indexi+1);
	}

	// convert from 3n to 3 sets of n variables (for reaction terms)
	for(indexi=0;indexi<g_nreactors;indexi++){
		f_reac_data->LocX = (double)indexi/ (double)g_nreactors;
	  if(indexi==0){

		  for(indexj=0;indexj<g_nsv;indexj++){
				P[indexj] = f_reac_data->interim[indexj];
				C[indexj] = PRin[indexi*g_nsv + indexj];
		  }

		  Vp = g_vcolon_l/g_nreactors;
		  Vc = g_vcolon_l/g_nreactors;

		  aflag = reacnsImplicit(C,P,O,f_reac_data,Vp,Vc);

		  for(indexj=0;indexj<g_nsv;indexj++){
				PRout[indexi*g_nsv + indexj] = O[indexj];
			}
		//-----------------
		} else {
			for(indexj=0;indexj<g_nsv;indexj++){
				P[indexj] = PRin[(indexi-1)*g_nsv + indexj];
				C[indexj] = PRin[indexi*g_nsv + indexj];
		  }

		  Vp = g_vcolon_l/g_nreactors;
		  Vc = g_vcolon_l/g_nreactors;

		  aflag = reac3sImplicit(C,P,O,f_reac_data,Vp,Vc);

		  for(indexj=0;indexj<g_nsv;indexj++){
				PRout[indexi*g_nsv + indexj] = O[indexj];
			}


		}
	//-----------------

	}

	for(indexi=0;indexi<nsv;indexi++){
		Ith(ydot,indexi+1) = PRout[indexi];
	}

	free(PRin);
	free(PRout);
	free(C);
	free(P);
	free(O);

  return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
int GS_ODEsolver(double* IG, double* Sol, int nReactors,
  void* user_reac_data)
{
	int nsv = nReactors*g_nsv;
	int indexi=0;
	int flag;

	realtype reltol, t, tout;
  N_Vector y, abstol;
  void *cvode_mem;

  ReactionData cvODE_reac_data;
  cvODE_reac_data = (ReactionData)user_reac_data;

  y = abstol = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(nsv);
  abstol = N_VNew_Serial(nsv);

  /* Initial Condition */
  for(indexi=0;indexi<nsv;indexi++){
   Ith(y,indexi+1) = IG[indexi];
}

  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
	for(indexi=0;indexi<nsv;indexi++){
  Ith(abstol,indexi+1) = ATOLs; //fmax(ATOLs,IG[indexi]*ATOLs);
  }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, GS_func, 0, y);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);

  /* Set pointer to reaction information */
  flag = CVodeSetUserData(cvode_mem,cvODE_reac_data);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, nsv);

  /* Solve ODE */
  t=0.0;
	tout = cvODE_reac_data->dt;
	flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

	/* Return Solution */
  for(indexi=0;indexi<nsv;indexi++){
	  Sol[indexi] = Ith(y,indexi+1);
  }

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int GradostatScheme(double **PSolution, double *Input, int nReactors,
  void* user_reac_data)
{
	 int indexi,indexj,indexk;
	 int aflag;
	 int nEQ = nReactors*g_nsv;

	 double* F = malloc(nEQ * sizeof (double));
	 aflag = Clean1DArray(F,nEQ);
	 double* G = malloc(nEQ * sizeof (double));
	 aflag = Clean1DArray(G,nEQ);

	 ReactionData GS_reac_data;
	 GS_reac_data = (ReactionData)user_reac_data;

	 indexk=0;
	 for(indexi=0;indexi<nReactors;indexi++){
		 for(indexj=0;indexj<g_nsv;indexj++){
			 F[indexk] = PSolution[indexi][indexj];
			 indexk++;
		 }
	 }

	 for(indexj=0;indexj<g_nsv;indexj++){
		 GS_reac_data->interim[indexj] = Input[indexj];
	 }

	 aflag = GS_ODEsolver(F,G,nReactors,GS_reac_data);

	 indexk=0;
	 for(indexi=0;indexi<nReactors;indexi++){
		 for(indexj=0;indexj<g_nsv;indexj++){
			 PSolution[indexi][indexj] = G[indexk];
			 indexk++;
		 }
	 }

	 free(F);
	 free(G);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double LengthInterpolator(double a, double c)
{
        double b;

        b = (a+c)/2.0;

return b;
} //end function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double Flux(double v, double c)
{
	double conv;
	double flux;

	flux = v;
	conv = flux*c;

return conv;
} //end function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double MinMod(double a, double b)
{
	double fluxLimiter;
	int sa, sb;

	if (a >= 0) {sa = 1;} else if (a < 0) {sa = -1; a = a*-1;}
	if (b >= 0) {sb = 1;} else if (b < 0) {sb = -1; b = b*-1;}

fluxLimiter = 0.5*(sa + sb)*fmin(a,b);

return fluxLimiter;
} //end function
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double FluxLimiter(double a1, double a2, double a3, double a4, double a5)
{
double Dp,Dn;
double uprime;

// I,C,B,D,J or K,I,C,B,D from discrete stencil

Dp = MinMod(a1+2*a2+a3, a2 + 2*a3 + a4);
Dn = MinMod(a2+2*a3+a4, a3 + 2*a4 + a5);

uprime = MinMod(a2 - a3 - 0.5*Dp, a3 - a4 + 0.5*Dn);

return(uprime);

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
static int func(N_Vector u, N_Vector f, void *user_data)
{
	int indexi,indexj,indexk;
	int aflag;
	int nsv = g_nsv;

	realtype *udata, *fdata;

	ReactionData data;
	data = (ReactionData)user_data;

	double tfac;
	double delt;
	double* interim = malloc(nsv * sizeof (double));
	aflag = Clean1DArray(interim,nsv);

	tfac = data->tfac;
	delt = data->dt/2.0;

	for(indexj=0;indexj<nsv;indexj++){
	  interim[indexj] = data->interim[indexj];
    }

	udata = NV_DATA_S(u);
	fdata = NV_DATA_S(f);

	double* rdata = malloc(nsv * sizeof (double));

	aflag = reac(udata,rdata,data);

    for(indexi=0;indexi<nsv;indexi++){
    fdata[indexi] = interim[indexi] + tfac*delt*rdata[indexi] - udata[indexi];
    }

	free(interim);
	free(rdata);
	//~ free(udata);
	//~ free(fdata);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// SUNDIALS
int SD_nl_solver(void* user_reaction_data, N_Vector u, double* IG)
{
	int indexi,indexj,indexk;
	int nsv = g_nsv;

	ReactionData data;
	data = (ReactionData)user_reaction_data;

	realtype fnormtol, scsteptol;
	N_Vector u1,s,c;
	int glstr, mset, flag;
	void *kmem;

	u1 = s = c = NULL;
	kmem = NULL;
	glstr = KIN_NONE;

	u1= N_VNew_Serial(nsv);
	s = N_VNew_Serial(nsv);
	c = N_VNew_Serial(nsv);

	// set u1 to be initial u
	for (indexi=0;indexi<nsv;indexi++){
		Ith(u1,indexi+1) = fmax(0,IG[indexi]);
	}
	//PauseMessage();

	N_VConst_Serial(1.0,s);
	N_VConst_Serial(1.0,c);

	//~ for(indexi=0;indexi<nsv;indexi++){
		//~ printf("%lf\n",Ith(u1,indexi+1));
	//~ }
	//~ PauseMessage(1544);

	fnormtol = FTOL; scsteptol=STOL;

	kmem = KINCreate();

	flag = KINSetUserData(kmem,data);
	flag = KINSetConstraints(kmem,c);
	flag = KINSetFuncNormTol(kmem,fnormtol);
	flag = KINSetScaledStepTol(kmem,scsteptol);
	flag = KINInit(kmem,func,u);
	flag = KINDense(kmem,nsv);

	N_VScale_Serial(1.0,u1,u);
	glstr = KIN_NONE;
	mset = 1;
	flag = KINSetMaxSetupCalls(kmem,mset);
	flag = KINSol(kmem,u,glstr,s,s);

	N_VDestroy_Serial(u1);
	N_VDestroy(s);
	N_VDestroy(c);
	KINFree(&kmem);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int Non_lin_solver(double* Sol,double* IG,int psize,void* user_reac_data)
{
	int i,nls_flag;
	ReactionData nls_reac_data;
	nls_reac_data = (ReactionData)user_reac_data;

	N_Vector u;
	u=NULL;
	u = N_VNew_Serial(psize);

	for(i=0;i<psize;i++){
		nls_reac_data->interim[i] = IG[i];
	}

	SD_nl_solver(nls_reac_data,u,IG);

	for(i=0;i<psize;i++){
		Sol[i] = Ith(u,i+1);
	}

	//~ printf("\n%d\n",gLocation);
	//~ gLocation++;
	//~ printf("\nSolution:\n");
	//~ for(i=0;i<psize;i++){
		//~ printf("%f %f\n",IG[i],Sol[i]);
	//~ }
	//~
	//~ PauseMessage(1602);

	N_VDestroy_Serial(u);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int Liotta_num_method(double* Am,double* Bm,double* Cm,double* Dm,
     double* Em, double* Fm,double* Gm,double* Hm,double* Im,
     double* Jm, double* Km,void* user_reac_data)
{
	int indexi,indexj,indexk;
	int aflag;
	int nsv = g_nsv;
	double fr_mpd = g_fr_mpd;

	// Reaction Terms
	ReactionData ho_nM_reac_data;
	ho_nM_reac_data = (ReactionData)user_reac_data;

	double delt = ho_nM_reac_data->dt/2;
	double delx = ho_nM_reac_data->dx;
        double gamma;
        double vl;

	gamma = delt/delx;

	double trap_factor;

	trap_factor = 3./8;

	aflag = reac(Gm,Gm,ho_nM_reac_data);
	aflag = reac(Hm,Hm,ho_nM_reac_data);



for(indexj=0;indexj<nsv;indexj++){
if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
Am[indexj]= ho_nM_reac_data->interim[indexj] = \
  0.5*(Bm[indexj]+Cm[indexj]) +
  0.125*( FluxLimiter(Im[indexj],Cm[indexj],Bm[indexj],Dm[indexj],Jm[indexj]) -
        FluxLimiter(Km[indexj],Im[indexj],Cm[indexj],Bm[indexj],Dm[indexj]) ) -
	gamma*( Flux(vl,Fm[indexj]) - Flux(vl,Em[indexj]) ) +
	delt*(trap_factor*Gm[indexj] + trap_factor*Hm[indexj]);
	}

   double* IG = malloc(nsv * sizeof (double)); aflag = Clean1DArray(IG,nsv);

	for(indexj=0;indexj<nsv;indexj++){
	 //~ if(indexj%2==0){
		IG[indexj] = Am[indexj];
	 //~ } else {
		//~ IG[indexj] = 0.5*(Bm[indexj] + Cm[indexj]);
	 //~ }
	}

ho_nM_reac_data->tfac = 0.25;
aflag = Non_lin_solver(Am,IG,nsv,ho_nM_reac_data);






free(IG);

return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int Liotta_numScheme(double** PSolution, double* Input, int ngrids,
    void* user_reac_data)
{
/*--------------------------------------------------------------------------
 * Numerical Scheme based on: Liota et al. Central Schemes for Balance Laws
 * of Relaxation Type. SIAM Journal on Numerical Analysis. Vol 38. No 4.
 * (2001)
 *
 * Breaks down to an edge centred implicit-explicit finite difference scheme
 * with the following stencil:
 *
 *  ------ ----- --A-- ----- -----
 *  |     |     |     |     |     |
 *  |     |	    E     F     |     |
 *  |     |     G     H	    |     |
 *  J-----D-----B-----C-----I-----K
 *
 * where A is the solution at the current time steped, determined using 5 
 * solutions from the previous time step (J,D,B,C,I,K) and 4 implicitly 
 * calculated intermediate time steps (1/3 time step = G,H; 
 * 1/2 time step = E,F). 
 */ 
		
	int indexi,indexj,indexk;
	int aflag;
	int nsv = g_nsv;
	double fr_mpd = g_fr_mpd;
	
	double** iSolution = Make2DDoubleArray(ngrids-1,nsv);
	aflag = Clean2DArray(iSolution,ngrids-1,nsv);

  double* Am = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Am,nsv);
  double* Bm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Bm,nsv);
  double* Cm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Cm,nsv);
  double* Dm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Dm,nsv);
  double* Em = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Em,nsv);
  double* Fm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Fm,nsv);
  double* Gm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Gm,nsv);
  double* Hm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Hm,nsv);
  double* Im = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Im,nsv);
  double* Jm = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Jm,nsv);
  double* Km = malloc(nsv * sizeof (double)); aflag = Clean1DArray(Km,nsv);
  double* iEm = malloc(nsv * sizeof (double));aflag = Clean1DArray(iEm,nsv);
  double* iFm = malloc(nsv * sizeof (double));aflag = Clean1DArray(iFm,nsv);
  double* iGm = malloc(nsv * sizeof (double));aflag = Clean1DArray(iGm,nsv);
  double* iHm = malloc(nsv * sizeof (double));aflag = Clean1DArray(iHm,nsv);
  double vl;
  double gamma;

	// Reaction Terms
	ReactionData ho_nS_reac_data;
	ho_nS_reac_data = (ReactionData)user_reac_data;
	
	double delt = ho_nS_reac_data->dt/2;
	double delx = ho_nS_reac_data->dx;
        gamma = delt/delx;
    
  double lcolon = g_lcolon_m;

int oc = 0;

//------------------------------------------------------------------------------
//Part 1:
  
for(indexi=0;indexi<ngrids-1;indexi++){
  ho_nS_reac_data->LocX = indexi*delx/lcolon;
	
  if (indexi==0){
		
    for (indexj=0;indexj<nsv;indexj++){
      if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	Jm[indexj] = Input[indexj];
	Dm[indexj] = Input[indexj];
	Bm[indexj] = PSolution[indexi][indexj];
	Cm[indexj] = PSolution[indexi+1][indexj];
	Im[indexj] = PSolution[indexi+2][indexj];
	Km[indexj] = PSolution[indexi+3][indexj];
				
	iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj], \
                                                             Jm[indexj]));

	iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                               Cm[indexj], \
                                                               Bm[indexj], \
                                                               Dm[indexj], \
                                                               Jm[indexj]));

	iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                             Im[indexj], \
                                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj]));

	iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                               Im[indexj], \
                                                               Cm[indexj], \
                                                               Bm[indexj], \
                                                               Dm[indexj]));

				
    } // end for indexj
				
    ho_nS_reac_data->tfac = 0.5;
    aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
    ho_nS_reac_data->tfac = 0.333;
    aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
    ho_nS_reac_data->tfac = 0.5;
    aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
    ho_nS_reac_data->tfac = 0.333;
    aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
	
		
    } else if (indexi==1){
		
      for (indexj=0;indexj<nsv;indexj++){
	if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
	
	Jm[indexj] = Input[indexj];
	Dm[indexj] = PSolution[indexi-1][indexj];
	Bm[indexj] = PSolution[indexi][indexj];
	Cm[indexj] = PSolution[indexi+1][indexj];
	Im[indexj] = PSolution[indexi+2][indexj];
	Km[indexj] = PSolution[indexi+3][indexj];
				
	iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
		                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj], \
                                                             Jm[indexj]));

	iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                               Cm[indexj], \
                                                               Bm[indexj], \
                                                               Dm[indexj], \
                                                               Jm[indexj]));

	iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                             Im[indexj], \
                                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj]));

	iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                               Im[indexj], \
                                                               Cm[indexj], \
                                                               Bm[indexj], \
                                                               Dm[indexj]));
				
	} // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
		  
	} else if (indexi==ngrids-3){
		
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = PSolution[indexi-2][indexj];
	    Dm[indexj] = PSolution[indexi-1][indexj];
	    Bm[indexj] = PSolution[indexi][indexj];
	    Cm[indexj] = PSolution[indexi+1][indexj];
	    Im[indexj] = PSolution[indexi+2][indexj];
	    Km[indexj] = fmax(2*Im[indexj] - Cm[indexj],0);
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

            iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));
				
            iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

            iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj],
                                                                   Im[indexj],
                                                                   Cm[indexj],
                                                                   Bm[indexj],
                                                                   Dm[indexj]));
				
            } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
		  
		
	} else if (indexi==ngrids-2){
			
		
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = PSolution[indexi-2][indexj];
	    Dm[indexj] = PSolution[indexi-1][indexj];
	    Bm[indexj] = PSolution[indexi][indexj];
	    Cm[indexj] = PSolution[indexi+1][indexj];
	    Im[indexj] = fmax(2*Cm[indexj] - Bm[indexj],0);
	    Km[indexj] = fmax(3*Cm[indexj] - 2*Bm[indexj],0);
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
			                                           Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj],\
                                                                   Im[indexj],\
                                                                   Cm[indexj],\
                                                                   Bm[indexj],\
                                                                   Dm[indexj]));
				
	   } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
		
	} else {
	
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = PSolution[indexi-2][indexj];
	    Dm[indexj] = PSolution[indexi-1][indexj];
	    Bm[indexj] = PSolution[indexi][indexj];
	    Cm[indexj] = PSolution[indexi+1][indexj];
	    Im[indexj] = PSolution[indexi+2][indexj];
	    Km[indexj] = PSolution[indexi+3][indexj];
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
	  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
		 
	} // end if
		 
		
    aflag = Liotta_num_method(Am,Bm,Cm,Dm,Em,Fm,Gm,Hm,Im,Jm,Km,ho_nS_reac_data);
	
	for(indexj=0;indexj<nsv;indexj++){
		iSolution[indexi][indexj] =Am[indexj];
	} 
	
} // end for index i
//------------------------------------------------------------------------------
//Part 2:

for(indexi=0;indexi<ngrids-1;indexi++){
  ho_nS_reac_data->LocX = indexi+1*delx/lcolon;
	
  if (indexi==0){
		
    for (indexj=0;indexj<nsv;indexj++){
      if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
      Jm[indexj] = Input[indexj];
      Dm[indexj] = Input[indexj];
      Bm[indexj] = iSolution[indexi][indexj];
      Cm[indexj] = iSolution[indexi+1][indexj];
      Im[indexj] = iSolution[indexi+2][indexj];
      Km[indexj] = iSolution[indexi+3][indexj];
				
      iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                           Cm[indexj], \
                                                           Bm[indexj], \
                                                           Dm[indexj], \
                                                           Jm[indexj]));

      iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj], \
                                                             Jm[indexj]));

      iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                           Im[indexj], \
                                                           Cm[indexj], \
                                                           Bm[indexj], \
                                                           Dm[indexj]));

      iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                             Im[indexj], \
                                                             Cm[indexj], \
                                                             Bm[indexj], \
                                                             Dm[indexj]));
				
        } // end for indexj
				
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
				
	} else if (indexi==1){
		
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = Input[indexj];
	    Dm[indexj] = iSolution[indexi-1][indexj];
	    Bm[indexj] = iSolution[indexi][indexj];
	    Cm[indexj] = iSolution[indexi+1][indexj];
	    Im[indexj] = iSolution[indexi+2][indexj];
	    Km[indexj] = iSolution[indexi+3][indexj];
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
		  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
				
	} else if (indexi==ngrids-4){
	
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = iSolution[indexi-2][indexj];
	    Dm[indexj] = iSolution[indexi-1][indexj];
	    Bm[indexj] = iSolution[indexi][indexj];
	    Cm[indexj] = iSolution[indexi+1][indexj];
	    Im[indexj] = iSolution[indexi+2][indexj];
	    Km[indexj] = fmax(2*Im[indexj] - Cm[indexj],0);
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
	  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);

	} else if (indexi==ngrids-3){
		
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = iSolution[indexi-2][indexj];
	    Dm[indexj] = iSolution[indexi-1][indexj];
	    Bm[indexj] = iSolution[indexi][indexj];
	    Cm[indexj] = iSolution[indexi+1][indexj];
	    Im[indexj] = fmax(2*Cm[indexj] - Bm[indexj],0);
	    Km[indexj] = fmax(2*Im[indexj] - Cm[indexj],0);
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
	  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
				
	} else if (indexi==ngrids-2){
			
	
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = iSolution[indexi-2][indexj];
	    Dm[indexj] = iSolution[indexi-1][indexj];
	    Bm[indexj] = iSolution[indexi][indexj];
	    Cm[indexj] = fmax(2*Bm[indexj] - Dm[indexj],0);
	    Im[indexj] = fmax(2*Cm[indexj] - Bm[indexj],0);
	    Km[indexj] = fmax(3*Cm[indexj] - 2*Bm[indexj],0);
	 			
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
	  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
				
	} else {
		
	  for (indexj=0;indexj<nsv;indexj++){
	    if(indexj%2==0){vl = fr_mpd;} else {vl = 0;}
				
	    Jm[indexj] = iSolution[indexi-2][indexj];
	    Dm[indexj] = iSolution[indexi-1][indexj];
	    Bm[indexj] = iSolution[indexi][indexj];
	    Cm[indexj] = iSolution[indexi+1][indexj];
	    Im[indexj] = iSolution[indexi+2][indexj];
	    Km[indexj] = iSolution[indexi+3][indexj];
				
	    iEm[indexj] = Bm[indexj] - 0.5*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj], \
                                                                 Jm[indexj]));

	    iGm[indexj] = Bm[indexj] - 0.333*vl*gamma*(FluxLimiter(Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj], \
                                                                   Jm[indexj]));

	    iFm[indexj] = Cm[indexj] - 0.5*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                 Im[indexj], \
                                                                 Cm[indexj], \
                                                                 Bm[indexj], \
                                                                 Dm[indexj]));

	    iHm[indexj] = Cm[indexj] - 0.333*vl*gamma*(FluxLimiter(Km[indexj], \
                                                                   Im[indexj], \
                                                                   Cm[indexj], \
                                                                   Bm[indexj], \
                                                                   Dm[indexj]));
				
	  } // end for indexj
		  
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Em,iEm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Gm,iGm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.5;
	aflag = Non_lin_solver(Fm,iFm,nsv,ho_nS_reac_data);
	ho_nS_reac_data->tfac = 0.333;
	aflag = Non_lin_solver(Hm,iHm,nsv,ho_nS_reac_data);
				

		 
	} // end if
		 
		
    aflag = Liotta_num_method(Am,Bm,Cm,Dm,Em,Fm,Gm,Hm,Im,Jm,Km,ho_nS_reac_data);

	for(indexj=0;indexj<nsv;indexj++){
		PSolution[indexi+1][indexj] =Am[indexj];
	} 
} // end for index i

//------------------------------------------------------------------------------
// Interpolate staggared solution back to the first x location	
	
for (indexj=0;indexj<nsv;indexj++){
  if (indexj%2 == 0){
    PSolution[0][indexj] =2./3*Input[indexj] + 2./6*PSolution[1][indexj]; 
  } else {
    PSolution[0][indexj] = 2*PSolution[1][indexj] - PSolution[2][indexj];
  }
}

//------------------------------------------------------------------------------
// free memory

free(Am);
free(Bm);
free(Cm);
free(Dm);
free(Em);
free(Fm);
free(Gm);
free(Hm);
free(Im);
free(Jm);
free(Km);
free(iEm);
free(iFm);
free(iGm);
free(iHm);

	for(indexi=0;indexi<ngrids-1;indexi++){
		free(iSolution[indexi]);
	}
	free(iSolution);

	
return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int LoadMealPlan(double** mealPlan)
{
	int 			  indexi,indexj;
	int 			  aflag;
	const int 	bsz = 300; 
	char 			  buf[bsz];
	//~ 
	double* u_data_time = malloc(MMD * sizeof (double)); 
	aflag = Clean1DArray(u_data_time, MMD);
	double* u_data_prebiotic = malloc(MMD * sizeof (double));
	aflag = Clean1DArray(u_data_prebiotic, MMD);
	double* u_data_probiotic = malloc(MMD * sizeof (double));
  aflag = Clean1DArray(u_data_probiotic, MMD);
  //~ 
	FILE *f_umeal_inputs;	
	f_umeal_inputs = fopen("InputFiles/user_diet_compuGUT.txt","r");
	
	  if (f_umeal_inputs != NULL){
		indexi = 0;
	    while (fgets(buf, bsz, f_umeal_inputs) != NULL) {
	      if (buf[0] != '#' && buf[0]!= '\n') {
		  sscanf(buf,"%lf %lf %lf", &u_data_time[indexi], \
		                            &u_data_prebiotic[indexi], \
		                            &u_data_probiotic[indexi]);
		  indexi++;
	      }
	    }	
	  } else {
	    perror("fopen");
	  }
	  
	fclose(f_umeal_inputs);
	g_MealData = indexi;
	// printf("%d\n",g_MealData);
	//~ 
	for (indexi=0;indexi<g_MealData;indexi++){
		mealPlan[indexi][0] = u_data_time[indexi]/24;
		mealPlan[indexi][1] = u_data_prebiotic[indexi];
		mealPlan[indexi][2] = u_data_probiotic[indexi];
	}
	//~ 
	free(u_data_time);
	free(u_data_prebiotic);
	free(u_data_probiotic);	
	
	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int PreGut(double* input, double** PSol, double** mealPlan, 
  double time, double dt)
{
	double time_h = 24*time;
	//printf("%f\n",time_h);
	
	int indexi,indexj,indexk;
	
	int sdiv = g_sdiv;
	
	double q = g_fr_lpd;
	double V = g_vsintest_l/sdiv;
	double D = q/V;
	
	double Fin;
	double Bin;
	double mealOn = time_h-floor(time_h);
	double svol = 1.5; // volume of stomach in litres
	
	indexi = floor(time_h);
	//~ // printf("%d\n",indexi);
	
	if (mealOn < MEALLENGTH){
		Fin = mealPlan[indexi][1]/q/MEALLENGTH*24;
		Bin = mealPlan[indexi][2]/q/MEALLENGTH*24/g_nsd;
	} else {
		Fin = 0;
	}
	
	double Fprime;
	double Bprime;
	double InF;
	double InB;
	double OutF;
	double OutB;
	
	// printf("%f %f %f %f\n",mealPlan[indexi][1],q,MEALLENGTH,Fin);
	// PauseMessage(4371);

	for(indexi=0;indexi<sdiv;indexi++){
	  if(indexi==0){
			InF = Fin;
			InB = Bin;
			OutF = PSol[indexi][18];
			OutB = PSol[indexi][20];
		} else {
			InF = PSol[indexi-1][18];
			InB = PSol[indexi-1][20];
			OutF = PSol[indexi][18];
			OutB = PSol[indexi][20];
		}
		Fprime = D*(InF-OutF);
		Bprime = D*(InB-OutB);
		PSol[indexi][18] = PSol[indexi][18] + dt/sdiv*Fprime;
		PSol[indexi][20] = PSol[indexi][20] + dt/sdiv*Bprime;
	}
		
	input[18] = PSol[sdiv-1][18];
	
	indexj=20;
	for(indexi=0;indexi<sdiv;indexi++){
		for(indexk=0;indexk<g_nsd-1;indexk++){
		indexj=indexj+2;
		PSol[indexi][indexj] = PSol[indexi][20];
		}
	}
		
	indexj=20;	
	for(indexk=0;indexk<g_nsd;indexk++){
		input[indexj] = PSol[sdiv-1][indexj];
		indexj=indexj+2;

	}	

		
	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
static int PGUT_interpolator(double** mealPlan, double mealTime, 
  double SimTime, double** OutputMeal, int MaxIter, double dt)
{
	int indexi,indexj,indexk,aflag;
	int sdiv = g_sdiv;
	
	double q = g_fr_lpd;
	double V = g_vsintest_l/sdiv;
	double D = q/V;
	
	double Vstomach = 1.5; // Volume of Stomach in L
	
	//~ printf("Small Intestine dilution rate: %f\n",D);
	
	double mealLength = MEALLENGTH/24;
	int nsteps = (int)SimTime*1000;
  double dte = SimTime/nsteps;
  double dta = 1.0/24.0;
  
  double * time = malloc(nsteps * sizeof (double));
  aflag = Clean1DArray(time, nsteps);
  double * CI = malloc(nsteps * sizeof (double));
  aflag = Clean1DArray(CI, nsteps);
  double * e1 = malloc(nsteps * sizeof (double));
  aflag = Clean1DArray(e1, nsteps);	
  double * e2 = malloc(nsteps * sizeof (double));
  aflag = Clean1DArray(e2, nsteps);
  double * e3 = malloc(nsteps * sizeof (double));
  aflag = Clean1DArray(e3, nsteps);
  
  
  double iTime;
  double d1;
  double d2;
  double d3;
  
  double teststat;
  double Fin;
  
  
  FILE *f7;
  f7=fopen("PreColonData.txt","w");
  
  for(indexi=0;indexi<nsteps;indexi++){
	  time[indexi] = dte*indexi;
	  iTime = fmod(time[indexi], mealTime);
	  //~ printf("%f\n",iTime);
	  	  
	  teststat = fmod(iTime,(dta+EPSI));
	  //~ printf("%f\n",teststat);
	  //~ PauseMessage(4831);
	  
	  if (teststat <= mealLength){
		  CI[indexi] = mealPlan[(int)(24.0*iTime+EPSI)][1];
		} else {
			CI[indexi] = 0;
		}
		
		Fin = (CI[indexi]/mealLength)/q;
		
    if (indexi==0){
			d1 = D*Fin;
			d2 = 0;
			d3 = 0;
		} else {
			d1 = D*(Fin-e1[indexi-1]);
			d2 = D*(e1[indexi-1]-e2[indexi-1]);
			d3 = D*(e2[indexi-1]-e3[indexi-1]);
		
		  e1[indexi] = e1[indexi-1]+dte*d1;
		  e2[indexi] = e2[indexi-1]+dte*d2;
		  e3[indexi] = e3[indexi-1]+dte*d3;
		}
		
		fprintf(f7,"%f\t%f\t%f\t%f\t%f\n",time[indexi], \
		                                  Fin, \
		                                  e1[indexi], \
		                                  e2[indexi], \
		                                  e3[indexi]);
		
	  
	}
  fclose(f7);
  
  
  free(time);
  free(CI);
  free(e1);
  free(e2);
  
  indexj=0;
  for(indexi=0;indexi<nsteps;indexi++){
		double d1 = indexi*dte;
		double d2 = indexj*dt;
		
		if(d2-d1<dte){
		  if(indexi==0){	
		    OutputMeal[indexj][18] = e3[indexi];
			} else {
			  OutputMeal[indexj][18] = (d2-d1)*(e3[indexi+1]-e3[indexi])/dte + e3[indexi];
			}
		indexj=indexj+1;
		}
		
	}
  
  //~ for(indexi=0;indexi<MaxIter;indexi++){
		//~ printf("%f\n",OutputMeal[indexi][18]);
	//~ }
  free(e3);
	
	return(0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// END OF FILE -----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

