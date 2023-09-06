/**
 * @file        Cell.cpp
 * @brief       Simulating the Eigenleben of a cell
 *
 * @author      Jan-Piotraschke
 * @date        2023-Aug-29
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */

#include <Cell.h>

Cell::Cell()
    : reltol(1e-4),
      abstol(1e-4),
      t(0.0),
      tout(0.001),
      y(N_VNew_Serial(NEQ)),
      cvode_mem(CVodeCreate(CV_BDF)),  // CVODE solver: CV_BDF for stiff problems; CV_ADAMS for nonstiff problems
      A(SUNDenseMatrix(NEQ, NEQ)),
      startTime(0.0),
      endTime(10.0),
      numberOfPoints(101)
{
    // Initialize y
    NV_Ith_S(y, 0) = 0.0; // y(0) = 0
    NV_Ith_S(y, 1) = 1.0; // y'(0) = 1

    /*
    For more information on the following functions, visit
    https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#user-callable-functions#
    */

    // Initialize the integrator memory and specify the user's right hand
    // side function in y'=f(t,y), the inital time T0, and the initial
    // dependent variable vector y.
    CVodeInit(cvode_mem, simulate_sine, t, y);

    // Set the scalar relative tolerance and scalar absolute tolerance
    // cvode_mem – pointer to the CVODE memory block returned by CVodeCreate()
    CVodeSStolerances(cvode_mem, reltol, abstol);

    LS = SUNLinSol_Dense(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    // Initialize SBML model simulation parameters.
    sbmlModelFilePath = PROJECT_PATH + "/sbml-model/BIOMD0000000613_url.xml";
}

Cell::~Cell() {
    free_memory();
}


// ========================================
// Public Functions
// ========================================

double Cell::update(realtype tout){
    CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    // Update v0 by the tenth of the sine wave which oscillates between 0 and 1
    double v0 = 0.1 * (0.5 * (1 + NV_Ith_S(y, 0)));

    return v0;
}


void Cell::perform_sbml_simulation() {
    rr = new rr::RoadRunner();

    // Load the SBML model.
    rr->load(sbmlModelFilePath);

    // Set up the integrator.
    rr->getIntegrator()->setValue("relative_tolerance", 1e-6);
    rr->getIntegrator()->setValue("absolute_tolerance", 1e-6);

    // Simulate the model.
    rr::SimulateOptions options;
    options.start = startTime;
    options.duration = endTime - startTime;
    options.steps = numberOfPoints - 1;

    // Print the result of the simulation.
    std::cout << *rr->simulate(&options) << std::endl;

    // Don't forget to free the memory.
    delete rr;
}


void Cell::free_memory() {
    // Free memory
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS); /* Free the linear solver memory */
    SUNMatDestroy(A); /* Free the matrix memory */
    N_VDestroy(y);
}



// ========================================
// Private Functions
// ========================================

// The ODE system
int Cell::simulate_sine(
    realtype t,
    N_Vector y,
    N_Vector ydot,
    void *user_data
) {
    realtype sine = NV_Ith_S(y, 0);
    realtype cose = NV_Ith_S(y, 1);

    // Store the data in the ydot vector
    NV_Ith_S(ydot, 0) = cose;
    NV_Ith_S(ydot, 1) = -sine;

    return 0;
}
