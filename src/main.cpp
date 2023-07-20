// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <boost/filesystem.hpp>
#include <cvode/cvode.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

#include <2DTissue.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

// The ODE system
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    realtype yval = NV_Ith_S(y, 0);
    NV_Ith_S(ydot, 0) = -0.5 * yval;
    return 0;
}


int main()
{
    int step_count = 30;
    bool save_data = false;

    // Path to the 3D mesh file
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";
    // std::string mesh_path = PROJECT_PATH.string() + "/meshes/sphere.off";

    // for (int particle_count = 200; particle_count <= 200; particle_count += 100) {
    //     _2DTissue _2dtissue(save_data, mesh_path, particle_count, step_count, 0.01);  // Initialize the 2DTissue object

    //     _2dtissue.start();

    //     std::clock_t start = std::clock();

    //     while(!_2dtissue.is_finished()) {
    //         System data = _2dtissue.update();
    //     }
    //     std::cout << _2dtissue.get_order_parameter() << '\n';

    //     std::clock_t end = std::clock();
    //     double duration = (end - start) / (double) CLOCKS_PER_SEC;
    //     std::cout << "Time taken: " << duration << " seconds" << '\n';
    // }

    /*
    Sundials test code
    */
    SUNContext sunctx;
    void* comm = nullptr;
    int flag = SUNContext_Create(comm, &sunctx);
    if (flag != 0) {
        // handle error
    }

    N_Vector y = N_VNew_Serial(1, sunctx); 
    realtype t = 0.0;
    realtype reltol = 1e-4, abstol = 1e-4;  // Set tolerances

    // Initialize y
    NV_Ith_S(y, 0) = 1.0;

    /*
    For more information on the following functions, visit
    https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#user-callable-functions#
    */
    // Call CVodeCreate to create the solver memory and specify the
    // The function CVodeCreate() instantiates a CVODE solver object and specifies the solution method.
    // CV_BDF for stiff problems.
    // CV_ADAMS for nonstiff problems 
    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);

    // Initialize the integrator memory and specify the user's right hand
    // side function in y'=f(t,y), the inital time T0, and the initial
    // dependent variable vector y.
    CVodeInit(cvode_mem, f, t, y);

    // Set the scalar relative tolerance and scalar absolute tolerance
    // cvode_mem – pointer to the CVODE memory block returned by CVodeCreate()
    CVodeSStolerances(cvode_mem, reltol, abstol);

    // Call CVDense to specify the CVDENSE dense linear solver
    SUNMatrix A = SUNDenseMatrix(1, 1, sunctx); // Create dense matrix for use in linear solves
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx); // Create dense linear solver for use in linear solves
    CVodeSetLinearSolver(cvode_mem, LS, A); // Attach the matrix and linear solver

    // Integrate over the interval, while the test function returns success
    realtype tout = 1.0;
    while (t < 10.0)
    {
        int flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (flag < 0)
        {
            std::cerr << "Error in integration" << std::endl;
            return -1;
        }
        std::cout << "At t = " << t << ", y = " << NV_Ith_S(y,0) << std::endl;
        tout += 1.0;
    }

    // Free integrator memory
    CVodeFree(&cvode_mem);
    SUNContext_Free(&sunctx);

    // Free vector
    N_VDestroy(y);

    return 0;
}
