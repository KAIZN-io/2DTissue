// Cell.h
#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

// Differential Equation Simulation
#include <rr/rrRoadRunner.h>
#include <rr/rrExecutableModel.h>

#include <cvode/cvode.h>
// #include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

class Cell {
public:
    Cell();
    ~Cell(); // Destructor to manage memory cleanup

    double update(realtype tout);
    void perform_sbml_simulation();
    void free_memory();

private:
    static constexpr int NEQ = 2; // Number of equations
    const std::string PROJECT_PATH = PROJECT_SOURCE_DIR;

    // Differential Equation Simulation
    realtype reltol, abstol; // Tolerances
    realtype t; // Time
    realtype tout = 0.001; // Time for next output
    void* cvode_mem; // CVODE memory
    N_Vector y; // Variables
    SUNMatrix A; // Dense SUNMatrix
    SUNLinearSolver LS; // Dense SUNLinearSolver object

    // SBML simulation
    rr::RoadRunner* rr;
    std::string sbmlModelFilePath;
    double startTime;
    double endTime;
    int numberOfPoints;

    static int simulate_sine(
        realtype t,
        N_Vector y,
        N_Vector ydot,
        void *user_data
    );
};
