#pragma once

#include <vector>
#include <Eigen/Dense>

class LocomotionHelperInterface {
public:
    virtual void calculate_forces_between_particles() = 0;
};
