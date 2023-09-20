#pragma once

class LocomotionHelperInterface {
public:
    virtual void calculate_forces_between_particles() = 0;
};

class OrientationHelperInterface {
public:
    virtual void calculate_average_n_within_distance() = 0;
};
