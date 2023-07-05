// author: @Jan-Piotraschke
// date: 2023-05-20
// license: Apache License 2.0
// version: 0.1.0

/**
* BACKGROUND:
*
* Cell-cell interactions (e.g. via cell-cell signaling) are important for the development of multicellular organisms.
*
* Cancer cells escape the cooperative rules.
*
* Interesting are syntethic oncogenes -> they are mutations that are not oncogenic on their own but,
* when they appear with other mutations, they drive tumor formation.
*
* The proliferation behavoir of cancer cells is still poorly understood in part because it is difficult to
* experimentally study the transmisson of Proliferative Factors from one cell to its neighbors.
*
* In this file we want to collect the different types of Cell-Cell Information Fluxes that we want to simulate.
*/

#include <Eigen/Dense>

#include <particle_simulation/cell_cell_interactions.h>


/**
* @brief: Cells can attract or repel each other as they move depending on their distance.
*/
Eigen::Vector3d repulsive_adhesion_motion(
    double k,
    double σ,
    double dist,
    double r_adh,
    double k_adh,
    const Eigen::Vector3d& dist_v
) {
    double Fij_rep = 0;
    double Fij_adh = 0;

    if (dist < 2*σ)
    {
        Fij_rep = (-k * (2 * σ - dist)) / (2 * σ);
    }

    if (dist >= 2*σ && dist <= r_adh)
    {
        Fij_adh = (k_adh * (2 * σ - dist)) / (2 * σ - r_adh);
    }

    double Fij = Fij_rep + Fij_adh;

    return Fij * (dist_v / dist);
}
