#pragma once

#include <iostream>
#include <vector>
#include <array>

#include <armadillo>

#include "atom_group.hpp"

struct svd_fitter {
    //void fit_base (const std::vector<std::array<float,3>>& data,
                   //const std::vector<std::array<float,3>>& ideal_base,
                   //arma::mat& rotation_out,
                   //arma::vec& translation_out);
    void fit_base (atom_group& group, const atom_group& ideal);
};
