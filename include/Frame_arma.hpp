#pragma once

#include <armadillo>

struct Frame {
    std::string name;
    int chain_id;
    std::string residue_name;
    int residue_number;
    arma::mat R;
    arma::vec t;
};

