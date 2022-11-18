#include "fitter_arma.hpp"

arma::vec centroid (const atom_group& group)
{
    arma::vec centroid = arma::zeros<arma::vec>(3);
    //sum of all the atoms coordinates:
    for (auto& atom : group.atoms) {
        centroid[0] += atom.position[0];
        centroid[1] += atom.position[1];
        centroid[2] += atom.position[2];
    }

    //divided by the sum of weights = number of atoms;
    centroid[0] /= group.atoms.size();
    centroid[1] /= group.atoms.size();
    centroid[2] /= group.atoms.size();

    return centroid;
}

void svd_fitter::fit_base (atom_group& group, const atom_group& ideal_base)
{

    arma::vec group_centroid = centroid(group);
    arma::vec ideal_base_centroid = centroid(ideal_base);

    arma::mat X = arma::mat(3, group.atoms.size());
    arma::mat Y = arma::mat(3, ideal_base.atoms.size());

    for (int i = 0; i < group.atoms.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            X(j,i) = group.atoms[i].position[j] - group_centroid[j];
    }

    for (int i = 0; i < ideal_base.atoms.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            Y(j,i) = ideal_base.atoms[i].position[j] - ideal_base_centroid[j];
    }
    arma::mat S = X*Y.t();
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U,s,V, S);
    arma::mat sigma = arma::eye(3,3);
    sigma(2,2) = arma::det(V*U.t());
    arma::mat R = V*sigma*U.t();
// The commented code is for curveplus fitting algo. 
/*    for (int i=0; i<group.atoms.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            X(j,i) = group.atoms[i].position[j] ;
    }

    for (int i=0; i<ideal_base.atoms.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            Y(j,i) = ideal_base.atoms[i].position[j];
    }
    double D = arma::det(S);
    arma::mat OHM = arma::mat(6,6,arma::fill::zeros);
    for (int i = 0 ; i < 3 ; i++) {
	for (int j = 3; j<6; j++) {
		OHM(i,j) = S(i,j-3);
		OHM(j,i) = S(i,j-3);
		}
	}
   arma::vec eigval ; 
   arma::mat eigvec ; 
   eig_sym(eigval,eigvec,OHM);
   arma::mat eig_vecr = arma::mat(6,6);
   eig_vecr = arma::real(eigvec);
   arma::vec h1 = arma::zeros<arma::vec>(3);
   arma::vec h2 = arma::zeros<arma::vec>(3);
   arma::vec h3 = arma::zeros<arma::vec>(3);
   arma::vec k1 = arma::zeros<arma::vec>(3);
   arma::vec k2 = arma::zeros<arma::vec>(3);
   arma::vec k3 = arma::zeros<arma::vec>(3);
   arma::mat R = arma::mat(3,3);
   arma::vec Sqr = arma::zeros<arma::vec>(1);
   double m = 1.41421356237 ; 
   for (int i=0; i<3; i++){
	h1[i] = m*eig_vecr(i,5);
	h2[i] = m*eig_vecr(i,4);
	k1[i] = m*eig_vecr(i+3,5);
	k2[i] = m*eig_vecr(i+3,4);
   	}
   h3 = cross(h1,h2);
   k3 = cross(k1,k2);
   R = k1*h1.t() + k2*h2.t() + k3*h3.t();
*/
   group.r =group_centroid - R.t()*ideal_base_centroid;
   group.R = R;
}
