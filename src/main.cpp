#include "simulation.hpp"
#include <iostream>

int main() {
    const std::size_t N = 10;
    const double dx = 0.1;

    std::vector<double> x(N), z(N), h(N), u(N);

    initialize(x, z, u, h, dx);

    // Calcul pente et reconstruction
    std::vector<double> dq, hL, hR;
    compute_slopes(h, dq);
    reconstruct_interfaces(h, dq, hL, hR);

    // Affichage test
    std::cout << "Interface i+1/2 : hL | hR\n";
    for (std::size_t i = 0; i < hL.size(); ++i) {
        std::cout << "i=" << i << " : " << hL[i] << " | " << hR[i] << "\n";
    }

    return 0;
}
