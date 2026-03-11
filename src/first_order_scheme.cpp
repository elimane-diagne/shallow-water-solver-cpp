#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

// Constants
const std::size_t N = 200;         // Number of spatial points
const double dx = 0.1;             // Spatial step
const double g = 9.81;             // Gravity
const double cfl = 1.0;            // CFL number
const int n_steps = 3;             // Number of time steps
const double h_eps = 1e-6;         // Small value to avoid division by zero

// Initialize the physical fields
void initialize(std::vector<double>& h, std::vector<double>& u, std::vector<double>& z) {
    for (std::size_t i = 0; i < N; ++i) {
        z[i] = 0.5;
        u[i] = 0.0;
        h[i] = 1.0 - z[i]; // water height adjusted to keep eta = 1.0
    }
}

// Compute maximum wave speed (for CFL condition)
double compute_max_wave_speed(const std::vector<double>& h, const std::vector<double>& u) {
    double s_max = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        double c = std::sqrt(g * std::max(h[i], h_eps));
        s_max = std::max(s_max, std::abs(u[i]) + c);
    }
    return s_max;
}

// Hydrostatic reconstruction + Rusanov flux + source term
void update_solution(const std::vector<double>& h, const std::vector<double>& u,
                     const std::vector<double>& z,
                     std::vector<double>& h_new, std::vector<double>& u_new,
                     double dt) {
    // Solid wall boundaries (ghost cells not implemented)
    h_new[0]     = h[0];
    u_new[0]     = u[0];
    h_new[N - 1] = h[N - 1];
    u_new[N - 1] = u[N - 1];

    for (std::size_t i = 1; i < N - 1; ++i) {
        // Hydrostatic reconstruction at i+1/2
        double zmax_plus = std::max(z[i], z[i + 1]);
        double hL_plus = std::max(0.0, h[i]     + z[i]     - zmax_plus);
        double hR_plus = std::max(0.0, h[i + 1] + z[i + 1] - zmax_plus);

        // Hydrostatic reconstruction at i-1/2
        double zmax_minus = std::max(z[i], z[i - 1]);
        double hL_minus = std::max(0.0, h[i - 1] + z[i - 1] - zmax_minus);
        double hR_minus = std::max(0.0, h[i]     + z[i]     - zmax_minus);

        // Velocities at interfaces (simple upwind)
        double uL_plus = u[i];
        double uR_plus = u[i + 1];
        double uL_minus = u[i - 1];
        double uR_minus = u[i];

        // Source terms
        double S_plus  = g * 0.5 * (hL_plus  + hR_plus)  * (z[i + 1] - z[i]);
        double S_minus = g * 0.5 * (hL_minus + hR_minus) * (z[i] - z[i - 1]);

        // Fluxes at i+1/2
        double lambda_plus = std::max(std::abs(uL_plus) + std::sqrt(g * hL_plus),
                                      std::abs(uR_plus) + std::sqrt(g * hR_plus));

        double Fh_plus = 0.5 * (hL_plus * uL_plus + hR_plus * uR_plus)
                       - 0.5 * lambda_plus * (hR_plus - hL_plus);

        double Fhu_plus = 0.5 * (hL_plus * uL_plus * uL_plus + 0.5 * g * hL_plus * hL_plus +
                                 hR_plus * uR_plus * uR_plus + 0.5 * g * hR_plus * hR_plus)
                        - 0.5 * lambda_plus * (hR_plus * uR_plus - hL_plus * uL_plus);

        // Fluxes at i-1/2
        double lambda_minus = std::max(std::abs(uL_minus) + std::sqrt(g * hL_minus),
                                       std::abs(uR_minus) + std::sqrt(g * hR_minus));

        double Fh_minus = 0.5 * (hL_minus * uL_minus + hR_minus * uR_minus)
                        - 0.5 * lambda_minus * (hR_minus - hL_minus);

        double Fhu_minus = 0.5 * (hL_minus * uL_minus * uL_minus + 0.5 * g * hL_minus * hL_minus +
                                  hR_minus * uR_minus * uR_minus + 0.5 * g * hR_minus * hR_minus)
                         - 0.5 * lambda_minus * (hR_minus * uR_minus - hL_minus * uL_minus);

        // Add source terms
        Fhu_plus  += S_plus;
        Fhu_minus += S_minus;

        // Update conserved variables
        h_new[i] = h[i] - (dt / dx) * (Fh_plus - Fh_minus);

        if (h_new[i] > h_eps) {
            u_new[i] = ((h[i] * u[i]) - (dt / dx) * (Fhu_plus - Fhu_minus)) / h_new[i];
        } else {
            u_new[i] = 0.0;
        }
    }
}

// Export results to file
void save_to_file(const std::vector<double>& h, const std::vector<double>& u,
                  const std::vector<double>& z, int step) {
    std::string filename = "results_step_" + std::to_string(step) + ".txt";
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return;
    }

    file << "# x\tz(x)\th(x)\teta(x)\tu(x)\n";
    for (std::size_t i = 0; i < N; ++i) {
        double x = i * dx;
        file << x << "\t" << z[i] << "\t" << h[i] << "\t"
             << h[i] + z[i] << "\t" << u[i] << "\n";
    }

    file.close();
}

// Main
int main() {
    std::vector<double> h(N), u(N), z(N);
    std::vector<double> h_new(N), u_new(N);

    initialize(h, u, z);
    save_to_file(h, u, z, 0);

    for (int step = 1; step <= n_steps; ++step) {
        double s_max = compute_max_wave_speed(h, u);
        double dt = cfl * dx / s_max;

        update_solution(h, u, z, h_new, u_new, dt);

        h = h_new;
        u = u_new;

        save_to_file(h, u, z, step);
    }

    std::cout << "Simulation completed. Results saved in text files.\n";
    return 0;
}
