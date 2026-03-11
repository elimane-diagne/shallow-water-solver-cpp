Simulation of the Shallow Water Equations using the Finite Volume Method (C++)

This project implements a **numerical solver for the 1D Shallow Water Equations** using the **Finite Volume Method**.
The objective is to simulate shallow water flows while preserving important physical properties such as **mass conservation**, **positivity of the water height**, and **steady states (lake at rest)**.

The solver includes a first-order finite volume scheme, hydrostatic reconstruction, and slope limiting techniques to improve numerical stability.

This type of model is widely used in **hydraulics, oceanography, and environmental fluid dynamics**.


## Main Features

 Numerical simulation of the **1D Shallow Water Equations**
   Finite Volume discretization**
   First-order numerical scheme**
   Slope limiter** for numerical stability
   Hydrostatic reconstruction** for well-balanced property
   Modular C++ implementation

## Mathematical Model

The Shallow Water Equations describe the evolution of a thin layer of fluid under gravity.

Continuity equation:

∂ₜ h + ∂ₓ (hu) = 0

Momentum equation:

∂ₜ (hu) + ∂ₓ (hu² + g h² / 2) = − g h ∂ₓ z


where:

* (h(x,t)) : water height
* (u(x,t)) : velocity
* (z(x)) : bottom topography
* (g) : gravitational acceleration

---

## Numerical Method

The numerical solution is obtained using the **Finite Volume Method**.

### Spatial discretization

The computational domain is divided into **control volumes**.
The conservative variables are updated using **numerical fluxes at cell interfaces**.

### Reconstruction

To improve accuracy and maintain stability, the solver includes:

Hydrostatic reconstruction**
Slope limiter

These techniques ensure:

preservation of steady states
avoidance of non-physical oscillations

## Project Structure

```
shallow-water-cpp
│
├── src
│   ├── main.cpp
│   ├── initialize.cpp
│   ├── reconstruction.cpp
│   ├── first_order_scheme.cpp
│   └── slope_limiter.cpp
│
├── include
│   └── simulation.hpp
│
├── README.md
└── .gitignore
```

Description of files:

* **main.cpp** : program entry point
* **initialize.cpp** : initialization of variables and domain
* **reconstruction.cpp** : hydrostatic reconstruction
* **first_order_scheme.cpp** : finite volume scheme implementation
* **slope_limiter.cpp** : slope limiting for stability
* **simulation.hpp** : simulation structures and declarations

---

## Requirements

The project requires:

* **C++ compiler** (g++ recommended)
* Standard C++ libraries

Tested with:

* g++ (GNU Compiler)

---

## Compilation and Execution

Clone the repository:

```bash
git clone https://github.com/elimane-diagne/SHALLOW-WATER-Cpp
cd SHALLOW-WATER-Cpp
```

Compile the code:

```bash
g++ src/*.cpp -o simulation
```

Run the simulation:

```bash
./simulation
```

---

## Applications

Shallow water models are widely used for:

* river flow modeling
* flood simulations
* coastal oceanography
* environmental hydraulics
* tsunami propagation studies

---

## Future Improvements

Possible extensions of the project include:

* **Second-order numerical scheme**
* **2D shallow water solver**
* **Visualization tools**
* **Parallel computing (OpenMP or MPI)**
* **Adaptive mesh refinement**

---

## Author

**Elimane Diagne**

Master student in applied mathematics / numerical simulation.

---

## License

This project is intended for **academic and educational purposes**.
