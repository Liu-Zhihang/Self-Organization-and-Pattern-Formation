# Self-Organization and Pattern Formation

These are self-study notes for the course [**Self-Organization and Pattern Formation, Prof. Erwin Frey, LMU Munich, Winter Semester 2025/2026**](https://www.theorie.physik.uni-muenchen.de/lsfrey/teaching/winter_term_2025_2026/index.html). Prof. Frey prefers chalkboard lectures. To document my learning process, I've structured my notes into articles and written [Python code](https://github.com/Liu-Zhihang/Self-Organization-and-Pattern-Formation/tree/main/code) to deepen my understanding.

![Course Overview](cn/assets/images/01_001_853e3a67-9126-4927-907c-9f48eade8561.png)

**Course Playlist:** [YouTube Playlist](https://www.youtube.com/watch?v=apMmakW_090&list=PL2IEUF-u3gRfIw2cHFPqNHjNfQL5BxG5m)

**Official Course Link:** [LMU Munich - Self-Organization and Pattern Formation](https://www.theorie.physik.uni-muenchen.de/lsfrey/teaching/winter_term_2025_2026/index.html)


## Course Overview

Can we build predictive theories for self-organizing systems—spanning fluids, reaction–diffusion systems, soft/active matter, and living systems—without tracking every microscopic detail? The answer is yes: by using universality, symmetry, and conservation laws, and by adopting coarse-graining methods, we can reveal instability mechanisms, mesoscale physical laws, and scaling behavior.

The course covers spatially extended nonequilibrium systems that spontaneously generate order. It focuses on the onset of instabilities, pattern selection, dynamical processes, and the emergence of mesoscale physical laws. Methods include dynamical systems, continuum field theory, statistical physics, weakly nonlinear analysis, and differential geometry.

## Course Topics

### Dynamical systems & bifurcations
Linear/nonlinear stability, fixed points/attractors, excitability, synchronization, chaos; normal forms and amplitude equations.

- 1\. Self-Organization: From Bird Flocks to Living Cells
- 2\. Dynamical Systems: Phase Space, Stability, Ginzburg-Landau, Evolutionary Game Theory
- 3\. Dynamical Systems: Fixed Points, Saddle-Node & Pitchfork Bifurcations
- 4\. Dynamical Systems: Cusp Bifurcation, Hysteresis, Infection Model & Transcritical Bifurcation
- 5\. Dynamical Systems: Mass Conserving 2D, Turing Patterns, Stability Analysis, Jacobian Matrix
- 6\. Dynamical Systems: 2D Non-Conservative Systems & Oscillators
- 7\. Dynamical Systems: Poincaré-Bendixson Theorem, RhoGTPase Oscillator & Rock-Paper-Scissors
- 8\. Dynamical Systems: Replicator Dynamics & Spatial Extension, Lyapunov Function, Diffusion Equation

### Pattern formation
Phase separation and coarsening; interface dynamics; reaction–diffusion patterns and waves; excitable media; fronts, pulses, dendrites; thin-film dewetting and flow.

- 9\. Relaxational Dynamics: Ginzburg-Landau Model, Model A, Interface Dynamics & Surface Tension
- 10\. Relaxational Dynamics: Interface Dynamics & Curvature-Driven Flow
- 11\. Thermodynamics and Phase Separation in Liquid Mixtures, Maxwell Construction, Phase Diagrams
- 12\. Dynamics of Liquid Mixtures, Osmotic Pressure, Cahn-Hilliard Equation
- 13\. The Lattice Gas Model & Mean-Field Theory, Mixing Entropy, Flory-Huggins Theory
- 14\. Cahn-Hilliard Equation & Phase Separation Dynamics, Dispersion Relation, Fastest Growing Mode
- 15\. Nucleation Theory & Interface Dynamics, Gibbs-Thomson, Growth and Shrinkage of a Single Droplet
- 16\. Ostwald Ripening & Coarsening Dynamics, Lifshitz-Slyozov-Wagner (LSW) Theory

### Active field theories & broken detailed balance
Non-variational couplings and stresses; absence of free energy/Lyapunov functional; deterministic currents; universality and scaling near nonequilibrium phase transitions.

- 17\. Chemotaxis & Active Motion of Bacteria, Run-and-Tumble Model
- 18\. Keller-Segel Model, Chemotactic Instability, Finite-Time Blowup
- 19\. Scalar Active Matter, Motility-Induced Phase Separation, Active Model B
- 20\. Active Brownian Particles, Active Model B+, Bubbling Phase

### Instability mechanisms
Turing and Hopf instabilities; interface instabilities (Mullins–Sekerka); hydrodynamic instabilities (Rayleigh–Bénard, shear/viscous, Marangoni); morphoelastic buckling and wrinkling; Saffman–Taylor fingering.

- 21\. Turing Patterns, Swift-Hohenberg Equation, Amplitude Equation
- 22\. From Perturbation to Amplitude Equation, NWS Equation, Multiple Scale Analysis
- 23\. Newell-Whitehead-Segel Equation, Secondary Instabilities
- 24\. Eckhaus & Zigzag Instabilities, Phase Winding Solution
- 25\. The Complex Ginzburg-Landau Equation
- 26\. Stability of the Complex Ginzburg-Landau Equation
- 27\. Phase Diagram of Complex Ginzburg-Landau Equation
- 28\. Application of CGLE, May-Leonard Model, Evolutionary Games
- 29\. Reaction-Diffusion Waves, Front Spreading Velocity
- 30\. The Kuramoto Model & Synchronization Transition

### Geometry & confinement
Curvature and weakly non-flat geometry; boundary conditions; selection by confinement and heterogeneity.

- 31\. Mass-Conserving Reaction-Diffusion Systems, E. coli Min Oscillations
- 32\. Mass-Redistribution Instability, Adiabatic Elimination
- 33\. Geometric Construction of Stationary Patterns
- 34\. Steady-State Phase Diagram & Bulk Gradient Effects
- 35\. Bulk-Boundary Coupling & Nucleotide Exchange
- 36\. Stationary Solutions & Emergent Saturation Effect
- 37\. Linear Stability Analysis in Box Geometry
- 38\. Modeling the E. coli Min System: The Skeleton Model
- 39\. Linear Stability Analysis for Min Skeleton Model
- 40\. The Min Switch Model & Robustness

### Coarsening dynamics
Advanced topics on coarsening, phase field models, and interfacial dynamics.

- 41\. Model A/B/C, Solidification Dynamics & Phase Field Model
- 42\. Coarsening Dynamics in Mass-Conserving Systems
- 43\. Universal Growth Laws & Logarithmic Coarsening
- 44\. Scaling Analysis of the Coarsening Law
- 45\. From Peaks to Mesas: Interaction via Exponential Tails
- 46\. Arrested Coarsening via Weakly Broken Mass Conservation
- 47\. Linear Stability Analysis of Weakly Broken Mass Conservation
- 48\. The Mechanism of Mesa Splitting: Source-Driven Lateral Instability
- 49\. Interrupted Coarsening of Mesas & The Gibbs-Thomson Relation
- 50\. Surface Tension in Multicomponent Reaction-Diffusion Systems
- 51\. The Gibbs-Thomson Relation in Reaction-Diffusion Systems


## Usage

Each Python file corresponds to specific topics covered in the lecture series. The code serves as practical implementations of the theoretical concepts presented in the YouTube videos, developed as part of self-study and learning notes.

Here are some code output demonstrations:

**[Lecture 2: Dynamical Systems - Phase Transitions to Evolutionary Game Theory](en/2.%20Dynamical%20Systems%20-%20From%20Physical%20Phase%20Transitions%20to%20Evolutionary%20Game%20Theory.md)**

| Landau Free Energy & Phase Transition | Evolutionary Game Dynamics | Replicator Equation Flow |
|:---:|:---:|:---:|
| ![Landau](cn/assets/images/02_003_621bc300-3817-4cf6-bb3a-baaced59b8f5.png) | ![Game](cn/assets/images/02_004_89cf16e7-61e4-4535-836a-4314cd5c8b69.png) | ![Replicator](cn/assets/images/02_007_dc98adf0-a196-4959-9fec-b771eefa8a62.png) |

**[Lecture 3: Fixed Points, Saddle-Node & Pitchfork Bifurcations](en/3.%20Dynamical%20Systems%20-%20Fixed%20Points%2C%20Saddle-Node%20and%20Pitchfork%20Bifurcations.md)**

| Bifurcation Diagram | Phase Portrait |
|:---:|:---:|
| ![Bifurcation](cn/assets/images/03_006_d0be1090-a825-48c8-af21-2cf9cfa4f7cc.png) | ![Phase](cn/assets/images/03_009_70061829-22ca-4bb7-933a-2f611b073267.png) |

**[Lecture 4: Cusp Bifurcation, Infection Model & Transcritical Bifurcation](en/4.%20Dynamical%20Systems%20-%20Cusp%20Bifurcation%2C%20Infection%20Model%20and%20Transcritical%20Bifurcation.md)**

| Cusp Bifurcation Surface | Infection Dynamics |
|:---:|:---:|
| ![Cusp](cn/assets/images/04_008_2b441799-2822-4439-bcaa-67292e672b48.png) | ![Infection](cn/assets/images/04_011_6e149644-d36a-456c-97a5-43c30d64b8c3.png) |

**[Lecture 5: Two-Dimensional Mass-Conserving Systems](en/5.%20Dynamical%20Systems%20-%20Two-Dimensional%20Mass-Conserving%20Systems%20and%20the%20Jacobian%20Matrix.md)**

| Jacobian Analysis & Phase Space |
|:---:|
| ![Jacobian](cn/assets/images/05_010_8baf3743-4296-492f-b3df-a015f9b4407c.png) |

**[Lecture 6: Two-Dimensional Non-Conserving Systems & Oscillators](en/6.%20Dynamical%20Systems%20-%20Two-Dimensional%20Non-Conserving%20Systems%20and%20Oscillators.md)**

| Hopf Bifurcation | Limit Cycle | Van der Pol Oscillator |
|:---:|:---:|:---:|
| ![Hopf](cn/assets/images/06_004_e376e171-e6ec-4ed1-8659-5b8579388051.png) | ![Limit](cn/assets/images/06_007_0e3242d5-a770-411b-bffd-c151adb26aad.png) | ![VdP](cn/assets/images/06_009_6a812e94-584f-4e92-90ee-70d202be59dc.png) |

**[Lecture 7: Poincaré-Bendixson Theorem, RhoGTPase & Rock-Paper-Scissors](en/7.%20Dynamical%20Systems%20-%20Poincar%C3%A9-Bendixson%20Theorem%2C%20RhoGTPase%20Oscillator%20and%20Rock-Paper-Scissors.md)**

| RhoGTPase Oscillator | Rock-Paper-Scissors Dynamics | Heteroclinic Orbit |
|:---:|:---:|:---:|
| ![Rho](cn/assets/images/07_004_059881e5-e0fc-4efb-90b8-f1aa7cf5d15b.png) | ![RPS](cn/assets/images/07_006_aea062d8-475b-4948-82f3-0ba967f1ce27.png) | ![Hetero](cn/assets/images/07_011_6315e0dc-d8d7-4c2b-ae2e-cab04a874edf.png) |

**[Lecture 8: Replicator Dynamics & Lyapunov Functions](en/8.%20Dynamical%20Systems%20-%20Replicator%20Dynamics%20and%20Spatial%20Extension%20-%20Lyapunov%20Functions%2C%20Boson%20Condensation%20and%20the%20Diffusion%20Equation.md)**

| Lyapunov Function | Diffusion Equation |
|:---:|:---:|
| ![Lyapunov](cn/assets/images/08_008_6fc576fa-0149-454b-93f4-9dfe7af2de8f.png) | ![Diffusion](cn/assets/images/08_009_9eeceb55-ff8d-4bff-91ca-c333813ae8aa.png) |

**[Lecture 9: Ginzburg-Landau Theory & Allen-Cahn Equation (Model A)](en/9.%20Relaxational%20Dynamics%20-%20Ginzburg-Landau%20Theory%20and%20the%20Allen-Cahn%20Equation.md)**

<video controls width="100%">
  <source src="cn/assets/images/model_a_dynamics.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

*Model A (Allen-Cahn) dynamics: Spinodal decomposition, domain formation, and coarsening*

**[Lecture 10: Interface Dynamics & Curvature-Driven Flow](en/10.%20Relaxational%20Dynamics%20-%20Interface%20Dynamics%20and%20Curvature-Driven%20Flow.md)**

<video controls width="100%">
  <source src="cn/assets/images/droplet_shrinkage.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

*Curvature-driven shrinkage of a circular droplet following Allen-Cahn dynamics*

**[Lecture 11: Thermodynamics & Phase Separation in Liquid Mixtures](en/11.%20Thermodynamics%20and%20Phase%20Separation%20in%20Liquid%20Mixtures.md)**

| Binary Phase Diagram |
|:---:|
| ![Phase Diagram](cn/assets/images/11_011_5fb64d92-89ff-4a37-9344-b3f4286a500c.jpg) |

**[Lecture 12: Cahn-Hilliard Equation & Osmotic Pressure](en/12.%20Liquid%20Mixture%20Dynamics%2C%20Osmotic%20Pressure%20and%20the%20Cahn-Hilliard%20Equation.md)**

<video controls width="100%">
  <source src="cn/assets/images/CahnHilliard_vs_McRD.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

*Comparison: Cahn-Hilliard (Model B) vs Mass-Conserving Reaction-Diffusion dynamics*

**[Lecture 13: Lattice Gas Model & Free Energy](en/13.%20Lattice%20Gas%20Model%20and%20Free%20Energy.md)**

![Lattice Gas Simulation](cn/assets/images/lattice_gas_vs_free_energy.gif)

*Monte Carlo simulation of lattice gas model showing phase separation and free energy evolution*

**[Lecture 14: Phase Separation Dynamics & Dispersion Relation](en/14.%20Phase%20Separation%20Dynamics.md)**

| Crowd Oscillation Dynamics | Spectral Heatmap |
|:---:|:---:|
| ![Crowd](cn/assets/images/crowd_oscillation.gif) | ![Spectrum](cn/assets/images/14_008_d76b6a87-990d-46a0-9c1b-f0830b67a6fb.png) |

**[Lecture 15: Nucleation Theory & Interface Dynamics](en/15.%20Nucleation%20Theory%20and%20Interface%20Dynamics.md)**

| 3D Cahn-Hilliard Droplets | Evolution Slices |
|:---:|:---:|
| ![3D Droplet](cn/assets/images/ch3d_droplet_rotate.gif) | ![Slices](cn/assets/images/ch3d_slices.gif) |

| Concentration Profile | Nucleation Barrier |
|:---:|:---:|
| ![Profile](cn/assets/images/15_008_24e1273d-2e11-44fa-bcd3-cc2fa30e08fe.png) | ![Barrier](cn/assets/images/15_009_9639068f-e1d9-4bf2-80d5-5f5abf4f50fa.png) |


## Prerequisites

- Statistical mechanics and thermodynamics
- Dynamical systems theory
- Differential equations
- Basic knowledge of field theory (helpful but not required)


## Citation

If you find these notes or code useful, you may cite this repository using the following BibTeX entry:

```bibtex
@misc{liu2025selforganization,
  author       = {Liu, Zhihang},
  title        = {Self-Organization and Pattern Formation: Course Notes and Code},
  year         = {2025},
  howpublished = {\url{https://zhihangliu.cn/Self-Organization-and-Pattern-Formation/}},
  note         = {Course by Prof. Erwin Frey, LMU Munich, Winter 2025/2026}
}
```
