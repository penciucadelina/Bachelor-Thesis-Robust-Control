# Robust State Feedback Linearization: A Data-Driven Approach

![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)
![Control Systems](https://img.shields.io/badge/Control_Systems-Data_Driven-success.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

## 📌 Overview
This repository contains the complete MATLAB/Simulink implementation for my **BSc Thesis (Grade: 10/10, Valedictorian)** in Automation and Computer Science at the Technical University of Cluj-Napoca.

The project implements a **data-driven control algorithm** that identifies state feedback linearization functions (coordinate transformations, decoupling matrices) directly from input-state data, bypassing the need for explicit analytical models. It specifically addresses the challenges of **noisy experimental data** by employing advanced numerical techniques and signal processing.

📄 **Full Thesis Document:** [Read the PDF here](Licenta_Penciuc_Adelina_Final.pdf) *(Asigură-te că link-ul corespunde cu numele real al PDF-ului tău)*

## 🚀 Key Features & Algorithmic Implementation
*   **Data-Driven Identification:** Reconstructs the required nonlinear control functions using predefined dictionaries and algebraic data matrices.
*   **Noise Handling & Signal Processing:** 
    *   Implements **Savitzky-Golay filtering** for state smoothing.
    *   Replaces highly sensitive derivative-based formulations with an **Integral Formulation** to mitigate noise amplification.
*   **Numerical Optimization:** 
    *   Uses **Singular Value Decomposition (SVD)** to extract null-space solutions for ideal datasets.
    *   Applies **Tikhonov Regularization** to stabilize the numerical identification process for noisy data.
*   **Closed-Loop Validation:** Integrates the identified nonlinear functions with an **LQR controller** to ensure asymptotic stability in closed-loop simulations.
*   **Custom GUI:** Includes a bespoke MATLAB App Designer interface for automated data processing, metric visualization (RMSE, residual validation), and closed-loop performance analysis.

## 🛠️ Technologies & Tools
*   **Core:** MATLAB, Simulink (Modeling & Simulation)
*   **Algorithms:** SVD, Tikhonov Regularization, Least Squares Optimization, Numerical Integration.
*   **Control Theory:** State-Space Modeling, LQR, Feedback Linearization (SISO & MIMO systems).

## 📈 Results
The data-driven integral approach with Tikhonov regularization successfully stabilized both SISO and MIMO nonlinear systems under noisy conditions ($SNR = 25 dB$), keeping the control effort (RMS command) within feasible physical limits and achieving near-zero steady-state errors.
