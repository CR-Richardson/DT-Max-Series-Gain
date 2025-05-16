# DT-Max-Series-Gain
This code acompanies the paper [Strengthened stability analysis of discrete-time Lurie systems involving ReLU neural networks](https://proceedings.mlr.press/v242/richardson24a.html) where the experimental setup is detailed in the *numerical examples* section. The code computes the maximum series gain for which global asymptotic stability is verified using various criteria. Furthermore, the number of decision variables used by the criteria is also returned to compare the complexity. This is tested on a number of example discrete-time Lurie systems assumed to have repeated ReLU nonlinearities.  

### Authors:
* Carl R Richardson (cr2g16@soton.ac.uk)
* Matthew C Turner (m.c.turner@soton.ac.uk)

## Prerequisites
All the code is written in MATLAB. The LMI's are solved using the *Robust Control Toolbox* which must be installed as an add-on.

## Overview
The repository is organised as follows:
- `DT_Max_Series_Gain.m` The master script. It loops through each example, computing the maximum series gain (and # of decision variables) according to each criterion,  and displays the results.
- `DT_Examples.m` Defines the (A,B,C,D,Ts) parameters of the example systems.
- `DT_Circle.m` Implementation of the DT Circle Criterion - See Section 5 (Haddad & Bernstein, 1994).
- `DT_Circle_Like.m` Implementation of the DT Circle-like Criterion - See Theorem 12.
- `DT_Popov.m` Implementation of the DT Popov Criterion - See Section 6 (Haddad & Bernstein, 1994).
- `DT_Popov_Like.m` Implementation of the Relaxed (H=I) DT Popov-like Criterion - See Remark 15.
- `DT_Park.m` Implementation of the DT Park Criterion - See Lemma 2 and Theorem 1 (Park, 2019).
- `DT_Aizerman.m` Computes the Nyquist gain (based on the Aizerman Conjecture) for each discrete-time example.
- `Poster.pdf` Poster associated with the paper.

## Getting Started
- Add all files to the path.
- Run `DT_Max_Series_Gain.m` to repeat the experiments in the paper or select a subset of the examples by defining them in the *Ex_array* variable.
- Run `DT_Aizerman.m` to repeat the Nyquist gain calculations.

Note: `DT_Park.m` takes much longer to converge on a final value of alpha compared the other criteria.
