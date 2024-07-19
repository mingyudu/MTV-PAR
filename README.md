# MTV-PAR

This is the repository for the paper 

>Time-varying $\ell_0$ optimization for spike inference from multi-trial calcium recordings'

It contains the code for the simulation study and real data analysis of the Multi-Trial time-Varying Penalized Auto-Regression (MTV-PAR) method.

## Introduction

Optical imaging of genetically encoded calcium indicators is a powerful tool to record the activity of a large number of neurons simultaneously over a long period of time from freely behaving animals. 
However, determining the exact time at which a neuron spikes and estimating the underlying firing rate from calcium fluorescence data remains challenging, especially for calcium imaging data obtained from a longitudinal study. 
We propose a Multi-Trial time-Varying $\ell_0$ Penalized Auto-Regression (MTV-PAR) method to jointly detect spikes and estimate firing rates by robustly integrating evolving neural dynamics across trials.
Our simulation study shows that the proposed method performs well in both spike detection and firing rate estimation. 
We demonstrate the usefulness of our method on calcium fluorescence trace data from two studies, with the first study showing differential firing rate functions between two behaviors and the second study showing evolving firing rate functions across trials due to learning. 

## Cite this work

If you find any of the source code in this repository useful for your work, please cite:

> Time-varying $\ell_0$ optimization for spike inference from multi-trial calcium recordings.
