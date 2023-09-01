---
title: Bachelor Thesis
date: 2023-08-30 11:33:00 +0800
categories: [Physics, AI]
tags: [Physics]
pin: true
math: true
mermaid: true
image:
    path: https://www.alcf.anl.gov/sites/default/files/styles/965x543/public/2019-09/Habib_hacc_outer_rim_1600x800.png?itok=2b5C0vZg
    alt: Image of lambda-CDM cosmology simulation of the early univerve
---
# Machine Learning Quintessence Dark Energy Potentials via Symbolic Regression and Bayesian Inference

Quintessence is a hypothesized scalar field forming a popular model for dark energy, the nonluminous energy density in the Universe responsible for the accelerated expansion of space.
This field can change in time and henceforth solve difficulties arising in the ΛCDM standard
model. The dynamics of the field are governed by its self-interaction potential. To infer the
form of this potential by observational data, we construct a suitable potential using machine
learning and derive the observable quantities then numerically. Symbolic regression is a form
of supervised learning, with the aim of finding the most likely symbolic representations of
mathematical models for a given data set.

We demonstrate how to build a symbolic regression machine learning pipeline, searching
for self-interaction potentials of the quintessence scalar field Lagrangian. To construct this
pipeline, we first implement a symbolic regression algorithm searching for classical Lagrangian
potentials and discuss its ability to solve the problem for artificial data. For upgrading this algorithm to search for quintessence potentials, we assembly an automatic process of evaluating
the individual potentials and comparing the resulting models to data of supernovae type Ia
with Bayesian statistics. We discuss the difficulties of automatically assigning Bayesian measures due to the high complexity of the fitness process as well as limited computing strength.
We conclude, that our symbolic regression pipeline is capable of finding suitable quintessence
potentials under consideration of sufficiently diverse function populations.

Read more [here](/assets/img/media/BA_Machine_learning_Quintessence_Dark_Energy_Potentials.pdf)