# stochastic_finite_model
R software for fitting and simulating the stochastic finite (SF) model (Pridham 2023).

The SF model uses longitudinal data to estimate an interaction network and equilibrium (steady-state) behaviour. It is a linear approximation of a more generate stochastic (Wiener) process model.

Written for R version 4.1.1. Required packages: MASS and survival. See vignettes for functionality.

sf.R contains the functions needed to fit the model.

Please citation
Pridham, G. & Rutenberg, A. D. Network dynamical stability analysis of homeostasis reveals ‘mallostasis’: biological equilibria drifting towards worsening health with age. arXiv [q-bio.OT] (2023)

Bibtex
@ARTICLE{Pridham2023-yt,
  title         = "Network dynamical stability analysis of homeostasis reveals
                   ``mallostasis'': biological equilibria drifting towards
                   worsening health with age",
  author        = "Pridham, Glen and Rutenberg, Andrew D",
  month         =  jun,
  year          =  2023,
  url           = "http://arxiv.org/abs/2306.12448",
  archivePrefix = "arXiv",
  eprint        = "2306.12448",
  primaryClass  = "q-bio.OT",
  arxivid       = "2306.12448"
}
