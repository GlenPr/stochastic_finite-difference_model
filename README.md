# stochastic_finite_model
R software for fitting and simulating the stochastic finite (SF) model (Pridham 2023).

# **How do I use it?**
Download the folder. You can open sf.R in R and start fitting data now. I recommend you read the vignettes which are written in R notebook (.Rmd) using RStudio.

# **What it does**
The SF model uses longitudinal data to estimate an interaction network and equilibrium (steady-state) behaviour. It is a linear approximation of a more generate stochastic (Wiener) process model.

# **Requirements**
Written for R version 4.1.1. Required packages: MASS and survival. 

# **How to use**
See vignettes for functionality. Vignettes are written as Rstudio notebooks (version 2023-03-00) (.Rmd). 

The .html are included for completeness but have not been optimized (so they may be ugly or lack figures).

# **Key file**
sf.R is main file which contains the functions needed to fit the model.

# **How the math work?**
See the supplemental of Pridham and Rutenberg (2023).

# **Cite as**
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
