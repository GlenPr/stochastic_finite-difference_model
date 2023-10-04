# stochastic_finite_model
R software for fitting and simulating the stochastic finite (SF) model (Pridham and Rutenberg 2023).

# **Getting started**
Download the latest release (tags, right hand side of screen). You can open sf.R in R and start fitting data now. I recommend you read the vignettes which are written in R notebook (.Rmd) using RStudio. (You can read the .html instead but it is less convenient.)

You should start with sf_vignette_fitsf.Rmd. You will need to source sf.R. You'll also need to enter the location of this repository's contents, which is stored in `outputDir`.

# **What it does**
The SF model uses longitudinal data to estimate an interaction network and equilibrium (steady-state) behaviour. It is a linear approximation of a more generate stochastic (Wiener) process model.

# **Requirements**
Written for R version 4.1.1. RStudio is recommended for the vignettes.

Required packages: MASS and survival. 

Recommended: parallel

Optional: nloptr.

# **How to use**
See vignettes for functionality. Vignettes are written as Rstudio notebooks (version 2023-03-00) (.Rmd). 

The .html are included for completeness but have not been optimized (so they may be ugly or lack figures).

# **Key file**
sf.R is main file which contains the functions needed to fit the model.

# **Tips**
1. Although not the default, it is worth trying the linear regression estimator: `FitFun=FitLM.SF` e.g. `sf = FitSF(y=y,x=x,dt=dt,FitFun=FitLM.SF)`. The linear regression estimator has yet to be formally validated but likely has lower bias than the default estimator in most cases.
2. If you have reduced rank W the model will likely fail to fit, this is due to collinearity problems. You can (a) reduce the parameters used in the fit or (b) pre-process using principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y)-2)` will drop the last 2 principal components, thus reducing collinearity. You can also try fitting a diagonal model, this works pretty well if you preprocess with principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y),diagonalW=T)`.

# **How does the math work?**
See Pridham and Rutenberg (2023).

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

# **See also**
For more information, see:
Pridham, G. & Rutenberg, A. D. Dynamical network stability analysis of multiple biological ages provides a framework for understanding the aging process. arXiv [q-bio.QM] Preprint at http://arxiv.org/abs/2309.10005 (2023)
