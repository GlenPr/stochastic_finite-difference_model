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
1. Although not the default, it is worth trying the linear regression estimator: `FitFun=FitLM.SF` e.g. `sf = FitSF(y=y,x=x,dt=dt,FitFun=FitLM.SF)`. The linear regression was proven in a subsequent paper, and usually works a lot better.
2. If you have reduced rank W the model will likely fail to fit, this is due to collinearity problems. You can (a) reduce the parameters used in the fit or (b) pre-process using principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y)-2)` will drop the last 2 principal components, thus reducing collinearity. You can also try fitting a diagonal model, this works pretty well if you preprocess with principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y),diagonalW=T)`.

# **How does the math work?**
See Pridham and Rutenberg (2023).

# **Cite as**
Pridham, G. & Rutenberg, A. D. Network dynamical stability analysis reveals key ‘mallostatic’ natural variables that erode homeostasis and drive age-related decline of health. Sci. Rep. 13, 1–12 (2023). https://doi.org/10.1038/s41598-023-49129-7

Bibtex

@ARTICLE{Pridham2023-pi,
  title     = "Network dynamical stability analysis reveals key ``mallostatic''
               natural variables that erode homeostasis and drive age-related
               decline of health",
  author    = "Pridham, Glen and Rutenberg, Andrew D",
  journal   = "Sci. Rep.",
  publisher = "Nature Publishing Group",
  volume    =  13,
  number    =  1,
  pages     = "1--12",
  month     =  dec,
  year      =  2023,
  url       = "https://www.nature.com/articles/s41598-023-49129-7",
  language  = "en",
  issn      = "2045-2322, 2045-2322",
  doi       = "10.1038/s41598-023-49129-7"
}


# **See also**
It is wise to use the least-squares estimator which is proven in this paper:

Pridham, G. & Rutenberg, A. D. Dynamical network stability analysis of multiple biological ages provides a framework for understanding the aging process. J. Gerontol. A Biol. Sci. Med. Sci. (2024) doi:10.1093/gerona/glae021

Bibtex

@ARTICLE{Pridham2024-cy,
  title    = "Dynamical network stability analysis of multiple biological ages
              provides a framework for understanding the aging process",
  author   = "Pridham, Glen and Rutenberg, Andrew D",
  journal  = "J. Gerontol. A Biol. Sci. Med. Sci.",
  month    =  jan,
  year     =  2024,
  url      = "http://dx.doi.org/10.1093/gerona/glae021",
  keywords = "biological age; complexity; eigen analysis; systems biology",
  language = "en",
  issn     = "1079-5006, 1758-535X",
  pmid     = "38206765",
  doi      = "10.1093/gerona/glae021"
}
