# stochastic_finite_model
R software for fitting and simulating the stochastic finite (SF) model (Pridham and Rutenberg 2023).

# **Getting started**
Download the latest release (tags, right hand side of screen). Installation is unnecessary, you can open sf.R in R and start fitting data now. I recommend you read the vignettes which are written in R notebook (.Rmd) using RStudio. (You can read the .html instead but it is less convenient.)

You should start with sf_vignette_fitsf.Rmd. You will need to source sf.R. You'll also need to enter the location of this repository's contents, which is stored in `outputDir`.

# **What it does**
The SF model uses longitudinal data to estimate an interaction network and equilibrium (steady-state) behaviour. It is a linear approximation of a more generate stochastic (Wiener) process model.

# **Requirements**
Written for R version 4.1.1. RStudio is recommended for the vignettes. Tested on Windows, MacOS and Ubuntu.

Required packages: MASS and survival. 

Recommended: parallel

Optional: nloptr.

# **How to use**
See vignettes for functionality. Vignettes are written as Rstudio notebooks (version 2023-03-00) (.Rmd). Vignettes run quickly and each block of code should run in < 10 minutes on a typical personal computer.

The .html are included for completeness but have not been optimized (so they may be ugly or lack figures).

# **Key file**
sf.R is main file which contains the functions needed to fit the model.

# **Tips**
1. Although not the default, it is worth trying the linear regression estimator: `FitFun=FitLM.SF` e.g. `sf = FitSF(y=y,x=x,dt=dt,FitFun=FitLM.SF)`. The linear regression was proven in a subsequent paper, and usually works a lot better.
2. If you have reduced rank W the model will likely fail to fit, this is due to collinearity problems. You can (a) reduce the parameters used in the fit or (b) pre-process using principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y)-2)` will drop the last 2 principal components, thus reducing collinearity. You can also try fitting a diagonal model, this works pretty well if you preprocess with principal component analysis e.g. `sf = FitSF(y=y,x=x,dt=dt,Npc=ncol(y),diagonalW=T)`.

# **Kidney directory**
See [Pridham2024-ql] for details. The directory contains all necessary files for that publication together with a vignette on how to compute and use the natural variable by using the coefficient files.

# **How does the math work?**
See Pridham and Rutenberg (2023). The supplemental is quite detailed.

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

An application that helps explain how the model works:
Pridham, G., Tennankore, K. K., Rockwood, K., Worthen, G. & Rutenberg, A. D. Systems-level health of patients living with end-stage kidney disease using standard lab values. arXiv [q-bio.QM] (2024)

@ARTICLE{Pridham2024-ql,
  title         = "Systems-level health of patients living with end-stage
                   kidney disease using standard lab values",
  author        = "Pridham, Glen and Tennankore, Karthik K and Rockwood,
                   Kenneth and Worthen, George and Rutenberg, Andrew D",
  month         =  may,
  year          =  2024,
  url           = "http://arxiv.org/abs/2405.20523",
  archivePrefix = "arXiv",
  eprint        = "2405.20523",
  primaryClass  = "q-bio.QM",
  arxivid       = "2405.20523"
}

# **Common errors**
-Individuals are sorted by their age, so make sure age captures changes in time e.g. if you measure someone at age 20 then 3 months later then the followup age should be 20.25 (NOT 20!)
-If you use Reshape() with an id column that can't be coerced into numeric, it will crash. I'll remove this in the future.
-If you use Reshape() with character sex it will fail to perform sex-specific normalization. In this case you should binarize sex or, (i) set sexSpecific=FALSE or (ii) set preprocess=FALSE.
-namespace issues might also crop up (for me it's Impute), you can always add a post hoc wrapper for this or just load the package last.
