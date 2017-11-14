# N15 LIM

This project integrates N15 values into Linear Inverse Models for ecosystem modeling. Below are the actual scripts, available in multiple formats. A major goal of this project is to facilitate the use of this tool within the oceanographic community. We welcome any questions or concerns, especially regarding how to get started with Linear Inverse Modeling.

1. [N15 Setup File](https://github.com/tbrycekelly/N15-LIM/blob/master/SetMatricesN15RW.ipynb): Reads in the spreadsheets and prepares the matricies and vectors used in the LIM. Also loads results from the Forward model for use in the LIM.
2. [Model Initialization Script](https://github.com/tbrycekelly/N15-LIM/blob/master/RunN15InverseRW.ipynb): Gets the run-time enviroment ready and handles all the bunr-in stages as well as the actual model runs. Saves all results at the end.
3. [A Modified MCMC Sampling Algorithm](https://github.com/tbrycekelly/N15-LIM/blob/master/xsampleN15.r): Is called by the __Model Initialization Script__, not meant to be edited from one model to the next. Adapted from Van den Meersche's _xsample_ script.
4. [With these additional functions](https://github.com/tbrycekelly/N15-LIM/blob/master/ExternalFunctions.ipynb): Includes the function that updates the N15-related equations during the random walk. Also an ideal file for accessory functions.


---
### Multiple Formats

In general, the user accessible scripts are available as with jupyter notebook files or as straight R code.

* N15 Setup File: [Notebook](https://github.com/tbrycekelly/N15-LIM/blob/master/SetMatricesN15RW.ipynb), [R code](https://github.com/tbrycekelly/N15-LIM/blob/master/SetMatricesN15RW.r), [Matlab](https://github.com/tbrycekelly/N15-LIM/blob/master/SetMatricesN15RW.m)
* Model Initialization Script: [Notebook](https://github.com/tbrycekelly/N15-LIM/blob/master/RunN15InverseRW.ipynb), [R code](https://github.com/tbrycekelly/N15-LIM/blob/master/RunN15InverseRW.R)
* A Modified MCMC Sampling Algorithm: [R code](https://github.com/tbrycekelly/N15-LIM/blob/master/xsampleN15.r)
* With these additional functions: [Notebook](https://github.com/tbrycekelly/N15-LIM/blob/master/ExternalFunctions.ipynb), [R code](https://github.com/tbrycekelly/N15-LIM/blob/master/ExternalFunctions.R)


### About Notebook Files

Since not everyone is familiar with the Jupyter Notebook platform, let us say a few words about it. Jupyter, formally iPython Notebook, is an interactive scripting environment for use in data analysis, script development and, importanly, for sharing scripts and workflows with others. While designed with Python in mind, the platform is nearly language agnostic with ready backends (i.e. kernels) for numerous languages: R, MatLab, C, Java, Clojure, Lisp, Fortran, etc. Do yourself a favor and check it out, it is certainly worth the time to learn about this valuable, open source resource.

---

### Citation

Stukel M.R., Decima M., Kelly T.B. A new approach for incorporating 15N isotopic data into linear inverse ecosystem models with Markov Chain Monte Carlo sampling. 