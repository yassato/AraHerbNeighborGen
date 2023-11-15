# export coefficients from LASSO results

import glmnet_python
from glmnet import glmnet

# Import relevant modules and setup for calling glmnet
import scipy, importlib, pprint, matplotlib.pyplot as plt, warnings
from glmnet import glmnet; from glmnetPlot import glmnetPlot
from glmnetPrint import glmnetPrint; from glmnetCoef import glmnetCoef; from glmnetPredict import glmnetPredict
from cvglmnet import cvglmnet; from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot; from cvglmnetPredict import cvglmnetPredict

import pandas
import numpy as np
import gzip

import pickle
def save_object(obj, filename):
    with gzip.open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


# load glmnet results

fpath = "./output/"
fname = "ScoreS3JPN_glmnetLassoMAF5"

fn = fpath + fname + ".pkl.gz"
with gzip.open(fn, "rb") as f: 
        fit = pickle.load(f)

out = pandas.DataFrame(fit["lambdau"],fit["df"])
df_fn = fpath + fname + "_DFlambda.csv"
out.to_csv(df_fn)

beta = pandas.DataFrame(fit["beta"])
beta_fn = fpath + fname + "_coefALL.csv"
beta.to_csv(beta_fn)

