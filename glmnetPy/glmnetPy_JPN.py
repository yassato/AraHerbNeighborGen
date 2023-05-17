import gc

import glmnet_python
from glmnet import glmnet

# Import relevant modules and setup for calling glmnet
import scipy, importlib, pprint, matplotlib.pyplot as plt, warnings
from glmnet import glmnet; from glmnetPlot import glmnetPlot
from glmnetPrint import glmnetPrint; from glmnetCoef import glmnetCoef; from glmnetPredict import glmnetPredict
from cvglmnet import cvglmnet; from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot; from cvglmnetPredict import cvglmnetPredict

import pandas as pd
import numpy as np
import gzip

import pickle
def save_object(obj, filename):
    with gzip.open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


#load data
with open("../data/glmnetPyJPNselfMAF5.pkl", 'rb') as f:
    self = pickle.load(f)

with open("../data/glmnetPyJPNneiS1MAF5.pkl", 'rb') as f:
    nei = pickle.load(f)

with open("../data/glmnetPyJPNphenoMAF5map.pkl", 'rb') as f:
    pheno = pickle.load(f)

with open("../data/glmnetPyJPNcovMAF5map.pkl", 'rb') as f:
    cov = pickle.load(f)


pheno = pd.DataFrame(pheno.values(), index=pheno.keys()).T

X = np.concatenate([np.transpose(cov), self, nei])
X = X.astype(scipy.float64)
X = np.transpose(X)
cov.shape[1]

#load pheno
y = scipy.array(pheno.sucker, dtype = scipy.float64)
y = np.log(y+1)
y = (y - np.mean(y))/np.std(y)

# set penalty.factors
wts = scipy.ones([1,X.shape[1]])
wts[0,range(cov.shape[1])] = 0

del self
del nei
gc.collect()

# lasso
fit = glmnet(x = X, y = y, family = "gaussian", alpha = 1.0, nlambda = 100, penalty_factor = wts, standardize=True)
save_object(fit, "../output/suckerS1JPN_glmnetLassoMAF5.pkl.gz")
del fit

# change y and run lasso again
y = scipy.array(pheno.chewer, dtype = scipy.float64)
y = np.log(y+1)
y = (y - np.mean(y))/np.std(y)

fit = glmnet(x = X, y = y, family = "gaussian", alpha = 1.0, nlambda = 100, penalty_factor = wts, standardize=True)
save_object(fit, "../output/chewerS1JPN_glmnetLassoMAF5.pkl.gz")
del fit

# change y and run lasso again
y = scipy.array(pheno.richness, dtype = scipy.float64)
y = np.log(y+1)
y = (y - np.mean(y))/np.std(y)

fit = glmnet(x = X, y = y, family = "gaussian", alpha = 1.0, nlambda = 100, penalty_factor = wts, standardize=True)
save_object(fit, "../output/richnessS1JPN_glmnetLassoMAF5.pkl.gz")
del fit

# change y and run lasso again
y = scipy.array(pheno.PxPr, dtype = scipy.float64)
y = np.log(y+1)
y = (y - np.mean(y))/np.std(y)

fit = glmnet(x = X, y = y, family = "gaussian", alpha = 1.0, nlambda = 100, penalty_factor = wts, standardize=True)
save_object(fit, "../output/PxPrS1JPN_glmnetLassoMAF5.pkl.gz")
del fit

# change y and run lasso again
y = scipy.array(pheno.Score, dtype = scipy.float64)
# y = np.log(y+1)
# y = (y - np.mean(y))/np.std(y)

fit = glmnet(x = X, y = y, family = "gaussian", alpha = 1.0, nlambda = 100, penalty_factor = wts, standardize=True)
save_object(fit, "../output/ScoreS1JPN_glmnetLassoMAF5.pkl.gz")
del fit

