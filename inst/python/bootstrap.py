import numpy as np
import pandas as pd
from BOAmodel import *
import random
from collections import defaultdict

def bootstrapBOA(supp, maxlen, N, alpha_1, beta_1, alpha_2, beta_2, prior_type, alpha_l, beta_l, lmd, nu, Niteration, Nchain, trainProp, print, df, Y, bootstrap=True, reps=1, sampleSize=-1):
    Y = np.array(Y)
    lenY = len(Y)
    # Store results
    allRules = dict()
    allStats = pd.DataFrame(columns = ['accuracy', 'tpr', 'fpr'], index = range(0,reps))
    #allTestData = dict()
    allIndices = dict()

    ## Floating point error
    old_settings = np.seterr(all='ignore')  #seterr to known value
    np.seterr(over='ignore')

    for i in range(0, reps):
        if bootstrap:
            sample_index = np.random.choice(lenY, sampleSize) # Bootstrap
            bootstrap_df = df.iloc[sample_index]
            bootstrap_Y = Y[sample_index]
            allIndices[i] = sample_index
        else:
            bootstrap_df = df
            bootstrap_Y = Y
            sampleSize=lenY
        train_index = random.sample(range(sampleSize),int(trainProp*sampleSize))
        test_index = [i for i in range(sampleSize) if i not in train_index]
        #allTestData[i] = test_index
        #allTestData[i] = [bootstrap_df.iloc[test_index], bootstrap_Y[test_index]]
        BOAmodel = BOA(bootstrap_df.iloc[train_index],bootstrap_Y[train_index], print)
        BOAmodel.generate_rules(supp,maxlen,N)
        BOAmodel.set_parameters(alpha_1,beta_1,alpha_2,beta_2,prior_type,alpha_l,beta_l,lmd,nu)
        rules = BOAmodel.SA_patternbased(Niteration,Nchain)
        allRules[i] = rules

        # Return stats on all data:
        #Yhat = predict(rules,bootstrap_df)
        #TP,FP,TN,FN = getConfusion(Yhat,bootstrap_Y)
        # Return results for test data:
        Yhat = predict(rules,bootstrap_df.iloc[test_index])
        TP,FP,TN,FN = getConfusion(Yhat,bootstrap_Y[test_index])

        accuracy = float(TP+TN)/(TP+TN+FP+FN)
        tpr = float(TP)/(TP+FN)
        fpr = float(FP)/(FP+TN)
        allStats.loc[i] = [accuracy, tpr, fpr]
    return allRules, allIndices, allStats


#def BOA(A):
