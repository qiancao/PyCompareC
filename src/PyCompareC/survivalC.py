import numpy as np

def csign(X_i, Censor_i, X_j, Censor_j):
    if X_i > X_j:
        if Censor_i == 1 and Censor_j == 1:
            return 1
        if Censor_i == 1 and Censor_j == 0:
            return 0
        if Censor_i == 0 and Censor_j == 1:
            return 1
        if Censor_i == 0 and Censor_j == 0:
            return 0
    
    if X_i < X_j:
        if Censor_i == 1 and Censor_j == 1:
            return -1
        if Censor_i == 1 and Censor_j == 0:
            return -1
        if Censor_i == 0 and Censor_j == 1:
            return 0
        if Censor_i == 0 and Censor_j == 0:
            return 0
    
    if X_i == X_j:
        if Censor_i == 1 and Censor_j == 1:
            return 0
        if Censor_i == 1 and Censor_j == 0:
            return -1
        if Censor_i == 0 and Censor_j == 1:
            return 1
        if Censor_i == 0 and Censor_j == 0:
            return 0
    
    return 0

def TauXX(timeX, statusX):
    nobs = len(timeX)
    est_tXX = 0
    
    for i in range(nobs):
        for j in range(nobs):
            if j == i:
                continue
            est_tXX += csign(timeX[i], statusX[i], timeX[j], statusX[j]) ** 2
    
    return est_tXX / (nobs * (nobs - 1))

def TauXY(timeX, statusX, scoreY):
    nobs = len(timeX)
    est_tXY = 0
    
    for i in range(nobs):
        for j in range(nobs):
            if j == i:
                continue
            est_tXY += csign(timeX[i], statusX[i], timeX[j], statusX[j]) * np.sign(scoreY[i] - scoreY[j])
    
    return est_tXY / (nobs * (nobs - 1))

def VarTauXX(timeX, statusX):
    nobs = len(timeX)
    temp_s1 = temp_s2 = temp_s3 = 0
    
    for i in range(nobs):
        temp_s1_j = temp_s3_j = 0
        for j in range(nobs):
            if j == i:
                continue
            temp = csign(timeX[i], statusX[i], timeX[j], statusX[j]) ** 2
            temp_s1_j += temp
            temp_s3_j += temp ** 2
        
        temp_s1 += 4 * temp_s1_j ** 2
        temp_s3 += 2 * temp_s3_j
        temp_s2 += temp_s1_j
    
    temp_s2 = temp_s2 / (nobs * (nobs - 1)) * temp_s2 * (2 * nobs - 3) * (-2)
    return (temp_s1 - temp_s3 + temp_s2) / (nobs * (nobs - 1) * (nobs - 2) * (nobs - 3))

def VarTauXY(timeX, statusX, scoreY):
    nobs = len(timeX)
    temp_s1 = temp_s2 = temp_s3 = 0
    
    for i in range(nobs):
        temp_s1_j = temp_s3_j = 0
        for j in range(nobs):
            if j == i:
                continue
            temp = csign(timeX[i], statusX[i], timeX[j], statusX[j]) * np.sign(scoreY[i] - scoreY[j])
            temp_s1_j += temp
            temp_s3_j += temp ** 2
        
        temp_s1 += 4 * temp_s1_j ** 2
        temp_s3 += 2 * temp_s3_j
        temp_s2 += temp_s1_j
    
    temp_s2 = temp_s2 / (nobs * (nobs - 1)) * temp_s2 * (2 * nobs - 3) * (-2)
    return (temp_s1 - temp_s3 + temp_s2) / (nobs * (nobs - 1) * (nobs - 2) * (nobs - 3))

def CovTauXXXY(timeX, statusX, scoreY):
    nobs = len(timeX)
    temp_s1 = temp_s2 = temp_s3 = temp_s4 = 0
    
    for i in range(nobs):
        temp_sXX_j = temp_sXY_j = temp_sXXXY_j = 0
        for j in range(nobs):
            if j == i:
                continue
            temp_XX = csign(timeX[i], statusX[i], timeX[j], statusX[j]) ** 2
            temp_XY = csign(timeX[i], statusX[i], timeX[j], statusX[j]) * np.sign(scoreY[i] - scoreY[j])
            temp_sXX_j += temp_XX
            temp_sXY_j += temp_XY
            temp_sXXXY_j += temp_XX * temp_XY
        
        temp_s1 += 4 * temp_sXX_j * temp_sXY_j
        temp_s3 += 2 * temp_sXXXY_j
        temp_s2 += temp_sXX_j
        temp_s4 += temp_sXY_j
    
    return (temp_s1 - temp_s3 + (2 * nobs - 3) * (-2) * temp_s2 * temp_s4 / (nobs * (nobs - 1))) / (nobs * (nobs - 1) * (nobs - 2) * (nobs - 3))

def CovTauXYXZ(timeX, statusX, scoreY, scoreZ):
    nobs = len(timeX)
    temp_s1 = temp_s2 = temp_s3 = temp_s4 = 0
    
    for i in range(nobs):
        temp_sXY_j = temp_sXZ_j = temp_sXYXZ_j = 0
        for j in range(nobs):
            if j == i:
                continue
            temp_XY = csign(timeX[i], statusX[i], timeX[j], statusX[j]) * np.sign(scoreY[i] - scoreY[j])
            temp_XZ = csign(timeX[i], statusX[i], timeX[j], statusX[j]) * np.sign(scoreZ[i] - scoreZ[j])
            temp_sXY_j += temp_XY
            temp_sXZ_j += temp_XZ
            temp_sXYXZ_j += temp_XY * temp_XZ
        
        temp_s1 += 4 * temp_sXY_j * temp_sXZ_j
        temp_s3 += 2 * temp_sXYXZ_j
        temp_s2 += temp_sXY_j
        temp_s4 += temp_sXZ_j
    
    return (temp_s1 - temp_s3 + (2 * nobs - 3) * (-2) * temp_s2 * temp_s4 / (nobs * (nobs - 1))) / (nobs * (nobs - 1) * (nobs - 2) * (nobs - 3))

