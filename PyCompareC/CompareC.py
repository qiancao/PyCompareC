import numpy as np
from survivalC import TauXX, TauXY, VarTauXX, VarTauXY, CovTauXXXY, CovTauXYXZ

def estC(timeX, statusX, scoreY):
    tau_xx = TauXX(timeX, statusX)
    tau_xy = TauXY(timeX, statusX, scoreY)
    return (tau_xy / tau_xx + 1) / 2

def vardiffC(timeX, statusX, scoreY, scoreZ):
    t11 = TauXX(timeX, statusX)
    t12 = TauXY(timeX, statusX, scoreY)
    t13 = TauXY(timeX, statusX, scoreZ)
    
    var_t11 = VarTauXX(timeX, statusX)
    var_t12 = VarTauXY(timeX, statusX, scoreY)
    var_t13 = VarTauXY(timeX, statusX, scoreZ)
    
    cov_t1112 = CovTauXXXY(timeX, statusX, scoreY)
    cov_t1113 = CovTauXXXY(timeX, statusX, scoreZ)
    cov_t1213 = CovTauXYXZ(timeX, statusX, scoreY, scoreZ)
    
    matrix_1 = np.array([[var_t12, cov_t1112], [cov_t1112, var_t11]])
    vector_1 = np.array([1/t11, -t12/t11**2])
    est_varCxy = 0.25 * np.dot(vector_1, np.dot(matrix_1, vector_1))
    
    matrix_2 = np.array([[var_t13, cov_t1113], [cov_t1113, var_t11]])
    vector_2 = np.array([1/t11, -t13/t11**2])
    est_varCxz = 0.25 * np.dot(vector_2, np.dot(matrix_2, vector_2))
    
    matrix_3 = np.array([[cov_t1213, cov_t1113], [cov_t1112, var_t11]])
    vector_3 = np.array([1/t11, -t13/t11**2])
    est_cov = 0.25 * np.dot(vector_1, np.dot(matrix_3, vector_3))
    
    est_vardiff_c = est_varCxy + est_varCxz - 2 * est_cov
    
    return {
        "est.vardiff_c": est_vardiff_c,
        "est.varCxy": est_varCxy,
        "est.varCxz": est_varCxz,
        "est.cov": est_cov
    }

def compareC(timeX, statusX, scoreY, scoreZ):
    est_c = np.array([estC(timeX, statusX, scoreY), estC(timeX, statusX, scoreZ)])
    est_diff_c = est_c[0] - est_c[1]
    
    tmpout = vardiffC(timeX, statusX, scoreY, scoreZ)
    zscore = est_diff_c / np.sqrt(tmpout["est.vardiff_c"])
    pval = 2 * (1 - (0.5 * (1 + np.math.erf(abs(zscore) / np.sqrt(2)))))
    
    return {
        "est.c": est_c,
        "est.diff_c": est_diff_c,
        "est.vardiff_c": tmpout["est.vardiff_c"],
        "est.varCxy": tmpout["est.varCxy"],
        "est.varCxz": tmpout["est.varCxz"],
        "est.cov": tmpout["est.cov"],
        "zscore": zscore,
        "pval": pval
    }

