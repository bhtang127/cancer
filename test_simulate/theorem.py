import numpy as np
import matplotlib.pyplot as plt

class GeneType():
    
    def __init__(self, Gene, alpha):
        # alpha = {division_rate:division, pasym, psym, dsym, 
        #          death_rate:death, v, compete:c, K, adjust,
        #          mul_division:dd, delta_psym:dp, delta_mute:dm,
        #          l_ratio:l/c1+l, pow:alpha in growth func}
        self.G = ""
        self.alpha = alpha
        self.mutate(Gene)
        return
    
    def mutate(self, types):
        
        alpha = {key:self.alpha[key] for key in self.alpha}
        for mut in types:
            
            if len(self.G) >= 3:
                return self.get_attr(), self.get_alpha()
        
            if mut == "S":
                alpha["division"] *= alpha["dd"]
            elif mut == "F":
                r_before = (alpha["psym"]-alpha["dsym"])*alpha["division"]-alpha["death"]
                
                alpha["psym"] += alpha["dp"]
                alpha["pasym"] -= alpha["dp"]
                
                r_after = (alpha["psym"]-alpha["dsym"])*alpha["division"]-alpha["death"]
                
            else:
                alpha["v"] += alpha["dm"]

            self.G += mut
        
        self.alpha = alpha
        return self.get_attr(), self.get_alpha()
    
    def get_attr(self):
        alpha = self.alpha
        
        self.va = alpha["v"] * alpha["adjust"] / 2
        self.vs = alpha["v"] * alpha["adjust"]
        self.v = self.va * alpha["pasym"] + self.vs * alpha["psym"]
        self.K = alpha["K"]
    
        self.b = self.alpha["psym"] * self.alpha["division"]
        self.d = self.alpha["dsym"] * self.alpha["division"] + self.alpha["death"]
        self.r = self.b - self.d
        
        return {"b":self.b, "d":self.d, "r":self.r, "K":self.K, "G":self.G,
                "va":self.va, "vs":self.vs, "v":self.v}
    
    def get_alpha(self):
        return self.alpha

class growth():
    
    def __init__(self, beta):
        # beta is alpha in main class
        self.beta = beta
    
    def growth_rate(self,x):
        beta = self.beta
        l = beta["division"]*beta["l_ratio"]
        c1 = beta["division"]*(1-beta["l_ratio"])*np.exp(1)
        
        if x > beta["K"]:
            div = l
        else:
            div = c1*np.exp(-1/(1-(x/beta["K"])**beta["pow"])) + l
            
        return beta["psym"] * div - beta["dsym"] * div - beta["death"]
    
    def get_init_time(self):
        beta = self.beta
        x_size = 1; t = 0.1
        while True:
            x_size *= np.exp(self.growth_rate(x_size)*0.1)
            if x_size >= beta["eps"]*beta["K"]:
                return t
            t += 0.1
    
    def get_dynamics(self, ts):
        beta = self.beta
        x_size = beta["eps"] * beta["K"]
        dt = ts[1]-ts[0]
        res = np.zeros_like(ts)
        for i, t in enumerate(ts):
            res[i] = x_size
            x_size *= np.exp(self.growth_rate(x_size)*dt)
        return res
    
    
class tissue():
    
    def __init__(self, alpha, pfunc, 
                 T, period, Gene="",
                 check = True):
        self.eps = alpha["eps"]
        self.ptrans = pfunc
        self.alpha = alpha
        self.T = T
        self.period = period
        self.para = GeneType(Gene, alpha).get_attr()
        self.ts = np.linspace(0, self.T, num=self.period)
        self.solution = {}
        if check:
            self.self_check()
    
    def self_check(self):
        self.excess_survive_rate = (self.para["d"]/self.para["b"]) * (self.eps*self.para["K"])**(-1)
        self.mutation_before_epsK = 1 - np.exp(
                        -self.para["v"]*self.eps*self.para["K"]/self.para["r"] * 
                         (np.log(self.para["r"]*self.eps*self.para["K"]/self.para["b"]) + 0.58)
                        )
        print("Excess Survival ratio:", self.excess_survive_rate/(self.para["r"]/self.para["b"]))
        print("P(mutiation in [0,tao]):", self.mutation_before_epsK)
#         init_T = growth(self.alpha).get_init_time()
        
#         print("Approximation Consistency: ", 
#               1-np.exp(-self.para["r"]*self.eps*self.para["K"]*\
#                        np.exp(-self.para["r"]*init_T)/self.para["b"]))
        
#         print("Init time to grow to size eps*K: {0:.2f} years".format(init_T/52))
#         self.T = self.T - growth(self.alpha).get_init_time()
#         self.ts = np.linspace(0, self.T, num=self.period)
    
    def density(self, ts, para):
        e_rGt = np.exp(-para["r"]*ts)
        return (para["r"]**2 * self.eps*para["K"])/para["b"] * e_rGt * \
                np.exp(-para["r"]*self.eps*para["K"] * e_rGt/para["b"])
    
    def get_Gene(self, G, alpha):
        return GeneType(G, alpha)
    
    def mutate(self, types):
        self.para, self.alpha = self.get_Gene(self.para["G"],self.alpha).mutate(types)
    
    def vGfunc(self, ts, alpha, para):
        return para["v"] * growth(alpha).get_dynamics(ts)
        
    def muiGUfunc(self, ts, mtype, para, alpha):
        GU, _ = self.get_Gene(para["G"], alpha).mutate(mtype)
        vG = self.vGfunc(ts, alpha, para)
        fGU = self.density(ts, GU)
        ds = (ts[-1] - ts[0]) / (len(ts)-1)
        coeff = self.ptrans(para["G"], mtype)*GU["r"] / GU["b"]
        return coeff * np.convolve(vG*ds, fGU)[:len(ts)]
    
    def evolve(self, ts, muiGU, qGU):
        ds = (ts[-1] - ts[0]) / (len(ts)-1)
        qGU_rev = qGU[::-1]
        aux = np.convolve(muiGU*ds, 1-qGU_rev)[:len(muiGU)]
        aux = aux[::-1]
        return np.exp(-aux)
        
    def solveAUX(self, ts, para, alpha):
        if para["G"]:
            print(para["G"], end="  ")
#         if len(para["G"]) == 3 and "S" in para["G"] and "F" in para["G"]:
#             self.solution[para["G"]] = np.zeros_like(ts)
#             return np.zeros_like(ts)
#         elif len(para["G"]) == 3:
#             self.solution[para["G"]] = np.ones_like(ts)
#             return np.ones_like(ts)
        if len(para["G"]) == 3:
            self.solution[para["G"]] = np.zeros_like(ts)
            return np.zeros_like(ts)
        else:
            qG = np.ones_like(ts)
            for U in "SFM":
                paraU,alphaU = self.get_Gene(para["G"], alpha).mutate(U)
                muiGU = self.muiGUfunc(ts, U, para, alpha)
                qG *= self.evolve(ts, muiGU, self.solveAUX(ts, paraU, alphaU))
            self.solution[para["G"]] = qG
            return qG
    
    def solve(self):
        print("Begin Solving")
        self.solveAUX(self.ts,self.para,self.alpha)
        print("")
        print("Done")
        return self.solution
        
    
def pGUfunc(G, U):
    if len(G) == 0:
        pdict = {"F":1, "S":0, "M":0}
    elif len(G) == 1:
        pdict = {"F":645, "S":643, "M":1197}
    elif len(G) == 2:
        if G[1] == "F":
            pdict = {"F":2, "S":19, "M":15}
        if G[1] == "S":
            pdict = {"F":5, "S":21, "M":30}
        if G[1] == "M":
            pdict = {"F":2, "S":9, "M":12}
    else:
        pdict = {"F":1, "S":1, "M":1}
        
    return pdict[U]/(pdict["S"]+pdict["F"]+pdict["M"])

def pGUfunc2(G, U):
    F, S, M = 0, 0, 0
    for c in G:
        if c == "F":
            F += 1
        if c == "S":
            S += 1
    pdict = {"F":4/(1+F), "S":7/(1+S), "M":4}
#     pdict = {"F":4, "S":7, "M":7}
    return pdict[U]/(pdict["F"]+pdict["S"]+pdict["M"])