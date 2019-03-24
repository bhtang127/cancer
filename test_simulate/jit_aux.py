# from numba import jit

# @jit
def Likelihoods(D,Y):
    l=np.log(sp.stats.norm.cdf(Y*D))
    dl=Y*sp.stats.norm.pdf(D)/sp.stats.norm.cdf(Y*D)
    ddl=-sp.stats.norm.pdf(D)**2/sp.stats.norm.cdf(Y*D)**2-Y*D*sp.stats.norm.pdf(D)/sp.stats.norm.cdf(Y*D)
    return l, dl, np.diag(ddl)

# @jit
def MaxLikelihood(K,Y,tol):
    n=np.shape(Y)[0]
    D=np.zeros(n)
    #DOld=np.ones(n)
    Lst=[]
    grads=[]
    values=[]
    gr=np.ones(n)
    Ds=[]
    k=0
    while (np.linalg.norm(gr)>tol) and k<1000:
        k +=1
        #DOld=D
        Lik=Likelihoods(D,Y)
        W=-Lik[2]
        WHalf=sp.linalg.sqrtm(W)
        L=np.linalg.cholesky(np.diag(np.ones(n))+WHalf @ K @ WHalf)
        b=W @ D+Lik[1]
        d=np.linalg.solve(L,WHalf @ K @ b)   
        a=b- WHalf @ (np.linalg.solve(L.transpose(),d))
        gr=Lik[1]-np.linalg.solve(K+0.0000001*np.diag(np.ones(n)),D)
        grads.append(np.linalg.norm(gr))
        D=np.dot(K,a)                      
        values.append(-1/2*np.dot(a,D)+np.sum(Lik[0]))
        Lst.append(-1/2*np.dot(a,D)+np.sum(Lik[0])-np.sum(np.diag(L)))
        Ds.append(D)
    return D,Lst,grads,values, Ds

# @jit
def EPApprox(K,Y,tol):
    n=np.shape(Y)[0]
    nuTilde=np.zeros(n)
    tauTilde=np.zeros(n)
    nuOld=np.ones(n)
    tauOld=np.ones(n)
    Sigma=np.copy(K)
    mu=np.zeros(n)
    k=0
    while ((np.linalg.norm(nuTilde-nuOld,np.inf)+np.linalg.norm(tauTilde-tauOld,np.inf))>tol) and k<50:
            nuOld=np.copy(nuTilde)
            tauOld=np.copy(tauTilde)
            k +=1
            
            for i in range(n):
                tauMinusI=Sigma[i,i]**(-2)-tauTilde[i]
                nuMinusI=Sigma[i,i]**(-2)*mu[i]-tauTilde[i]
                muMinusI=nuMinusI/tauMinusI
                sigmaMinusI2=1/tauMinusI
                zI=Y[i]*muMinusI/np.sqrt(1+sigmaMinusI2)
                #ZhatI=sp.stats.norm.cdf(zI)
                muHatI=muMinusI+Y[i]*sigmaMinusI2*sp.stats.norm.pdf(zI)/(sp.stats.norm.cdf(zI)*np.sqrt(1+sigmaMinusI2))
                sigmaHatI2=sigmaMinusI2-(zI+sp.stats.norm.pdf(zI)/sp.stats.norm.cdf(zI))*(sigmaMinusI2**2*sp.stats.norm.pdf(zI)/((1+sigmaMinusI2)*sp.stats.norm.cdf(zI)))
                deltaTauTilde=1/sigmaHatI2-tauMinusI-tauTilde[i]
                tauTilde[i]=tauTilde[i]+deltaTauTilde
                #print(deltaTauTilde)
                nuTilde[i]=(1/sigmaHatI2)*muHatI-nuMinusI
                sI=Sigma[:,i].reshape((n,1))
                Sigma -=(1/((1/deltaTauTilde)+Sigma[i,i]))*(sI @ sI.transpose())
                mu= Sigma @ nuTilde
            print(((np.linalg.norm(nuTilde-nuOld,np.inf)+np.linalg.norm(tauTilde-tauOld,np.inf))))
            STilde=np.diag(tauTilde)
            STildeHalf=sp.linalg.sqrtm(STilde)
            L=np.linalg.cholesky(np.diag(np.ones(n))+STildeHalf @ K @ STildeHalf)
            V=np.linalg.solve(L.transpose(),STildeHalf @ K)
            Sigma=K-((V.transpose()) @ V)
            mu= Sigma @ nuTilde
    return nuTilde,tauTilde,k
            
# @jit
def PosteriorMeanVariance(K,Y,tol,x,X,C,l):  
    d=np.shape(x)[0]           
    N=np.shape(Y)[0]
    #K=np.zeros
    #K=np.zeros((N,N))
    Kstar=np.zeros(N)
    Kstarprime=np.zeros((N,d))
#        for i,j in product(range(N),range(N)):
#            K[i,j]=C**2*np.exp(-1/(2*l**2)*(np.linalg.norm(X[i,:]-X[j,:]))**2)
    DHat=MaxLikelihood(K+0.0*np.diag(np.ones(N)),Y,tol)         
    for i in range(N):
        Kstar[i]=C**2*np.exp(-1/(2*l**2)*(np.linalg.norm(X[i,:]-x))**2) 
    for i,j in product(range(N),range(d)):
        Kstarprime[i,j]=-1/(l**2)*Kstar[i]*(x[j]-X[i,j])                          
    kStar=C**2
    Lik=Likelihoods(DHat[0],Y)
    PMean=np.dot(Kstar,Lik[1])
    PmeanPrime=np.dot(Kstarprime.transpose(),Lik[1])
    W=-Lik[2]
    L=np.linalg.cholesky(np.diag(np.ones(N))+np.sqrt(W) @ K @ np.sqrt(W))
    #A=np.linalg.solve((K-np.linalg.inv(Lik[2])),Kstar)
    v=np.linalg.solve(L,np.sqrt(W) @ Kstar)
    vPrime=np.linalg.solve(L, np.sqrt(W) @ Kstarprime)
    PVariance=kStar-np.dot(v,v)
    PVariancePrime=-2* (vPrime.transpose() @ v)
    ObjectivePrime=sp.stats.norm.pdf(PMean/np.sqrt(1+PVariance))*(PmeanPrime*(1+PVariance)**(-1/2)-1/2*PMean*PVariancePrime*(1+PVariance)**(-3/2))
    return PMean,PVariance,sp.stats.norm.cdf(PMean/np.sqrt(1+PVariance)), ObjectivePrime,K


# @jit
def FitGradient(L_K,Y,tol,X,C,l,Const,LearningRate,iterations):
    N=np.shape(Y)[0]
    d=np.shape(Y)[1]
    x=np.copy(X[N-2,:])
    #xold=np.ones(d)
    Jacobian=np.zeros((d,d))
    values=np.zeros(d)
    k=0
    while  k<iterations:
            k +=1
            #print(k)
            print((values**2).sum())
            #xold=np.copy(x)
            for i in range(d):
                Post=PosteriorMeanVariance(L_K[i],Y[:,i],tol,x,X,C[i],l[i])
                values[i]=Post[2]-Const[i]
                Jacobian[i,:]=Post[3]
                #print (Jacobian,np.linalg.det(Jacobian))
            gradient=np.sum(np.diag(values) @ Jacobian ,axis=0)    
            x -=(2*LearningRate)*gradient
            
            x[0]=max(0.0501,min(x[0],0.06))/0.0501
            x[1]=max(1,min(x[1],2))
            x[2]=max(1e-3,min(x[2],2e-3))*1e3
    return x, values,gradient
        
# @jit
def ChooseBestParam(Y,tol,X,X_new,C,l,Const):
    d=np.shape(Y)[1]           
    N=np.shape(Y)[0]
    Perf=np.zeros(N)
    for i in range(N):
        for j in range(d):
            Post=PosteriorMeanVariance(Y[:,j],tol,X_new[i,:],X,C[j],l[j])
            Perf[i] +=np.abs(Post[2]-Const[j])
    index=np.argmin(Perf)
    return index,X_new[index,:],Perf[index]
        
    