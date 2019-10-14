from __future__ import print_function
import numpy as np
from numpy.linalg import norm
from scipy.sparse.linalg import cg,spsolve,LinearOperator,eigs,LinearOperator,eigsh
from scipy.sparse import eye


def linesearch2(ot,psik,dk,x,w,nbiter):
	Minit,Bar,Costinit,ginit,Hess,t = ot.compute(x,psik,w,hessian=None,parallel=True)
	assert(np.min(Minit)>= 1e-16)
	alpha = 1.0
	nEval = 0
	Mcur,Bar,Costcur,gcur,Hess,t = ot.compute(x,psik+alpha*dk,w,hessian=None,parallel=True)
	nEval += 1
	while norm(ginit) < norm(gcur) or np.min(Mcur) <=0.0 :
		alpha *= 0.5
		Mcur,Bar,Costcur,gcur,Hess,t = ot.compute(x,psik+alpha*dk,w,hessian=None,parallel=True)
		nEval += 1
		if alpha < 1e-5:
			if nbiter:
				return alpha,nEval
			else:
				return alpha
			
	if nbiter:
		return alpha,nEval
	else:
		return alpha

def optimalTransportHess(ot,psiini,x,w,gradTol,iterMax,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	g = gradTol*2 + 1
	gradList = []
	i = 0
	ind = 0
	indmax = 20
	itt = 0
	print('   itt \t      CF             ||Grad||          alpha \t    flag        hidden  \n--------------------------------------------------------------------------------------')
	alpha = -1
	nEval = 0
	M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
	Hess += 1e-6*eye(n)
	gradList.append(norm(g))
	nEval += 1
	print('  %3i    %2.6e      %1.6e     %1.3e'%(itt,Cost,norm(g),alpha))

	while norm(g) > gradTol and itt < iterMax:
		g = g - np.mean(g)
		normG = norm(g)
		# sol = cg(Hess,-g/normG,tol=1e-10,maxiter=3*Hess.shape[0])
		# newtonDirection = sol[0]*normG

		sol = [None,0]
		newtonDirection = spsolve(Hess,-g)

		if(sol[1] != 0):
			print("Warning the CG solve did not converged using iterative solver")
			newtonDirection = spsolve(Hess,-g)

		if nbiter :
			alpha, nEvalLoc = linesearch2(ot,psi,newtonDirection,x,w,nbiter)
			nEval += nEvalLoc
		else:
			alpha = linesearch2(ot,psi,newtonDirection,x,w,nbiter)

		psi += alpha*newtonDirection
		M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
		gradList.append(norm(g))
		nEval += 1
		itt += 1
		print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	# if nbiter :
	# 	return psi,M,Bar,Cost,nEval,gradList
		
	# else:
	# print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	return psi,M,Bar,Cost


def optimalTransportHessINFO(ot,psiini,x,w,gradTol,iterMax,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	g = gradTol*2 + 1
	gradList = []
	fList = []
	alphaList = []
	i = 0
	ind = 0
	indmax = 20
	itt = 0
	print('   itt \t      CF             ||Grad||          alpha \t    flag        hidden  \n--------------------------------------------------------------------------------------')
	alpha = -1
	nEval = 0
	M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
	Hess += 1e-6*eye(n)
	gradList.append(norm(g))
	nEval += 1
	fList += [Cost]
	print('  %3i    %2.6e      %1.6e     %1.3e'%(itt,Cost,norm(g),alpha))

	while norm(g) > gradTol and itt < iterMax:
		g = g - np.mean(g)
		normG = norm(g)
		# sol = cg(Hess,-g/normG,tol=1e-10,maxiter=3*Hess.shape[0])
		# newtonDirection = sol[0]*normG

		sol = [None,0]
		newtonDirection = spsolve(Hess,-g)

		if(sol[1] != 0):
			print("Warning the CG solve did not converged using iterative solver")
			newtonDirection = spsolve(Hess,-g)

		if nbiter :
			alpha, nEvalLoc = linesearch2(ot,psi,newtonDirection,x,w,nbiter)
			nEval += nEvalLoc
		else:
			alpha = linesearch2(ot,psi,newtonDirection,x,w,nbiter)
		alphaList +=[alpha]
		psi += alpha*newtonDirection
		M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
		fList += [Cost]
		gradList.append(norm(g))
		nEval += 1
		itt += 1
		print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	# if nbiter :
	# 	return psi,M,Bar,Cost,nEval,gradList
		
	# else:
	# print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	return psi,M,Bar,Cost,alphaList,fList


def optimalTransportHessRNMINFO(ot,psiini,x,w,gradTol,iterMax,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	g = gradTol*2 + 1
	gradList = []
	fList = []
	alphaList = []
	i = 0
	ind = 0
	indmax = 20
	itt = 0
	print('   itt \t      CF             ||Grad||          alpha \t    flag        hidden  \n--------------------------------------------------------------------------------------')
	alpha = -1
	nEval = 0
	M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
	Hess += 1e-6*eye(n)
	gradList.append(norm(g))
	nEval += 1
	fList += [Cost]
	print('  %3i    %2.6e      %1.6e     %1.3e'%(itt,Cost,norm(g),alpha))

	while norm(g) > gradTol and itt < iterMax:
		g = g - np.mean(g)
		normG = norm(g)
		# sol = cg(Hess,-g/normG,tol=1e-10,maxiter=3*Hess.shape[0])
		# newtonDirection = sol[0]*normG

		sol = [None,0]
		Hess += eye(n)*norm(g)/np.sqrt(n)

		newtonDirection = spsolve(Hess,-g)

		if(sol[1] != 0):
			print("Warning the CG solve did not converged using iterative solver")
			newtonDirection = spsolve(Hess,-g)

		if nbiter :
			alpha, nEvalLoc = linesearch2(ot,psi,newtonDirection,x,w,nbiter)
			nEval += nEvalLoc
		else:
			alpha = linesearch2(ot,psi,newtonDirection,x,w,nbiter)
		alphaList +=[alpha]
		psi += alpha*newtonDirection
		M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
		fList += [Cost]
		gradList.append(norm(g))
		nEval += 1
		itt += 1
		print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	# if nbiter :
	# 	return psi,M,Bar,Cost,nEval,gradList
		
	# else:
	# print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	return psi,M,Bar,Cost,alphaList,fList



def optimalTransportDogLegTR(ot,psi,x,w,gradTol,iterMax=1000,initTRradius=1e-13):
	TRradius = initTRradius
	n = w.size
	sigma1 = 0.2;sigma2 = 0.8;beta1=0.25;beta2=2.0;
	itt = 0
	g = 2*gradTol
	while (norm(g) > gradTol) and (itt < iterMax):
		itt+=1
		M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian='psi')
		# compute the cauchy direction pseudo solution of the sub problem
		# Hess need to be pd thus we remove potentally columns line 
		# import ipdb; ipdb.set_trace()
		maskDog = np.where(M>1e-16)[0]
		print(np.where(M<1e-16)[0])
		HessDog = Hess[np.ix_(maskDog,maskDog)]
		gradDog = g[maskDog]
		# import ipdb; ipdb.set_trace()
		tauk = 1 if gradDog.dot(HessDog*gradDog) < 0 else np.minimum(norm(gradDog)**3.0/(TRradius*gradDog.dot(HessDog*gradDog)),1)
		cauchyDir = - tauk*TRradius*gradDog/norm(gradDog)

		# dk = np.ones(n)*1.0/n
		dk = np.zeros(n)
		dk[maskDog] = cauchyDir

		notCool = True
		alpha = 1.0
		# while notCool:
		# 	M,narv,basa,bsb,bsa,ts = ot.compute(x,psi+alpha*dk,w,hessian=None)
		# 	# import ipdb; ipdb.set_trace()
		# 	if np.where(M<1e-16)[0].size == 0:
		# 		notCool = False
		# 	else:
		# 		alpha /=2.0
		
		# dk = alpha*dk
		psi += alpha*dk

		fkdk = Cost + np.dot(dk,g) + 0.5*dk.dot(Hess*dk)
		M,Bar,fkp1,g,Hess,t = ot.compute(x,psi,w,hessian=None)
		rk = (Cost - fkp1)/(Cost - fkdk)
		# print(rk)
		if rk < sigma1:
			TRradius = beta1*norm(dk)
		if sigma2 <= rk and (norm(dk) - TRradius) < 1e-13:
			TRradius = beta2*TRradius
		print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i TRadius : %1.3e alpha = %1.3e'%(itt,fkp1,norm(g),n-maskDog.size,TRradius,alpha))	
	return psi,M,Bar,fkp1





def optimalTransportMIX(ot,psiini,x,w,gradTol,iterMax):
	n = psiini.size
	psi = psiini
	g = gradTol*2 + 1
	i = 0
	ind = 0
	indmax = 20
	worstcase = 0.0
	while (norm(g) > gradTol) and (i < iterMax and ind < indmax) :
		M,Bar,Cost,g,Hess = ot.compute(x,psi,w,hessian=None)
		
		if (np.min(M) >= 1.0e-16):
			ind += 1
		else:
			ind = 0

		if (i == 0):
			alpha = 1e-3;
		else:
			dg = g - gp
			dpsi = psi - psip
			alpha = np.dot(dpsi,dg)/(norm(dg)**2.0)

		psip = psi.copy()
		gp = g

		psi -= alpha*g
		worstcase = max(np.sum(np.abs(M)<1e-14),worstcase)

		if( i % 10 == 0):
			print('BB MIX -- i:%3i -- Cost:%-02.6e -- ||g||_2:%1.6e -- alpha:%1.2e -- #hidden:%i ' %(i,Cost,norm(g),alpha,worstcase))
			worstcase = 0

		i += 1

	psi += alpha*g

	if (np.min(np.abs(M)) <= 1e-16):
		return psi,M,Bar,Cost
	else:	
		itt = 0
		print('   itt \t      CF             ||Grad||          alpha \t    flag        hidden  \n--------------------------------------------------------------------------------------')
		alpha = -1
		print('  %3i    %2.6e      %1.6e     %1.3e'%(itt,Cost,norm(g),alpha))
		M,Bar,Cost,g,Hess = ot.compute(x,psi,w,hessian='psi')		

		while norm(g) > gradTol and itt < iterMax//50:
			g = g - np.mean(g)
			normG = norm(g)
			sol = cg(Hess,-g/normG,tol=1e-10,maxiter=3*Hess.shape[0])
			newtonDirection = sol[0]*normG

			if(sol[1] != 0):
				print("Warning the CG solve did not converged using iterative solver")
				newtonDirection = spsolve(Hess,-g)

			alpha = linesearch2(ot,psi,newtonDirection,x,w)
			psi += alpha*newtonDirection
			M,Bar,Cost,g,Hess = ot.compute(x,psi,w,hessian='psi')
			itt += 1
			print("  %3i    %2.6e      %1.6e     %1.2e \t \t %i        %i   " % (itt,Cost,norm(g),alpha,sol[1],np.sum(np.abs(M)<1e-14)))
	return psi,M,Bar,Cost

			
def optimalTransportBB(ot,psiini,x,w,gradTol,iterMax):
	psi = psiini
	g = gradTol*2 + 1
	i = 0
	# import ipdb; ipdb.set_trace()
	while (norm(g) > gradTol and i < iterMax) :
		M,Bar,Cost,g,Hess,t = ot.compute(x,psi,w,hessian=None)
		if (i == 0):
			alpha = 1;
			# alpha = 0.002;
			print('Barzilai Borwein -- i:%3i -- Cost:%-02.6e -- ||g||:%1.6e -- alpha:%1.2e -- #hidden:%i ' %(i,Cost,norm(g),alpha,np.sum(np.abs(M)<1e-14)));
		else:
			dg = g - gp
			dpsi = psi - psip
			alpha = np.dot(dpsi,dg)/(norm(dg)**2.0)

		psip = psi.copy()
		gp = g

		psi = psi - alpha*g
		if( i % 1 == 0):
			print('Barzilai Borwein -- i:%3i -- Cost:%-02.6e -- ||g||:%1.6e -- alpha:%1.2e -- #hidden:%i ' %(i,Cost,norm(g),alpha,np.sum(np.abs(M)<1e-14)));

		i += 1

	return psi,M,Bar,Cost



def optimalTransportHessStab(ot,psiini,x,w,gradTol,iterMax):
	n = psiini.size
	psi = psiini
	M,Bar,Cost,Grad,Hess = ot.compute(x,psi,w,hessian='psi')		
	itt = 0
	NonEmpty = np.where(M>0)[0]
	nbEmpty = n-NonEmpty.size
	alpha = 1.
	NonConverged=False
	print('--iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  ||Grad||1 : %1.6e nbHidden : %3i '%(itt,Cost,norm(Grad),norm(Grad,ord=1),nbEmpty))	

	while norm(Grad) > gradTol and itt < iterMax:
		g = Grad[NonEmpty]
		Mean = np.mean(g)
		g -= Mean
		normG = norm(g)
		choice='Newton' # by default choose a Newton update
		if NonConverged :
			choice='Newton'
			alpha=1.
			NonConverged=False
		if norm(g) < 1.e-2*Mean*np.sqrt(g.size) :
			# the vector g is too close from its mean value, pick the gradient update
			choice='Gradient'
		if choice=='Newton' :
			def mv(x): # Define the matrix-vector product corresponding to the extraction of the Hessian
				z=np.zeros(n)
				z[NonEmpty]=x
				tmp=Hess.dot(z)
				return tmp[NonEmpty]
			A=LinearOperator(dtype=float,shape=(NonEmpty.size,NonEmpty.size),matvec=mv)			
			sol = cg(A,-g/normG,tol=1e-5) # solve the linear system on the dof nonEmpty
			SolCG = sol[0]*normG #renormalization to avoid the absolute convergence of cg
			#if (sol[1] != 0) : # Newton method did not converge, revert to Gradient method
			#	choice='Gradient'
			if choice=='Newton' : 
				Dir = np.zeros(n)
				Dir = - Grad
				Dir[NonEmpty] = SolCG - np.mean(SolCG) - Mean
				alpha=min(1.,2*alpha)
		if choice=='Gradient'	: # We choose the gradient update
			alpha=2*alpha
			Dir = - Grad
		Mtmp,Bar,Costtmp,gtmp,Hess = ot.compute(x,psi+alpha*Dir,w,hessian='psi')
		NonEmptytmp = np.where(Mtmp>0)[0]
		nbEmptytmp = n-NonEmptytmp.size
		while (norm(gtmp) > norm(Grad)*(1+1.e-5)) or (nbEmptytmp > nbEmpty) :
			alpha *= 0.5
			Mtmp,Bar,Costtmp,gtmp,Hess = ot.compute(x,psi+alpha*Dir,w,hessian='psi')
			NonEmptytmp = np.where(Mtmp>0)[0]
			nbEmptytmp = n-NonEmptytmp.size
			if alpha < 1e-6:
				if choice=="Newton" : alpha=0.
				NonConverged=True
				
		psi += alpha*Dir
		M=Mtmp
		Cost=Costtmp
		Grad=gtmp
		NonEmpty = NonEmptytmp
		nbEmpty = nbEmptytmp
		itt += 1
		if( itt % 1 == 0):
				print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  ||Grad||1 : %1.6e nbHidden : %3i alpha : %1.3e'%(itt,Cost,norm(Grad),norm(Grad,ord=1),nbEmpty,alpha)+' '+choice)	
		#assert False

	return psi,M,Bar,Cost


def optimalTransportHessLM(ot,psiini,x,w,gradTol,iterMax,verbose=True,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	# import matplotlib.pyplot as plt 
	# plt.plot(x[:,0],x[:,1],'ro',ms=2)
	# plt.show()
	# import pdb
	# pdb.set_trace()
	
	gradList = []
	nEval = 0
	M,Bar,Cost,Grad,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
	nEval += 1
	c=1.e-2*np.max(np.abs(Hess.diagonal()))		
	itt = 0
	NonEmpty = np.where(M>0)[0]
	nbEmpty = n-NonEmpty.size
	if verbose :
		print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i c : %1.3e'%(itt,Cost,norm(Grad),nbEmpty,c))	
		gradList.append(norm(Grad))

	while (norm(Grad) > gradTol or nbEmpty >0 ) and itt < iterMax   :
		import time; ti = time.time()
		# c/=2.7
		c/=1.7
		Hess2=Hess+c*eye(n)
		Dir = spsolve(Hess2,-Grad) # solve the linear system on the dof nonEmpty
		Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi+Dir,w,hessian='psi',parallel=parallel)
		nEval += 1
		NonEmptytmp = np.where(Mtmp>0)[0]
		nbEmptytmp = n-NonEmptytmp.size
		epsilon=1.e-1
		predicted_decrease=np.dot(Dir,Grad)
		#print(norm(gtmp),norm(Grad),Costtmp,Cost,Cost+epsilon*predicted_decrease)
		while ((norm(gtmp) > norm(Grad)) and (Costtmp > Cost+epsilon*predicted_decrease) ) or (nbEmptytmp > nbEmpty) :
			c = 3.1*c
			# c = 3.4*c
			Hess2=Hess+c*eye(n)
			Dir = spsolve(Hess2,-Grad) # solve the linear system on the dof nonEmpty
			Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi+Dir,w,hessian='psi',parallel=parallel)
			nEval += 1
			NonEmptytmp = np.where(Mtmp>0)[0]
			nbEmptytmp = n-NonEmptytmp.size
			if c > 1e15:
				Dir =np.zeros(Dir.shape)
				itt=iterMax
				Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
				nEval += 1
				print('WARNING OF A BREAK')
				if nbiter :
					return psi,M,Bar,Cost,nEval,gradList
				else :
					return psi,M,Bar,Cost
		psi += Dir
		M=Mtmp
		Cost=Costtmp
		Grad=gtmp
		Hess=Hesstmp
		NonEmpty = NonEmptytmp
		nbEmpty = nbEmptytmp
		itt += 1
		if (verbose) and ( itt % 1 == 0):
				print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i c : %1.3e s=%3.2f'%(itt,Cost,norm(Grad),nbEmpty,c,time.time()-ti))	
				gradList.append(norm(Grad))
	
	if nbiter :
		return psi,M,Bar,Cost,nEval,gradList
	else :
		return psi,M,Bar,Cost



def optimalTransportRNM(ot,psiini,x,w,gradTol,iterMax,verbose=True,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	# import matplotlib.pyplot as plt 
	# plt.plot(x[:,0],x[:,1],'ro',ms=2)
	# plt.show()
	# import pdb
	# pdb.set_trace()
	
	gradList = []
	nEval = 0
	M,Bar,Cost,Grad,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
	nEval += 1
	c=1.e-2*np.max(np.abs(Hess.diagonal()))		
	itt = 0
	NonEmpty = np.where(M>0)[0]
	nbEmpty = n-NonEmpty.size
	if verbose :
		print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i c : %1.3e'%(itt,Cost,norm(Grad),nbEmpty,c))	
		gradList.append(norm(Grad))

	while (norm(Grad) > gradTol or nbEmpty >0 ) and itt < iterMax   :
		import time; ti = time.time()
		# c/=2.7
		c/=1.7
		Hess2=Hess+np.linalg.norm(Grad)*eye(n)
		Dir = spsolve(Hess2,-Grad) # solve the linear system on the dof nonEmpty
		DirFull = Dir.copy()
		Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi+Dir,w,hessian='psi',parallel=parallel)
		nEval += 1
		NonEmptytmp = np.where(Mtmp>0)[0]
		nbEmptytmp = n-NonEmptytmp.size
		epsilon=1.e-1
		predicted_decrease=np.dot(Dir,Grad)
		#print(norm(gtmp),norm(Grad),Costtmp,Cost,Cost+epsilon*predicted_decrease)
		TookGrad = False
		while ((norm(gtmp) > norm(Grad)) and (Costtmp > Cost+epsilon*predicted_decrease) ) or (nbEmptytmp > nbEmpty) :			
			if np.linalg.norm(Dir)/np.linalg.norm(DirFull) < 5e-5:
				Dir = -5e-5*Grad
				TookGrad = True
				break
			Dir *= .5 
			Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi+Dir,w,hessian=None,parallel=parallel)
			nEval += 1
			NonEmptytmp = np.where(Mtmp>0)[0]
			nbEmptytmp = n-NonEmptytmp.size
		psi += Dir
		Mtmp,Bar,Costtmp,gtmp,Hesstmp,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
		nEval += 1
		M=Mtmp
		Cost=Costtmp
		Grad=gtmp
		Hess = Hesstmp
		NonEmpty = NonEmptytmp
		nbEmpty = nbEmptytmp
		itt += 1
		if (verbose) and ( itt % 1 == 0):
				print('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i step : %1.3e s=%3.2f N : %i'%(itt,Cost,norm(Grad),nbEmpty,np.linalg.norm(Dir)/np.linalg.norm(DirFull),time.time()-ti,0 if TookGrad else 1))	
				gradList.append(norm(Grad))
	
	if nbiter :
		return psi,M,Bar,Cost,nEval,gradList
	else :
		return psi,M,Bar,Cost


def optimalTransportGRNM(ot,psiini,x,w,gradTol,iterMax,verbose=True,nbiter=False,parallel=False):
	n = psiini.size
	psi = psiini
	epsilon = gradTol
	M,Bar,Cost,Grad,Hess,t = ot.compute(x,psi,w,hessian=None,parallel=parallel)
	phi = Cost
	s = 0.0
	l = 0.0
	phiHat = phi
	xHat = psi
	k = .5
	m0 = 1e-3
	M0 = 1e4
	k0 = m0*k
	mBar = .2
	MBar = 100
	goTo3 = False
	goTo7 = False
	itt = -1
	# print(k0)
	# assert(False)

	while True:
		itt+=1
		M,Bar,Cost,Grad,Hess,t = ot.compute(x,psi,w,hessian=None,parallel=parallel)
		NonEmpty = np.where(M>0)[0]
		nbEmpty = n-NonEmpty.size
		print(('  iteration : %3i  CF : %2.6e  ||Grad||2 : %1.6e  nbHidden : %3i'%(itt,Cost,norm(Grad),nbEmpty)),end=' ')	
		if np.linalg.norm(Grad) < epsilon and not goTo3 :
			# print("step 1:")
			return psi,M,Bar,Cost,nEval,gradList
		if np.linalg.norm(Grad) < k0 and not goTo3 :
			# print("step 2:")
			goTo7 = True
		else :
			goTo7 = False
			# print("step 3:")
			# A = 0*Hess
			s += 1.
			t = 1./s
			t = 1e-4
			# descente de gradient normalisee
			# print("step 4:")
			xHat = (psi -t*Grad/np.linalg.norm(Grad)).copy()
			M,Bar,Cost,Grad,Hess,t = ot.compute(x,xHat,w,hessian=None,parallel=parallel)
			if Cost < phiHat :
				print("step 5:")
				phiHat = Cost
				goTo3 = False
				psi = xHat.copy()
			else :
				print("step 6:")
				psi = xHat.copy()
				goTo3 = True

		if goTo7:
			M,Bar,Cost,Grad,Hess,t = ot.compute(x,psi,w,hessian='psi',parallel=parallel)
			Hess2=Hess+np.linalg.norm(Grad)*eye(n)
			xHat =  (psi + spsolve(Hess2,-Grad)).copy()
			M,Bar,Cost,Grad,Hess,t = ot.compute(x,xHat,w,hessian=None,parallel=parallel)
			if np.linalg.norm(Grad) > np.linalg.norm(Grad)**1.5 :
				t = 0.5*m0/M0
				l+= 1
				m0 = mBar*l**(-.1)
				M0 = MBar*l**(.1)
				# k0 = m0*k
				xHat = (psi -t*Grad/np.linalg.norm(Grad)).copy()
				print("step 7: k0= "+str(k0))
			else:
				print("step 8:")
			M,Bar,Cost,Grad,Hess,t = ot.compute(x,xHat,w,hessian=None,parallel=parallel)
			phiHat = Cost
			psi = xHat.copy()
	
	if nbiter :
		return psi,M,Bar,Cost,nEval,gradList
	else :
		return psi,M,Bar,Cost