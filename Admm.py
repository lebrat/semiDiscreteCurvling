import scipy.sparse as sp
import numpy as np
from scipy.sparse.linalg import factorized,cg,spsolve
import matplotlib.pyplot as plt
from  time import time

class admm() :
	def __init__(self,shape,As,projector,gamma,B=None,metric=None,verificators=None,solver='factorize') :
		self.metric=sp.eye(shape[0]) if metric is None else metric
		self.verificators=len(As)*[None] if verificators is None else verificators
		self.solver=solver
		self.check_init(shape,As,projector,gamma,B,self.metric,self.verificators,self.solver)
		self.Alist=As
		self.gamma=gamma
		self.ATlist=[A.T.tocsr() for A in self.Alist]
	#	print('transposition done')
		M=self.metric.copy()
		for (A,AT,g) in zip(self.Alist,self.ATlist,self.gamma) :
			M=M+g*g*AT.dot(A)
	#	print('Matrix M assembly done')
		self.M= M if B==None else sp.bmat([[M,B.T],[B,None]])
		if self.solver=='factorize' :
			self.M_factorized=factorized(self.M.tocsc())
		self.Alist=[A.tocsr() for A in self.Alist]
	#	print('Matrix M factorization done')
		self.projector=projector
		self.shape_x=(shape[0],shape[1])
		self.indices_y=[]
		begin=0
		for  A in self.Alist :
			self.indices_y.append((begin,begin+A.shape[0]))
			begin+=A.shape[0]
		self.shape_y=(begin,self.shape_x[1])

	def gA(self,x) : 
		y=np.zeros(self.shape_y)
		for (ind,A,g) in zip(self.indices_y,self.Alist,self.gamma) :
			for i in range(self.shape_y[1]) :
				y[ind[0]:ind[1],i]=g*A.dot(x[:,i])
		return y	
		
	def gAT(self,y) : 
		x=np.zeros(self.shape_x)
		for (ind,AT,g) in zip(self.indices_y,self.ATlist,self.gamma) :
			for i in range(self.shape_x[1]) :			
				x[:,i]+=g*AT.dot(y[ind[0]:ind[1],i])
		return x
		
	def solveM(self,z,b) :
		x=np.zeros(self.shape_x)
		mu=None if b is None else np.zeros(b.shape) 
		for i in range(self.shape_x[1]) :
			if self.solver=='factorize':			
				tmp=self.M_factorized(z[:,i]) if b is None else self.M_factorized(np.concatenate((z[:,i],b[:,i])))
			elif self.solver=='spsolve' :
				tmp=spsolve(self.M,z[:,i]) if b is None else spsolve(self.M,np.concatenate((z[:,i],b[:,i])))
			x[:,i]=tmp[:self.shape_x[0]]
			if b is not None :
				mu[:,i]=tmp[self.shape_x[0]:]
		return (x,mu)
		
	def project(self,z,list_bnd) :
		y=np.zeros(self.shape_y)
		for (ind,projector,g,bnd) in zip(self.indices_y,self.projector,self.gamma,list_bnd) :
			y[ind[0]:ind[1],:]=projector(z[ind[0]:ind[1],:],g,bnd)
		return y
		
	def verif(self,z) :
		result=[]
		for (ind,verif,g) in zip(self.indices_y,self.verificators,self.gamma) :
			if verif is not None :
				result.append(verif(z[ind[0]:ind[1]],g))
			else :
				result.append(None)
		return result
		
			
	def solve(self,z,list_bound,b=None,nitermax=100,xInit=None,lagr=None,tol=1.e-6,verbose=0,verificators=None) :
		niter=0
		lagrange=np.zeros(self.shape_y) if lagr is None else lagr
		x=np.zeros(z.shape) if xInit is None else np.copy(xInit)
		self.check_solve(z,list_bound,b,xInit,lagrange)
		tmp0=self.gA(z)
		tmp=self.project(tmp0.copy(),list_bound)

		if np.linalg.norm(tmp-tmp0) <tol*np.linalg.norm(tmp0) :
			x=np.copy(z)
			tab = self.verif(tmp0)
			tab=[t/l*100 for (t,l) in zip(tab,list_bound)]
			if verbose > 0:
				print('Immediate convergency |z-proj(z)| : %1.2e  |z| : %1.2e, bnds=[%2.1f%% (%1.2e),%2.1f%%  (%1.2e)] '%(np.linalg.norm(tmp-tmp0),np.linalg.norm(tmp0),tab[0],tab[1],list_bound[0],list_bound[1]))
			return z,(True,0,np.zeros(self.shape_y),tab)
		if verbose > 1 :
			tab = self.verif(self.gA(x))
			tab=[t/l*100 for (t,l) in zip(tab,list_bound)]
			print(' niter : %4i  distance : %1.2e --- evol(L) :%1.2e ---- evol(x) :%1.2e bnds=[%2.1f%% (%1.2e),%2.1f%%  (%1.2e)]'%(niter, np.linalg.norm(x-z),0.,0.,tab[0],list_bound[0],tab[1],list_bound[1]))
		while True :				
				niter+=1
			 #update y
				y=self.project(self.gA(x) +lagrange,list_bound)
				#print y,self.gA(x),list_bound
				#return
				#update x
				xOld=np.copy(x)
				x,mu=self.solveM(self.gAT(y - lagrange)+self.metric.dot(z),b)

				#update lagrange
				lOld=np.copy(lagrange)
				lagrange+=self.gA(x)-y
				# calcul des criteres de convergence
			
				erreurl=np.linalg.norm(lagrange-lOld)/(np.linalg.norm(lagrange)+1.e-12)				
				erreurx=np.linalg.norm(x-xOld)/(np.linalg.norm(x)+1.e-12)
				tab = self.verif(self.gA(x))
				tab=[t/l*100 for (t,l) in zip(tab,list_bound)]
				if verbose > 2 :
					print(' niter : %4i  distance : %1.2e --- evol(L) :%1.2e ---- evol(x) :%1.2e bnds=[%2.1f%% (%1.2e),%2.1f%%  (%1.2e)]'%(niter, np.linalg.norm(x-z),erreurl,erreurx,tab[0],list_bound[0],tab[1],list_bound[1]))
				if verbose >1:
					if (niter % (nitermax//10)) == 0 :
						if verbose ==1 : 
							print(' niter : %4i  distance : %1.2e --- evol(L) :%1.2e ---- evol(x) :%1.2e bnds=[%2.1f%% (%1.2e),%2.1f%%  (%1.2e)]'%(niter, np.linalg.norm(x-z),erreurl,erreurx,tab[0],list_bound[0],tab[1],list_bound[1]))
						#plt.clf()
						#plt.plot(z[:,0],z[:,1])
						#plt.plot(x[:,0],x[:,1])
						#plt.pause(0.01)
				if niter>nitermax+1 :
					if verbose > 0:
						print((' niter : %4i MAXITER REACHED distance : %1.2e --- Lagrange :%1.2e ---- x :%1.2e [%1.2e,%1.2e]'%(niter, np.linalg.norm(x-z),erreurl,erreurx,tab[0],tab[1])).center(128,'#'))
					return x,(False,niter,erreurl,erreurx,lagrange,tab)
				if (erreurx<tol) and (erreurl<tol) :
					if verbose > 0:
						print((' niter : %4i CONVERGENCE --- distance : %1.2e --- Lagrange :%1.2e ---- x :%1.2e [%1.2e,%1.2e]'%(niter, np.linalg.norm(x-z),erreurl,erreurx,tab[0],tab[1])).center(128,'#'))
					return x,(True,niter,erreurl,erreurx,lagrange,tab)
	def check_init(self,shape,As,projector,gamma,B,metric,verificators,solver) :
		for A in As :
			if not A.shape[1]==shape[0] :
				raise ValueError("All the matrices must have the shape[0] = %3i (found = %3i)"%(shape[0],A.shape[1]))
		if not len(As)==len(projector):
			raise ValueError("number of matrices (%3i) and of projectors must match (%3i)"%(len(As),len(projector)))
		if not len(As)==len(gamma):		
			raise ValueError("number of matrices (%3i) and of multipliers must match (%3i)"%(len(As),len(gamma)))
		if not len(As)==len(verificators):		
			raise ValueError("number of matrices (%3i) and of verificators must match (%3i)"%(len(As),len(verificators)))
		if not metric.shape==(shape[0],shape[0]) :
			raise ValueError('The metric does not have the correct size, found : '+str(metric.shape)+' must have '+str((shape[0],shape[0])))
		if not solver in ['factorize','spsolve','cholesky']:
			raise ValueError('The solver parameter is not acceptable, found : '+str(solver)+' must have in [factorize,spsolve,cholesky]')
	def check_solve(self,z,list_bound,b,xInit,lagrange) :
		if not z.shape==self.shape_x :
			raise ValueError('z does not have the correct size, found : '+str(z.shape)+' must have :'+str(self.shape_x))
		if not len(self.Alist)==len(list_bound) :
			raise ValueError("number of matrices (%3i) and of bounds must match (%3i)"%(len(self.Alist),len(list_bound)))
		if not xInit.shape==self.shape_x :
			raise ValueError('xInit does not have the correct size, found : '+str(xInit.shape)+' must have :'+str(self.shape_x))
		if not lagrange.shape==self.shape_y :
			raise ValueError('lagrange does not have the correct size, found : '+str(lagrange.shape)+' must have :'+str(self.shape_y))

def proj_inf_inf(x,g,bnd) :
	return np.sign(x)*np.minimum(np.abs(x),g*bnd)
	
def proj_inf_2(x,g,bnd) :
	norm=np.linalg.norm(x,axis=1)+1e-16
	I= norm>g*bnd+1e-16
	for i in range(x.shape[1]):
		x[I,i] = x[I,i]/norm[I]*g*bnd
	return x
def proj_2_2(x,g,bnd) :
	norm=np.linalg.norm(x)
	if norm>g*bnd :
		x/=norm
		x*=(g*bnd)
	return x

def verif_inf_inf(x,g) :
	return np.max(np.abs(x))/g
def verif_inf_2(x,g) :
	return np.max(np.linalg.norm(x,axis=1))/g
def verif_2_2(x,g) :
	return np.linalg.norm(x)/g	


def mat_deriv(n) :
	return sp.spdiags([-np.ones(n), np.ones(n)], [0,1], n-1, n)

def mat_deriv_periodic(n) :
	D=sp.spdiags([-np.ones(n), np.ones(n)], [0,1], n, n).tolil()
	D[-1,0]=1
	return D
	
def mat_deriv2(n) :
	return sp.spdiags([-np.ones(n), 2*np.ones(n),-np.ones(n)], [0,1,2], n-2, n)	

def mat_deriv2_periodic(n) :
	D=sp.spdiags([-np.ones(n), 2*np.ones(n),-np.ones(n)], [-1,0,1], n, n).tolil()	
	D[0,-1]=-1
	D[-1,0]=-1	
	return D
	
def mat_equal_indices(equal_indices,n) :
	unchanged_ind=np.setdiff1d(np.arange(n),equal_indices)
	col=np.arange(len(unchanged_ind)).tolist()
	col.extend(len(equal_indices)*[len(unchanged_ind)])
	row=unchanged_ind.tolist()
	row.extend(equal_indices)
	value=np.ones(len(unchanged_ind)).tolist()
	value.extend(len(equal_indices)*[1.])
	return sp.coo_matrix((value,(col,row)),shape=(n-len(equal_indices)+1,n))
	
