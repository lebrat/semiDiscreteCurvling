import sys
import numpy as np
import rt 
import scipy.sparse as sps
from scipy.sparse.linalg import LinearOperator,cg

def calcul_masse(img,isQ1) :
	if not isQ1 :
		return np.sum(img.ravel())/(img.shape[0]*img.shape[1])
	else :
		Nx=img.shape[0]
		Ny=img.shape[1]
		mass=np.sum(img[1:-1,1:-1])
		mass+=0.5*(np.sum(img[0,1:-1])+np.sum(img[-1,1:-1])+np.sum(img[1:-1,0])+np.sum(img[1:-1,-1]))
		mass+=0.25*(img[0,0]+img[-1,0]+img[0,-1]+img[-1,-1])
		mass=mass/((Nx-2)*(Ny-2)+0.5*2*(Ny+Nx-4)+1.)
		return mass

class OT() :
	def __init__(self,u,box,isQ1=False):
		self.T = rt.Triangulation2()
		self.T.setImage(u.shape[0],u.shape[1],isQ1,u.ravel(),box[0],box[1],box[2],box[3])
	def compute(self,x,psi,w,hessian=None,parallel=False):
		self.T.setCoordinates(x[:,0],x[:,1])
		n=x.shape[0]
		bar=np.zeros(2*n)
		costlocal=np.zeros(n)
		mass=np.zeros(n)
		self.T.computeLaguerre(psi)
		self.nnzHess = self.T.computeTesselation()
		val=np.zeros(self.nnzHess)
		row=np.zeros(self.nnzHess,dtype=int)
		col=np.zeros(self.nnzHess,dtype=int)
		if parallel :
			self.T.IntegrationParall(bar,costlocal,mass,val,col,row)
		else:
			self.T.Integration(bar,costlocal,mass,val,col,row)	
		if hessian is None :
			Hess=None
		elif hessian=='psi' :
			Hess=sps.csr_matrix((val,(row,col)),(psi.shape[0],psi.shape[0]))	
		bar=bar.reshape((n,2))
		costlocal=-costlocal+2*x[:,0]*bar[:,0]+2*x[:,1]*bar[:,1]-(x[:,0]**2+x[:,1]**2)*mass
		I = np.where(np.abs(mass)>1e-20)
		bar[I,0]=bar[I,0]/mass[I]
		bar[I,1]=bar[I,1]/mass[I]
		cf = np.sum(costlocal)
		cf += np.dot(psi,mass)
		cf -= np.dot(psi,w)
		Grad = mass - w		
		return mass,bar,cf,Grad,Hess,costlocal
	def Tesselation(self,x,psi):
		"""
		Return a dict with key = Inside and val is a np.array where [0] are the Xcoord and [1] Ycoord
		"""
		self.T.setCoordinates(x[:,0],x[:,1])		
		self.T.computeLaguerre(psi)
		Laguerre = self.T.computeTess()
		result = dict([(int(val[0]),np.reshape(np.array(val[1:]),(2,(len(val)-1)/2),order='F'))for val in Laguerre])
		return result


	def Adjacency(self,x,psi):
		self.T.setCoordinates(x[:,0],x[:,1])		
		self.T.computeLaguerre(psi)
		nitroGlycerine = self.T.computeAdjacency()
		result = dict([(int(val[0]),np.array(val[1:]).tolist())for val in nitroGlycerine])
		return result

	def compute_local_cost_matrix(self,x,psi):
		
		self.T.setCoordinates(x[:,0],x[:,1])
		n=x.shape[0]
		momentX2=np.zeros(n)
		momentY2=np.zeros(n)
		momentXY=np.zeros(n)
		self.T.computeLaguerre(psi)
		self.T.ComputeLocalMatrix(momentX2,momentY2,momentXY)
		return momentX2,momentY2,momentXY
class PostProcess() :
	""" USAGE IS POSTPROCESS(ot,x,psi,w)
	# INPUT 
	#	x : either a array(n,2) or array(2*n) vector, coordinates of Dirac Masses
	#   psi : array(n) , power of the power diagramm
	#   w :  array(n) , weights of Dirac masses
	# OUTPUT
	# computes for a Laguerre tesselation several info
	# self.mass : array(n) of the masses of the Laguerre Tesselation
	# self.bar  : array(n,2) of the barycenters of the Laguerre Tesselation
	# self.G : value of the transport for the choice of psi
	# self.x : array(2*n), copy of x
	# self.Grad_P_G : array(n), derivative of self.G wrto psi
	# self.Grad_X_G : array(2*n), derivative of self.G wrto x
	# self.Hess_PP_G : sparse matrix (n,n), Hessian of self.G wrto psi
	# self.D_X_Popt(d) : Linear function, return array(n) if d is array(2*n)
	#			returns the derivative of the optimal psi wrto x in direction d if
	#			psi is assumed to be optimal for x 
	# self.Hess_XX_Gopt(d) : Linear function, return u=array(2*n) if d=array(2*n)
	#			returns u, the multplication of d by the Hessian of G_opt wrto X assuming that
	#			psi is optimal
	# self.Hess_G : Matrix (3*n,3*n)
	#			returns the Hessian of G wrto (psi,X) 
	# self.Hess_XX_G(d) : Matrix (2*n,2*n)
	#			returns the Hessian of G wrto (X) 
	"""
	def __init__(self,ot,x,psi,w,parallel=False):
		self.nbPoints=psi.shape[0]
		#import pdb;pdb.set_trace()
		self.x=np.copy(x).ravel(order='F')
		self.x2=np.reshape(x,(self.nbPoints,2),order='F')
		self.mass,self.bar,self.G,self.Grad_P_G,self.Hess_PP_G,self.G_local=ot.compute(self.x2,psi,w,hessian='psi')
		self.Grad_X_G=-((self.x2 - self.bar)*(2*self.mass).repeat(2).reshape((self.nbPoints,2))).ravel(order='F')
		val=np.zeros(6*ot.nnzHess)
		row=np.zeros(6*ot.nnzHess,dtype=int)
		col=np.zeros(6*ot.nnzHess,dtype=int)
		if parallel :
			ot.T.HessianComputeParall(val,col,row)
		else :
			ot.T.HessianCompute2(val,col,row)
		self.Hess_WX_G = sps.csr_matrix((-val[:2*ot.nnzHess],(row[:2*ot.nnzHess],col[:2*ot.nnzHess])),(self.nbPoints,2*self.nbPoints)) 
		self.Hess_XX_G = sps.csr_matrix((-val[2*ot.nnzHess:],(row[2*ot.nnzHess:],col[2*ot.nnzHess:])),(2*self.nbPoints,2*self.nbPoints))
		self.Hess_XW_G = self.Hess_WX_G.transpose(copy=True)

	def Hess_XX_Gopt(self,direction) :
		return self.Hess_XX_G.dot(direction)+self.Hess_XW_G.dot(self.D_X_Popt(direction))

	def D_X_Popt(self,direction) :
		u=-self.Hess_WX_G.dot(direction)
		normU=np.linalg.norm(u)
		if normU > 1.e-16 :
			return cg(self.Hess_PP_G,u/normU,tol=1e-12)[0]*normU
		else :
			return np.zeros(u.shape)
	def Hess_G(self,d):
		r=np.zeros(3*self.nbPoints)
		r[:self.nbPoints]=self.Hess_PP_G.dot(d[:self.nbPoints])+self.Hess_WX_G.dot(d[self.nbPoints:])
		r[self.nbPoints:]=self.Hess_XW_G.dot(d[:self.nbPoints])+self.Hess_XX_G.dot(d[self.nbPoints:])
		return r
