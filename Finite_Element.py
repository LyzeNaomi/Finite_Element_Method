from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt
from pylab import *
import numpy as np

#Function 1: Entering the polygonal domain. I.e we move from the individual triangles to enter the entire domain, B=vector containing all nodes with their boundary specifications as third entry.
def InputPoly():
        Results=input('Enter the path to the file containing the nodes of triangles\n');Sheet=open(Results+'.txt','r');z=[];Px=[];Py=[];D=[];N=[];B=[]
        readFile=open(Results+'.txt','r');sepFile=readFile.read().split('\n'); sepFile.remove('');readFile.close()
        for plotFair in sepFile:
            xa=plotFair.split(' ')
            for xb in xa:
                xc=xb.split(' ');z.append((xc[0]))
        p=int(((len(z))/3))
        for n in range(0,p):
            Px.append(float(z[3*n]));Py.append(float(z[3*n+1]));B.append([float(z[3*n]),float(z[3*n+1]),(z[3*n+2])])
            if (z[3*n+2])=='D':
                D.append([(float(z[3*n])),(float(z[3*n+1]))])
                for a in D:
                    if D.count(a) > 1:
                        D.remove(a)
            elif (z[3*n+2])=='N':
                N.append([(float(z[3*n])),(float(z[3*n+1]))])
                for b in N:
                    if N.count(b) > 1:
                        N.remove(b)
            for ap in B:
                if B.count(ap)>1:
                    B.remove(ap)
        Px.append(Px[0]);Py.append(Py[0]);x=[];y=[];s=[];m=int(input('Enter the number of refinements desired please. '))
        lamda=float(input('Enter the value of the constant lamda. ')); f1=input('Enter the source term please\nf(x,y)= ')
        U=input('Please enter the true solution\n U(x,y)=');#plt.plot(Px,Py);plt.show()
        while (Px!=[] and Py!=[]):
            x1=float(Px[0]);Px.remove(Px[0]);y1=float(Py[0]);Py.remove(Py[0]);x.append(x1);y.append(y1)
            if len(x)==len(y)==3:
                s.append([x,y]);x=[];y=[]
        print('Input OK !!!!')
        return(Px,Py,s,m,lamda,f1,D,N,B,U)
Px,Py,s,m,lamda,f1,D,N,B,U=InputPoly()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Function 2: Collecting all edges with Dirichlet or Neumann function defined on each edge.
def AllEdges():
    Results = input('Enter the path to the file containing the boundary edges and functions. \n');Sheet = open(Results + '.txt', 'r')
    z = [] ; Aedge=[];readFile = open(Results + '.txt', 'r');sepFile = readFile.read().split('\n');sepFile.remove('');readFile.close()
    for plotFair in sepFile:
        xa = plotFair.split(' ')
        for xb in xa:
            xc = xb.split(' ');z.append((xc[0]))
    p=int((len(z)/5))
    for n in range (0,p):
        Aedge.append([[float(z[5*n]),float(z[5*n+1])],[float(z[5*n+2]),float (z[5*n+3])],(z[5*n+4])])
    print('All edges='+str(Aedge))
    print('All edges collected finished !!!!')
    return Aedge
Aedge=AllEdges()
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 2: Collecting all boundary edges.
def Edges(s):
    E1=[];E=[]
    for e in s:
        for i in range(0,2):
            for j in range (i+1,3):
                xn=[e[0][i],e[1][i]] ; yn=[e[0][j],e[1][j]] ; E1.append([xn,yn])
    for plus in E1:
        if E1.count(plus)==1:
            E.append(plus)
    print('\n')
    print('Boundary edges='+str(E))
    return E
E=Edges(s)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 3: Creating Dirichlet and Neumann Edges De and Ne respectively.
def BEdges(B,E):
    De=[]; Ne=[]
    for i in range (0,(len(B))-1):
        for j in range (i+1, len(B)):
            if 'D' in B[i] and 'D' in B[j]:
                x=[[B[i][0],B[i][1]],[B[j][0],B[j][1]]] ; y=[[B[j][0],B[j][1]],[B[i][0],B[i][1]]]
                if x in E or y in E:
                    De.append(x)
                    for more in De:
                        if De.count(more)>1:
                            De.remove(more)
            if 'N' in B[i] and 'N' in B[j] or ('D' in B[i] and 'N' in B[j]) or ('N' in B[i] and 'D' in B[j]):
                xN = [[B[i][0], B[i][1]], [B[j][0], B[j][1]]];yN = [[B[j][0], B[j][1]], [B[i][0], B[i][1]]]
                if xN in E or yN in E:
                    Ne.append(xN)
                    for excess in Ne:
                        if Ne.count(excess)>1:
                            Ne.remove(excess)
    for a in De:
        if a[0][0]>a[1][0] or a[0][1]>a[1][1]:
            De[De.index(a)]=[a[1],a[0]]
    for b in Ne:
        if b[0][0]>b[1][0] or b[0][1]>b[1][1]:
            Ne[Ne.index(b)]=[b[1],b[0]]
    print('Creating Dirichlet and Neumann Edges finished !!!!')
    return De, Ne
De, Ne=BEdges(B,E)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 4:  Refinement.
def Ref( s,m):
    p=0;v=[];Mx=[];My=[]
    while p<m:
        Init1=[0];Init2=[0]
        for a in s:
            for j in range (0,2):
                for k in range (j+1,3):
                    x2=(a[0][k]+a[0][j])/2;y2=(a[1][k]+a[1][j])/2;Mx.append(x2);My.append(y2)
            Mx=[Mx[0],Mx[2],Mx[1]];My=[My[0],My[2],My[1]];s1=[a[0][0],Mx[0],Mx[2]],[a[1][0],My[0],My[2]]
            s2=[Mx[0], Mx[1], Mx[2]],[My[0],My[1],My[2]];s3=[Mx[0],a[0][1],Mx[1]],[My[0],a[1][1],My[1]]
            s4=[Mx[2],Mx[1],a[0][2]],[My[2],My[1],a[1][2]];si=[s1,s2,s3,s4];Mx=[];My=[]
            Vec1=Init1+s1[0]+s2[0]+s3[0]+s4[0];Vec2=Init2+s1[1]+s2[1]+s3[1]+s4[1];Init1=Vec1;Init2=Vec2
            for b in si:
                v.append(b)
            si=[]
        Vec1.remove(0);Vec2.remove(0);p+=1;s=v;v=[];Vec1.append(Vec1[0]);Vec2.append(Vec2[0])
    print('Refinement finished !!!!')
    return (s,Vec1,Vec2)
(s,Vec1,Vec2)=Ref( s,m)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 5: Collecting all nodes in the domain, including boundary nodes and boundary coordinates of each node.
def Node(s):
    Nodes=[] ; All=[] ; BCoor=[] ; BNode=[]
    for o in s:
        if ([o[0][0], o[1][0]]) not in Nodes:
            Nodes.append([o[0][0], o[1][0]])
        if ([o[0][1], o[1][1]]) not in Nodes:
            Nodes.append([o[0][1], o[1][1]])
        if ([o[0][2], o[1][2]]) not in Nodes:
            Nodes.append([o[0][2], o[1][2]])
        All.append([o[0][0], o[1][0]]);All.append([o[0][1], o[1][1]]);All.append([o[0][2], o[1][2]])
    for f1 in Nodes:
        if All.count(f1) < 6:
            BNode.append(Nodes.index(f1))
    for f2 in BNode:
        BCoor.append(Nodes[f2])
    print('Number of Nodes='+str(len(Nodes)))
    return Nodes, All, BNode, BCoor
Nodes,All,BNode,BCoor=Node(s)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 6: Classifying boundary nodes in the boundary they belong.
def Boundary(D,N,De,BCoor):
    for bcoor in BCoor:
        for de in De:
            if (de[0][0]<= bcoor[0] <= de[1][0]) and (de[0][1]<= bcoor[1] <= de[1][1]):
                D.append(bcoor)
                for a in D:
                    if D.count(a)>1:
                        D.remove(a)
    for nt in BCoor:
        if nt not in D:
            N.append(nt)
            for b in N:
                if N.count(b)>1:
                    N.remove(b)
    D1=D ; N1=N
    print('Classification of nodes complete !!!!')
    return D1,N1
D1,N1= Boundary(D,N,De,BCoor)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 7: Creating Neumann edges after discretization that will be used for integration.
def NeumannE(s,D1,N1):
    Ne = []
    for o in s:
        f1=[o[0][0],o[1][0]] ; f2=[o[0][1],o[1][1]] ; f3=[o[0][2],o[1][2]]
        if (f3 in N1 and f2 in N1) or (f2 in N1 and f3 in N1) or (f3 in D1 and f2 in N1) or (f2 in D1 and f3 in N1):
            if f3[0]<=f2[0] and f3[1]<=f2[1]:
                Ne.append([f3,f2])
            elif f2[0]<=f3[0] and f2[1]<=f3[1]:
                Ne.append([f2, f3])
        if (f3 in N1 and f1 in N1) or (f1 in N1 and f3 in N1) or (f3 in D1 and f1 in N1) or (f1 in D1 and f3 in N1):
            if f3[0]<=f1[0] and f3[1]<=f1[1]:
                Ne.append([f3,f1])
            elif f1[0]<=f3[0] and f1[1]<=f3[1]:
                Ne.append([f1,f3])
        if (f2 in N1 and f1 in N1) or (f in N1 and f2 in N1) or (f2 in D1 and f1 in N1) or (f1 in D1 and f2 in N1):
            if f1[0]<=f2[0] and f1[1]<=f2[1]:
                Ne.append([f1, f2])
            elif f2[0]<=f1[0] and f2[1]<=f1[1]:
                Ne.append([f2, f1])
    for w in Ne:
        if Ne.count(w)>1:
            Ne[:]=(value for value in Ne if value!=w)
    Ne2=Ne
    print('Neumann edges finished !!!!')
    return Ne2
Ne2=NeumannE(s,D1,N1)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 8: Collecting global node numbering for each triangle in the domain, NNum= vector containing global nodes numbering of each triangle.
def NNumbre(s,Nodes):
    NNum = []
    for o in s:
        NNum.append([Nodes.index([o[0][0],o[1][0]]),Nodes.index([o[0][1],o[1][1]]),Nodes.index([o[0][2],o[1][2]])])
    print('Nodes numbering finished !!!!')
    return NNum
NNum=NNumbre(s, Nodes)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 9: Collecting interior nodes both coordinates and global numbering.
def INodes(s,All,Nodes):
    INode=[];INCoor=[]
    for f in Nodes:
        if All.count(f)==6 or All.count(f)>6:
            INode.append(Nodes.index(f))
    for f1 in INode:
        INCoor.append(Nodes[f1])
    return INode,INCoor
INode,INCoor=INodes(s,All,Nodes)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 10: Plotting discretization
def Plotting(Vec1,Vec2):
    Q=int(input('Please enter 1 if the image of the refinment wants to be viewed or any other number otherwise  '))
    if Q==1:
        plt.plot(Vec1,Vec2);plt.savefig("Figure.png", transparent = True);plt.show()
    return()
#Plotting(Vec1,Vec2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 11: Building the Elements Stiffness matrices
def ElementM(s):
    r = 1;K = []
    for t in s:
        Js = [[t[0][1] - t[0][0], t[0][2] - t[0][0]], [t[1][1] - t[1][0], t[1][2] - t[1][0]]];Js = reshape(Js, [len(Js), len(Js)]);Js = np.linalg.det(Js)
        Kr = zeros([3, 3], float)
        Kr[0][0] = (1 / (2 * Js)) * (((t[1][1] - t[1][2]) ** 2) + ((t[0][2] - t[0][1]) ** 2));Kr[1][1] = (1 / (2 * Js)) * (((t[1][2] - t[1][0]) ** 2) + ((t[0][2] - t[0][0]) ** 2))
        Kr[2][2] = (1 / (2 * Js)) * (((t[1][1] - t[1][0]) ** 2) + ((t[0][1] - t[0][0]) ** 2)); Kr[0][1] = Kr[1][0] = (1 / (2 * Js)) * (((t[1][2] - t[1][0]) * (t[1][1] - t[1][2])) + ((t[0][2] - t[0][0]) * (t[0][1] - t[0][2])))
        Kr[0][2] = Kr[2][0] = (1 / (2 * Js)) * (((t[1][0] - t[1][1]) * (t[1][1] - t[1][2])) + ((t[0][0] - t[0][1]) * (t[0][1] - t[0][2])))
        Kr[1][2] = Kr[2][1] = (1 / (2 * Js)) * (((t[1][2] - t[1][0]) * (t[1][0] - t[1][1])) + ((t[0][2] - t[0][0]) * (t[0][0] - t[0][1])));K.append(Kr)#;print('K'+str(r)+'='+str(Kr));r+=1
    print('Building Element Stiffness finished !!!!')
    return (Js,K)
Js,K=ElementM(s)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 12: Building Global Stiffness Matrix
def GlobalM(K,Nodes,NNum):
    A1 = zeros([len(Nodes), len(Nodes)], float)
    for a in NNum:
        r=NNum.index(a)
        for i in range (0,3):
            for i1 in range (0,3):
                j=a.index (a[i]) ; j1=a.index(a[i1])
                A1[a[i]][a[i1]]=A1[a[i]][a[i1]] + K[r][j][j1]
    print('Building global Stiffness finished !!!!')
    return A1
A1= GlobalM(K,Nodes,NNum)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 13: Computing the entries in the mass matrix: \lamda u , where \lamda \in R.
def MassMatrix(lamda,s):
    k=[]
    for t in s:
        Js = [[t[0][1] - t[0][0], t[0][2] - t[0][0]], [t[1][1] - t[1][0], t[1][2] - t[1][0]]];Js = reshape(Js, [len(Js), len(Js)]);Js = np.linalg.det(Js)
        ka=[1/12,1/24,1/24,1/24,1/12,1/24,1/24,1/24,1/12];ka=np.reshape(ka,[3,3]);ka=lamda*Js*ka
        k.append(ka)
    return k
k=MassMatrix(lamda,s)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 14: Assembling the global mass matrix.
def GlobalMass(k,Nodes,NNum):
    A2 = zeros([len(Nodes), len(Nodes)], float)
    for a in NNum:
        r=NNum.index(a)
        for i in range (0,3):
            for i1 in range (0,3):
                j=a.index (a[i]) ; j1=a.index(a[i1])
                A2[a[i]][a[i1]]=A2[a[i]][a[i1]] + k[r][j][j1]
    return A2
A2= GlobalMass(k,Nodes,NNum)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 15: Adding matrices A1 and A2 to form A:
def Matrix(A1,A2):
    A=A1+A2
    return A
A=Matrix(A1,A2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 16: Calculating the element load vectors using dblquad.
from sympy import *
from sympy.utilities import lambdify
from numpy import *
from scipy.integrate import dblquad
def ElementL(s,f1):
    eta, phi = symbols('eta phi');v1=[1-phi-eta, phi,eta];F=[]
    for t in s:
        fi=[]
        Js = [[t[0][1] - t[0][0], t[0][2] - t[0][0]], [t[1][1] - t[1][0], t[1][2] - t[1][0]]];Js = reshape(Js, [len(Js), len(Js)]);Js = np.linalg.det(Js)
        x = (t[0][1] - t[0][0]) * phi + (t[0][2] - t[0][0]) * eta + t[0][0]; y = (t[1][1] - t[1][0]) * phi + (t[1][2] - t[1][0]) * eta + t[1][0];f=eval(f1)
        for h in v1:
            produit=h*f*Js ; prdt=lambdify((phi,eta),(h*f*Js)) ; g=dblquad ( prdt, 0,1, lambda phi:0, lambda  phi: 1-phi) ;  fi.append(g[0])
        F.append(fi)
    print('Building Element load finished !!!!')
    return F
F = ElementL(s,f1)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 17: Assembling global load vector.
def GlobalL(F,Nodes,NNum):
    b=np.zeros([len(Nodes),1],float)
    for a in NNum:
        r=NNum.index(a)
        for i in range (0,3):
                j=a.index (a[i])
                b[a[i]]=b[a[i]] + F[r][j]
    print('Building first load finished  !!!!')
    return b
b=GlobalL(F,Nodes,NNum)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Function 18: Computing Neumann function at Neumann boundary nodes ; gn=Neumann function evaluated at Neumann nodes.
from scipy.integrate import quad
def Neumann (Ne2,Aedge):
    xi=symbols('xi'); v=[1-xi,xi];FN=[]
    for a in Ne2:
        for a1 in Aedge:
                if a1[0][0]<=a[0][0] and a1[0][1]<=a[0][1] and a[1][0]<=a1[1][0] and a[1][1]<=a1[1][1]:
                        Fi=[]
                        l=sqrt(((a[0][0]-a[1][0])**2)+((a[0][1]-a[1][1]))**2) ; x=a[0][0]+xi*(a[1][0]-a[0][0]) ; y=a[0][1]+xi*(a[1][1]-a[0][1]) ; gn= eval(a1[2])
                        for h in v:
                            product=lambdify((xi),(h*gn*l)) ; g=quad(product, 0,1) ; Fi.append(g[0])
                        FN.append(Fi)
    print('Neumann Contribution finished !!!!')
    return FN
FN = Neumann(Ne2,Aedge)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 19: Assembling Neumann boundary condition; First, we classify the Neumann edges using nodes position in the domain.
def NeuNumb (Ne2,Nodes,N1,b,FN):
    Nindex=[]
    for a in Ne2:
            Nindex.append([Nodes.index(a[0]),Nodes.index(a[1])])
    print(Nindex)
    for i in Nindex:
        r=Nindex.index(i)
        for j in range (0,2):
            j1=i.index(i[j])
            b[i[j]]=b[i[j]] + FN[r][j1]
    return Nindex,b
Nindex,b=NeuNumb(Ne2,Nodes,N1,b,FN)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Function 20: Computing the Dirichlet boundary values and correcting the right hand side.
def Dirichlet(D1,Nodes,b,A,De,Aedge):
    Dnumb=[] ; Many=[]
    for f in D1:
        Dnumb.append(Nodes.index(f))
    for a in Dnumb:
        h = Nodes[a]
        for a2 in Aedge:
            if [a2[0],a2[1]] in De:
                if a2[0][0]<=h[0]<=a2[1][0] and  a2[0][1]<=h[1]<=a2[1][1] :
                    x = h[0];y = h[1] ; g = eval(a2[2])
                    for i in range (0,len(Nodes)):
                        b[i]=b[i] - ((A[i][a])*g)
    return b, Dnumb
b, Dnumb=Dirichlet(D1,Nodes,b,A,De,Aedge)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 21: Deleting Dirichlet nodes on both the stiffness matrix and the load vector
def ReduceMatrix(A,b,Dnumb):
    while Dnumb!=[]:
        New=[]
        A = np.delete (A,(Dnumb[0]),axis=0)
        A = np.delete(A, (Dnumb[0]), axis=1)
        b = np.delete(b, (Dnumb[0]), axis=0)
        Dnumb.remove(Dnumb[0])
        for a in Dnumb:
            b1=a-1;New.append(b1)
        Dnumb = New
    return A, b
A,b=ReduceMatrix(A,b,Dnumb)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Function 22: Copying the matrix outside.
def Writing(A,b):
    New=[]
    Sheet=open('Stiffness.txt','w')
    File=open('Load.txt','w')
    n1=(len(A[0]));r=0
    while r<n1:
        Sheet.write(" ".join(str(e) for e in A[r]))
        r+=1
        Sheet.write('\n')
    for c1 in b:
        for c in c1:
            New.append(c)
    File.write(" ".join(str(e) for e in New))
    File.write('\n')
    return b
b=Writing(A,b)


#Evaluate true function at internal nodes
def TrueSolution(Nodes,D1,U):
    UV=[] ; Sheet=open ('Solution obtained at non-Dirichlet nodes.txt','w')
    for a in Nodes:
        if a not in D1:
            x=a[0];y=a[1]
            u=eval(U)
            UV.append(u)
    n2=len(UV)
    #UV=reshape(UV,[n2,1])
    Sheet.write(" ".join(str(e) for e in UV))
    Sheet.write('\n')
    return UV
UV=TrueSolution(Nodes,D1 ,U)
