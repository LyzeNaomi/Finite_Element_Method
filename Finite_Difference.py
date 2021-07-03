from numpy import *
import scipy
import scipy.linalg
import numpy as np



def Input():
    print('PLEASE CHOOSE HOW YOU WANT TO COLLECT THE INPUT INFORMATION')
    q=int(input('\033[1;34mPRESS \033[1;34m\n\033[1;38m1:TO COLLECT THE INPUT MANUALLY\033[1;38m\n\n2:TO COLLECT THE INPUT FROM A FILE\n'))
    if q==1:

        a=float(input('Enter the domain of the function \n a='));b=float(input('\n b='))
        while(b==a or b<a):
            print('\033[1;31m!!!Please a should be less than b!!!!\033[1;31m')
            a=float(input('\033[1;38ma=\033[1;38m'));b=float(input('\n b='))
        c=float(input('\n c='));d=float(input('\n d='))
        while(d==c or d<c):
            print('\033[1;31m!!!Please c should be less than d!!!!\033[1;31m')
            c=float(input('\033[1;38mc=\033[1;38m'));d=float(input('\n d='))
        #print('The domain is \n (a,b)x(c,d)='+str((a,b))+'x'+str((c,d)))
        n=float(input('Enter the number of uniform partitions for the interval[a,b]=' +'\n n='))
        while(n==0 or n<0):
            print('\033[1;31m!!!n should be a natural number different from 0!!!\033[1;31m')
            n=float(input('\033[1;38mn=\033[1;38m'))
        m=float(input('Enter the number of uniform partitions for the interval[c,d]=' +'\n m='))
        while(m==0 or m<0):
            print('\033[1;31m!!!m should be a natural number different from 0!!!\033[1;31m')
            m=float(input('\033[1;38mn=\033[1;38m'))
        f1=input('Enter the source term please\n f(x,y)=')
        q=input('Enter the Dirichlet Boundary Condition please \n U(x,y)=')
        U1=input('Enter the true solution of U please \n U(x,y)=')
        n2=int(input('Enter the number of times you want to do the refinement\n n='))
        return(a,b,c,d,f1,q,U1,m,n,n2)
    else:
        Results=input('Enter the path to the file \n')
        Sheet=open(Results+'.txt','r')
        x=[];y=[]
        readFile=open(Results+'.txt','r')
        sepFile=readFile.read().split('\n')
        #print(sepFile)
        readFile.close()
        for plotFair in sepFile:
            xa=plotFair.split(',')
            for xb in xa:
                xc=xb.split(',')
                x.append((xc[0]))
        x.remove('')
        print(x)
        for a in x:
            y.append((a))
        print(y)
        n=len(y)
        a=float(y[0]);b=float(y[1]);c=float(y[2]);d=float(y[3]);m=int(y[4]);n=int(y[5]);f1=y[6];q=float(y[7]);U1=y[8];n2=int(y[9])
        return (m,n,a,b,c,d,f1,q,U1,n2)

(m,n,a,b,c,d,f1,q,U1,n2)=Input()

#Partitioning of omega(lexicographic)
def Hk(a,b,c,d,m,n):
    global h,k
    h=round(((b-a)/n),6)
    print('h='+str(h))
    k=round((((d-c)/m)),6)
    print('h='+str(h)+'\n'+'k='+str(k))
    return(h,k)
h,k=Hk(a,b,c,d,m,n)

#Computation of the system matrix
#i.Initialization and modification of interior nodes of l.
def Matrixl(a,c,m,n,h,k):
    x2=[]
    y2=[]
    for j in range (0,m+1):
        for i in range (0,n+1):
            y1=round((c+j*k),6)
            y2.append(y1)
            a1=round((a+i*h),6)
            x2.append (a1)

    i=n+1;j=m+1;l=zeros([(i*j),(i*j)],float);f=0;g=1;s=((n+1)*(m+1));w=h*h;w1=k*k;w2=2*(w+w1);w3=w*w1
    for j in range (0,s):
        for i in range(f,g):
            l[i][i]=w2/w3
            f=f+1;g=g+1

#ii.Modifying the diagonal elements of l which are different from zero

    for j in range (0,s):
        for i in range (j+1,s):
            if (round(x2[i]-x2[j],6)==h and y2[i]-y2[j]==0):
                l[i][j]=l[j][i]=(-1)*(h**(-2))
            elif (x2[i]-x2[j]==0 and round(y2[i]-y2[j],6)==k):
                l[i][j]=l[j][i]=(-1)*(k**(-2))
    return(l,x2,y2)
(l,x2,y2)=Matrixl(a,c,m,n,h,k)
#D.Defining the right hand side vector and evaluating the function at the interior nodes.
def Definf(m,n,x2,y2,f1):
    f=[]
    for j in range(1,m):
        for i in range(1,n):
            x=x2[i+j*(n+1)]
            y=y2[i+j*(n+1)]
            (eval(str(f1)))
            f.append(float(repr(eval(str(f1)))))
            #f.append(f1)
    return f
n3=len(f)
f=Definf(m,n,x2,y2,f1)

#Calculate true solution of u at the interior nodes.
def TrueSolnU(m,n,x2,y2,U1):
    U=[]
    for j in range (1,m):
        for i in range (1,n):
            x=x2[i+j*(n+1)]
            y=y2[i+j*(n+1)]
            (eval(str(U1)))
            U.append(float(eval(str(U1))))
    return U
U=TrueSolnU(m,n,x2,y2,U1)
U1=reshape(U,[n3,1])
#E.Define the Dirichlet B.C
#(i.) Initialization.
def Defing(m,n):
    g=[]
    for j in range (0,m+1):
        for i in range (0,n+1):
            g1=i+j*(n+1)
            g1=0
            g.append(g1)
    return g
g=Defing(m,n)
#(ii.) Modify the values on the boundary omega h bar.Function phy is x+y
def ModifBdry(m,n,x2,y2,q,g):
    for i in range (0,n+1):
        x=x2[i];y=y2[i]
        g[i]=float(repr(eval(str(q))))
    for j in range (1,m+1):
        i=j+(j+1)*n
        x=x2[i];y=y2[i]
        g[i]=float(repr(eval(str(q))))
    for i in range (0,n):
        j=i+(n+1)*m
        x=x2[j];y=y2[j]
        g[j]=float(repr(eval(str(q))))
    for j in range(1,n):
        i=j+(j*n)
        x=x2[i];y=y2[i]
        g[i]=float(repr(eval(str(q))))
    return (g)
g=ModifBdry(m,n,x2,y2,q,g)
#G.Correct the R.H.S vector to take care of the D.C on the boundary segments.
#(i.) Boundary segment 1:
def CorrectBdry(m,n,x2,y2,q,f,l):
    h1=[];z=0
    for j in range (1,m):
        for i in range (1,n):
            x=x2[i];y=y2[i];u=i+j*(n+1);q1=float(eval(str(q)));f[z]=f[z]-(q1*l[u][i]);h1.append(f[z])
            z+=1
    z=0;h2=[]
    #(ii.) Boundary segment 2:
    for j in range (1,m):
        for i in range (1,n):
            u=i+j*(n+1);u1=n+j*(n+1);x=x2[u1];y=y2[u1];q1=float(eval(str(q)));h1[z]=f[z]-(q1*l[u][u1]);h2.append(h1[z])
            z+=1
    z=0;h3=[]
    #(iii.) Boundary segment 3:
    for j in range(1,m):
        for i in range(1,n):
            u=i+j*(n+1);u1=j*(n+1);x=x2[u1];y=y2[u1];q1=float(eval(str(q)));h2[z]=f[z]-(q1*l[u][u1]);h3.append(h2[z])
            z+=1
    '''
     #(iv.) Boundary segment 4:
    for j in range (1,m):
        for i in range (1,n):
           u=i+j*(n+1);u1=j*(n+1);x=x2[u1];y=y2[u1];q1=float(eval(str(q)));h1[z]=f[z]-(q1*l[u][u1]);h2.append(h1[z])
           z+=1
    '''
    return (h3)
h3=CorrectBdry(m,n,x2,y2,q,f,l)
#deleting rows
def ReduceLToA(m,n,l):
    s=((n+1)*(m+1))
    for j in range (0,m+2):
        j=0
        l=np.delete (l,(j),axis=0)
        l=np.delete(l,(j),axis=1)
    k1=s-1-(m+2)
    for j in range (0,m+2):
        l=np.delete(l,(k1),axis=0)
        l=np.delete(l,(k1),axis=1)
        k1=k1-1
    for k1 in range (1,n-1):
        for j in range (0,2):
            l=np.delete(l,(k1*(n-1)),axis=0)
            l=np.delete(l,(k1*(n-1)),axis=1)
    return l
l=ReduceLToA(m,n,l)
def Writing(l,h3,U):
    Sheet=open('Matrix.txt','w')
    File=open('Matrix1.txt','w')
    Folio=open('Matrix2.txt','w')
    A=l;b=h3;s = set(A[0]);n1=(len(A[0]));r=0
    while r<n1:
        Sheet.write(" ".join(str(e) for e in A[r]))
        r+=1
        Sheet.write('\n')
    File.write(" ".join(str(e) for e in h3))
    File.write('\n')
    Folio.write(" ".join(str(e) for e in U))
    Folio.write('\n')

    return
Writing(l,h3,U)

#import Linear
