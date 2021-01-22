
# coding: utf-8

# # Sistemas Lineales. Metodos directos

# Ax=b; A es un array nxn   b es un array nx1



import numpy as np
import scipy.linalg as scla # El comando para resolver sistemas lineales es solve
A = np.array([[1,2,3],[2,4,1],[-1,-1,2]])
b = np.array([[1],[5],[3]])
scla.solve(A,b)


# ## Factorizacion LU




solucion = scla.lu(A)# 1 arry es L, 2 array es U





L = solucion[0]
L





L, U, P = scla.lu(A)





L





U





P





#Factorizacion de cholesky, el comando es cholesky
scla.cholesky(A) #la matriz no es definida positiva





def Newton_Rapshon(f,f1,x0,tol, Nmax):
    """"""
    import numpy as np
    
    for i in range(1,Nmax+1):
        x1 = x0-(f(x0)/f1(x0))
        error= np.abs(x0-x1)
        if error<tol:
            return(x1)
        else:
            x0 = x1
        if i >= Nmax:
            print('numero maximo alcanzado')
            return(x1)


# ## Metodo de Jacobi- metodo de Gauss-Seidel




#Vamos a programar antes la iteracion de punto fijo


# Programese una funcion llamada "mipuntofijo" que tenga como parametros de entrada g, x0, tol(tolerancia max), Nmax(nº maximo de iteraciones) y devulva el punto fijo o un mensaje de error




def mipuntofijo (g, x0, tol, Nmax):
    """cosas"""
    import numpy as np
    
    for i in range(1,Nmax+1):
        x1 = g(x0)
        error= np.fabs(x0-x1)
        if error<tol:
            return(x1)
        else:
            x0 = x1
        if i >= Nmax:
            print('numero maximo alcanzado')
            return(x1)
            




def g(x):
    return(2**(-x))






mipuntofijo(g,1,0.00001,5)


# ### Programacion metodo de Jacobi




get_ipython().magic('reset -f')


# Programaremos el metodo de jacobi utilizando la forma matricial: x^{m+1}=B*x^m+D^-1*b siendo B = D{-1^}*(-(L+U)) y C=D^-1*b

# Programese una funcion llamada "mijacobi" que tenga como parametros de entradas una matriz A, un vector b, una tolerancia de error tol, y un numero maximo de iteranciones Nmax. Debe imprimir por pantalla el mensaje de error adecuado si se excede Nmax o no se alanza la tolerancia de error debe devolver el numero de iteraciones utilizado, la aproximacion al error y ademas imprimirlo por pantalla con un mensaje adecuado.




import numpy as np
A = np.array([[1.2,2.8],[2.3,4.6]])
np.shape(A) #numero de filas y columnas
np.transpose(A)#transpuesta de la matriz
a = np.array([[1],[2]])#vector en colmuna
np.diagonal(A)#Nos proporciona la diagonal de A
b = np.diag(a)
l = np.tril(A)#triangulas inferior
u = np.triu(A)#triangular superior
producto = A.dot(a)#matriz por vector






print(b)





print(l)





print(u)





print(producto)





import numpy.linalg as npla





npla.inv(A)#calcula la inversa de A
npla.norm(A)#calcula la norma infinito
npla.norm(A,2)#norma 2





def mijacobi(A,b,tol,Nmax,x0):
    """"""
    import numpy as np
    import numpy.linalg as npla
    
    error = tol+1
    i = 0
    D = np.diag(np.diagonal(A))#saca la diagonal y hace una matriz diagonal
    L = np.tril(A)-D
    U = np.triu(A)-D
    B = np.dot(npla.inv(D),(-(L+U)))
    print(B)
    C = np.dot(npla.inv(D),b)
    
    
    while (error>=tol and i<=Nmax):
        i=i+1
        x = np.dot(B,x0)+C
        
        error = npla.norm(x-x0)
        
        if(i==Nmax):
            print('se ha alcanzado el numero maximo de iteracciones')
            return(x,i)
        if error < tol:
            print('el error es menor que la tolerancia.')
            return(x,i)
        else:
            x0 = x





import numpy as np
A = np.array([[10,2,1],[1,5,1],[2,3,10]])
b = np.array([[7],[-8],[6]])
I = np.array([[1,0,0],[0,1,0],[0,0,1]])
x0 = np.array([[0.7],[-1.6],[0.6]])





mijacobi(A,b,0.1,4,x0)





get_ipython().magic('reset -f')


# ### Metodo de la secante 

# x(n+1)=x(n)-[(x(n)-x(n-1))*f(x(n))]/[f(x(n))-f(x(n-1))]




def metodo_secante(F,x0,x1,tol,Nmax):
    """"""
    import numpy as np
    
    error = tol+1
    i = 0
    
    while (error>=tol and i<=Nmax):
        i = i+1
        x2 = x1-((x1-x0)*F(x1))/(F(x1)-F(x0))
        error= np.fabs(x0-x1)
        if(i==Nmax):
            print('se ha alcanzado el numero maximo de iteracciones')
            return(x2)
        if error < tol:
            print('el error es menor que la tolerancia.')
            return(x2,i)
        else:
            x0 = x1
            x1 = x2





import numpy as np
def F(x):
    return(np.cos(5*x)-3*x)





metodo_secante(F,1/6,1/3,0.00000001,100)


# # Simpson




def simp_comp(f,a,b,m):
    """"""
    import numpy as np
    x=0
    i=1
    h=(b-a)/(2*m)
    while i<=m:
        x=x+(f(a+h*(2*i-2))+4*f(a+h*(2*i-1))+f(a+h*(2*i)))
        i=i+1
        
    return(x*(h/3))







a=0
b=1
h=0.125
x=0
m=(b-a)/(2*h)
def f(x):
    return(x/((x+1)*(x+2)))





simp_comp(f,a,b,m)





def simpson(f,a,b):
    """"""
    import numpy as np
    h = (b-a)/2
    
    x = (h/3)*(f(a)+4*f((a+b)/2)+f(b))
    return(x)



# ### Prgramacion metodo Runge-Kuta 

# Metodo RK de orden 4 clasico (apuntes)




def miRK4(f,t0,y0,tf,h):
    """bla bla bla"""
    import numpy as np
    
    t1 = t0+h
    i = 0
    
    while(t1 <=tf):
        k1 = f(t0,y0)
        k2 = f(t0+0.5*h,y0+0.5*h*k1)
        k3 = f(t0+0.5*h,y0+0.5*h*k2)
        k4 = f(t0+h,y0+h*k3)
        y1 = y0+h*((1/6)*k1+(1/3)*k2+(1/3)*k3+(1/6)*k4)
        
        i = i +4
        t0 = t1
        t1 = t1+h
        y0 = y1
    if (t0)<tf:
        
        h = tf-(t0)
        k1 = f(t0,y0)
        k2 = f(t0+0.5*h,y0+0.5*h*k1)
        k3 = f(t0+0.5*h,y0+0.5*h*k2)
        k4 = f(t0+h,y0+h*k3)
        y1 = y0+h*((1/6)*k1+(1/3)*k2+(1/3)*k3+(1/6)*k4)
        
        y0 = y1
        i = i+4
            
    print('T final {0:-23.15e}, aproximacion {1:-23.15e}, numero de evaluaciones de funcion {2:23d}'.format(t0+h,y1,i)) #imprimir con formato   
    return(tf,y0,i)
        


# y'=-y; y(0)=1




def f(t,y):
    return(-y)





def f(t,y):
    return(-y)
t0 = 0
y0 = 1
tf = 1
h = 0.1
sol = miRK4(f,t0,y0,tf,h)





import numpy as np
np.exp(-sol[0])





np.exp(-sol[0])-sol[1]


# # Funciones libreria para aproximar soluciones de ecuaciones diferenciales

# El paquete scipy.integrate contienes funciones para resolver EDO: odeint, Isoda, vode, dopri5, dop853

# y'=-y; y(0)=1




get_ipython().magic('reset -f')





import scipy.integrate as scin
import numpy as np
import matplotlib.pyplot as plt




def f(x,y):
    return(-y)

def sol(x):
    return(np.exp(-x))





y00 = 0, 1
x0 = 0
y0 = 1
xf = 1





x = np.linspace(x0,xf,100)#discretizacion, dividir en intervalos
y = scin.odeint(f,y0,x)
y1 = sol(x)
plt.plot(x,y,'r')
plt.plot(x,y1,'b')
plt.show()


# Polinomio interpolador de lagrange forma de newtoncon las diferencias divididas
