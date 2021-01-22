
# coding: utf-8

# # Aproximacion de raices de ecuaciones

# Aproximar soluciones de ecuaciones f(x)=0, primero con funciones de libreria: NumPy, SciPy(modulo optimize)
# 

# ### Ejercicio 1: 

# Aproximar las raices de P(x)=x^3-2x+1 usando una funcion de la libreria NumPy




import numpy as np

#introducimos los coeficientes del polinomio que deben estar en una lista.

p = [1, 0, -2, 1]

r = np.roots(p)

print(r)





p = [1, 0, 2]

r = np.roots(p)

print(r)


# ### Ejercicio 2

# Aproximar las raices del polinomio anterior utilizando la libreria SciPy

# In[7]:


get_ipython().magic('reset -f')





import scipy.optimize as sco

#Para el polinomio necesitamos definir una funcion.

def p(x) :
    return(x**3-2*x+1)

sco.root(p,-1.6)#Necesita una aproximacion incial para arrancar
    


# ### Ejercicio 3: 

# Programese una funcion llamada 'mibisecsion' que tenga como parametros de entrada el intervalo [a,b], donde hay un unico 0, una tolerancia de error, un numero maximo de iteraciones y debe devolver la aproximacion a la raiz.




def mibiseccion(f,a,b,tol,Nmax):
    
    """descripcion del metodo de biseccion"""
    
    import numpy as np
    
    cont, error = 0, tol+1 
    
    if f(a)*f(b)>=0:
        
        print('intervalo inadecuado')
        
        return
    
    while (cont < Nmax) and (error > tol) :
        
        aprox = (a+b)/2
        
        if f(aprox) == 0:
           
            return(aprox)
        
        elif f(a)*f(aprox)<0:
            
            b = aprox
            
            eror = np.fabs(b-a)
        else:
            
            a = aprox
            
        eror = np.fabs(b-a)
        
        cont = cont +1
        
        if (cont==Nmax):
            
            print('se ha alcanzado el numero maximo de iteracciones')
            
    return(aprox)
    





def f(x):
    return(x**2-1)
    





mibiseccion(f, 0, 6, 0.0001, 3)





get_ipython().magic('reset -f')





import MisMatematicas3





def f(x):
    return(x**2-1)
    





MisMatematicas3.mibiseccion(f, 0, 6, 0.0001, 3)


# # Derivacion e integracion numerica

# ### Ejercicio 1: 

# Aproximar la derivada en 0.5 de la funcion sen(2x)*cos(3x) usando la formula de diferencias centradas y tomando h=0.1 ... h=10^-7.
# 




import numpy as np





def f(x):
    
    return(np.sin(2*x)*np.cos(3*x))





h=0.1
aprox1= (f(0.5+h)-f(0.5-h))/(2*h)





H = []
for i in range(1,8):
    H.append(1/10**i)
    











H1 = np.linspace(-7,-1,7)
H1





10**H1





h = 10**H1





aprox2 =  (f(0.5+h)-f(0.5-h))/(2*h)





aprox2


# ### Ejercicio 3 

# Aproxima las siguientes funciones




get_ipython().magic('reset -f')





# Apartado A





import numpy as np





def f(x):
    return((np.sin((x**2)+1)**3)*np.cos(x))





aprox = (1/3)*(f(1)+4*f(2)+f(3))





aprox





#apartado B





h=1/6
aprox2 = (h/3)*(f(0)+4*f(1/6)+f(2/6))+(h/3)*(f(2/6)+4*f(3/6)+f(4/6))+(h/3)*(f(4/6)+4*f(5/6)+f(6/6))





aprox2


# comprobamos con funcion de libreria




import scipy.integrate as sc_in





sc_in.quad(f,1,3)


# ### Ejercicio 4 

# Interpola y dibuja la siguiente funcion




get_ipython().magic('reset -f')





import numpy as np
import matplotlib.pyplot as plt




def g(x):
    return(((x**2)+1)*np.cos(x+1))





#apartado A





def miinterpolador(x):
    return(0.54+(-2.74)*(x-0))
    





miinterpolador(1.3)





g(1.3)





x = np.linspace(0,6)
y1 = g(x)
y2 = miinterpolador(x)
plt.plot(x,y1,'b')
plt.plot(x,y2,'g')
plt.show()





get_ipython().magic('reset -f')

