import sympy as sp

t = sp.Symbol('t')
x, y, z = sp.symbols('x y z')
E, G = sp.symbols('E G')
K = sp.Symbol('K')
L = sp.Symbol('L')

# Strain Values

def strain_x(u):
    return sp.diff(u,x)

def strain_y(v):
    return sp.diff(v,y)

def strain_z(w):
    return sp.diff(w,z)
    
def strain_yz(v,w):
    return 1/2*(sp.diff(v,z)+sp.diff(w,y))

def strain_xz(u,w):
    return 1/2*(sp.diff(u,z)+sp.diff(w,x))

def strain_xy(u,v):
    return 1/2*(sp.diff(u,y)+sp.diff(v,x))

# Stress Values

def stress_x(u):
    return E*strain_x(u)

def stress_y(v):
    return E*strain_y(v)

def stress_z(w):
    return E*strain_z(w)

def stress_yz(v,w):
    return 2*G*K*strain_yz(v,w)

def stress_xz(u,w):
    return 2*G*K*strain_xz(u,w)

def stress_xy(u,v):
    return 2*G*K*strain_xy(u,v)

# Integration by parts

def IntByParts(u,dv):
    v = sp.integrate(dv,x)
    I = ((u*v).subs(x,L)-(u*v).subs(x,0)) - sp.integrate(v*sp.diff(u,x),(x,0,L))
    return I

# Equation of Motion from Langrangian Density

w = sp.Function('w')(x,t)
w_t = sp.diff(w,t)
w_x = sp.diff(w,x)
w_xx= sp.diff(w,x,2)

def GovEq(L,q):
    eq_Lhs = sp.diff(L,w) - sp.diff((sp.diff(L,w_t)),t) - sp.diff((sp.diff(L,w_x)),x) + sp.diff((sp.diff(L,w_xx)),x,2) + q
    return eq_Lhs