import math
import matplotlib.pyplot as plt

class Resevoir:
    def __init__(self, h, hmin, h_outlet, n, k, L, Qin, A, g=9.82):
        self.n = n
        self.h = h #level above outlet
        self.hmin = hmin #min level above the outlet
        self.h_outlet = h_outlet
        self.k = k
        self.L = L
        self.Qin = Qin
        self.A = A #outlet cross sectional area
        self.g = g #gravity
   
    def Crosssection(self, h):
        Width = (h/self.k)**(1/self.n)
        return 2*self.L*Width
   
    def Q_out(self, h):
        #outflow
        if h - self.hmin < 0:
            return 0
        w = math.sqrt((2*self.g*(h - self.h_outlet)/(1 - (self.A/self.Crosssection(h))**2)))
        return w*self.A
   
    def dh_dt(self, t, h):
        #change in water height
        S = self.Crosssection(h)
        if S == 0:
            return 0
        return (self.Qin-self.Q_out(h))/S
    
    def Analytical_n1(self): # Note that this only works for n=1
        C0 = math.sqrt(2/self.g)*self.L/(self.A*self.k)
        C1 = self.Qin / (self.A * math.sqrt(2*self.g))
        
        return (
            2*C0*(
                ((self.h-self.h_outlet)**1.5 - (self.hmin-self.h_outlet)**1.5)/3
                    + C1*(self.h - self.hmin)/2
                        + (C1**2 + self.h_outlet) * (math.sqrt(self.h-self.h_outlet) - math.sqrt(self.hmin-self.h_outlet))
                            + (C1**3 + C1*self.h_outlet) * math.log((math.sqrt(self.h-self.h_outlet) - C1)/(math.sqrt(self.hmin-self.h_outlet) - C1))
            )
        )
    
def Euler_step(reservoir, t, h, dt):
    dh = reservoir.dh_dt(t, h)
    return h + dt * dh

def Heun_step(reservoir, t, h, dt):
    k1 = reservoir.dh_dt(t, h)
    k2 = reservoir.dh_dt(t+dt, h+dt*k1)
    return h + dt*(k1 + k2)*0.5

def RK4_step(resevoir, t, h, dt):
    k1 = resevoir.dh_dt(t, h)
    k2 = resevoir.dh_dt(t+dt/2, h + dt*k1/2)
    k3 = resevoir.dh_dt(t+dt/2, h + dt*k2/2)
    k4 = resevoir.dh_dt(t+dt, h + dt*k3)
    return h + dt*(k1 + 2*k2 + 2*k3 + k4)/6


def main():
    res = Resevoir(
    h=100,
    hmin=80,
    h_outlet = 70,
    n=1,
    k=1,
    L=1000,
    Qin=20,
    A=10
    )
    print(res.Analytical_n1())

    t = 0
    dt = 0.01
    h = res.h
    T = []
    H = []
    for _ in range(2500000):
        T.append(t)
        H.append(h)
        h = Euler_step(res, t, h, dt)
        t += dt
    plt.plot(T, H)
    plt.show()
    
if __name__ == "__main__":
    main()

