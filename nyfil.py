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
        return self.L*Width*2
   
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
        return (self.Qin - self.Q_out(h))/S
    
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
    time_step_for_equilibrium = 0.01
    iterations_for_equilibrium = 20000000
    list_of_time_steps = [0.1*2**i for i in range(5)]
    iterations = [int(2000000/(2**i)) for i in range(5)]
    res = Resevoir(
    h=100,
    hmin=80,
    h_outlet = 70,
    n=1,
    k=1,
    L=1000,
    Qin=20,
    A=1
    )

    t = 0
    dt = time_step_for_equilibrium
    h = res.h
    T = []
    H = []
    for _ in range(iterations_for_equilibrium):
        T.append(t)
        H.append(h)
        h = Euler_step(res, t, h, dt)
        t += dt
    equilibrium = H[-1]


    list_of_errors = []
    for i in range(len(list_of_time_steps)):
        t = 0
        dt = list_of_time_steps[i]
        h = res.h

        T = []
        H = []

        for _ in range(iterations[i]):
            T.append(t)
            H.append(h)
            h = Euler_step(res, t, h, dt)
            t += dt

        list_of_errors.append(abs(H[-1]-equilibrium))

        plt.plot(T, H)
        plt.show()
    starting_error = list_of_errors[0]
    theoretical_development =  [starting_error*(2**i) for i in range(5)]    

    print(list_of_errors)
    plt.plot(list_of_time_steps, list_of_errors, '-o', label = 'Numerical Errors')
    plt.plot(list_of_time_steps, theoretical_development, '--', label = 'Theoretical Errors')
    plt.title("Order of Accuracy")
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main()

