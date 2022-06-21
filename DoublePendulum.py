from vpython import *
#Web VPython 3.2

#for energy graphs 
TotalE = gcurve(color=color.blue, label="Total Energy")
Kinetic = gcurve(color=color.red, label="Kinetic Energy")
Potential = gcurve(color=color.green, label = "Potential Energy")

#initial conditions for first bob 
m1 = 5
l1 = 1
theta1 = 1.61
thetaDot1 = 3

x1 = l1*sin(theta1)
y1 = -l1*cos(theta1)

#initial conditions for second bob 
m2 = 3
l2 = 1 
theta2 =  1.54
thetaDot2 = 0.1

x2 = l1 * sin(theta1) + l2 * sin(theta2)
y2 = -l1 * cos(theta1) - l2 * cos(theta2) 

g = 9.8 

#making the pendulum 
pivot = sphere(pos = vector(0, l1, 0), radius = l1 * 0.0005)
bob1 = sphere(pos = pivot.pos + vector(x1, y1, 0), radius = 0.05*l1, color = color.yellow)  
arm1 = cylinder(pos = pivot.pos, axis = bob1.pos, color = color.white, length = 0.008, radius = 0.008)

bob2 = sphere(pos = pivot.pos + vector(x2, y2, 0), radius = 0.05*l2, color = color.purple, make_trail = True)
arm2 = cylinder(pos = bob1.pos, axis = bob2.pos - bob1.pos, color = color.white, length = 0.008, radius = 0.008)

#animate
dt = 0.0001 
time = 0.0

while (time <= 30): 
    
    rate(10000) 
    
    #constants obtained from solving the Euler-Lagrange equation for the given system 
    a = -(m1 + m2)* g * l1 * sin(theta1) - m2* l1 *l2 * thetaDot2**2 * sin(theta1-theta2)
    b =  (m1 + m2) * l1**2
    c =  m2 * l1 * l2 * cos(theta1-theta2)
    d = -m2 * g * l2 * sin(theta2) + m2 * l1 * l2 * thetaDot1**2 * sin(theta1-theta2)
    e =  m2 * l2**2
    f =  m2 * l1 * l2 * cos(theta1-theta2)
    
    #equations for angular acceleration 
    tdd2 = (d - a * f / b)/(e - c * f / b)
    tdd1 = a / b - c * (tdd2) / b
    
    #equations for angularvelocity 
    thetaDot2 += tdd2 * dt 
    thetaDot1 += tdd1 * dt
    
    #equations for angles 
    theta1 += thetaDot1 * dt
    theta2 += thetaDot2 * dt 
    
    time += dt 
    
    #updating positions 
    x1 = l1*sin(theta1)
    y1 = -l1*cos(theta1)
    
    x2 = l1 * sin(theta1) + l2 * sin(theta2)
    y2 = -l1 * cos(theta1) - l2 * cos(theta2) 
    
    bob1.pos = pivot.pos + vector(l1*sin(theta1),-l1*cos(theta1),0)
    bob2.pos = pivot.pos + vector(x2, y2, 0)
    arm1.axis = bob1.pos - pivot.pos
    arm2.pos = bob1.pos
    arm2.axis = bob2.pos - bob1.pos 
    
    T = 0.5*m1*l1**2*thetaDot1**2 + 0.5*m2*(l2**2*thetaDot2**2 + 2 * l1 * l2 * thetaDot1 * thetaDot2 * cos(theta1-theta2) + l1**2 * thetaDot1**2)
    U = -(m1+m2) * g* l1* cos(theta1) - m2*g*l2*cos(theta2)
    
    TotalE.plot(time, T + U)
    Kinetic.plot(time, T)
    Potential.plot(time, U)

print("done")