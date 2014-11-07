using PyPlot
using SDE
import Base.show

tmax = 1e1
dt = 1e-4
R = 1
init = [0.2,0.2,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0]

I = 1.0
M = 1.0
g = 9.8
mu = 0.1
wn_rp = 9.0
wn_z = 5.0
wn_xy = 2.0
Kp=M * diagm([wn_xy^2,wn_xy^2,wn_z^2])
Kdp=M * diagm([wn_xy,wn_xy,wn_z])
Ka=I * diagm([wn_rp^2,wn_rp^2])
Kda=I * diagm([wn_rp,wn_rp])

A = [
     zeros(3,3) eye(3) zeros(3,4);
     zeros(2,6) [0 g; -g 0] zeros(2,2);
     zeros(1,10);
     zeros(2,8) eye(2);
     zeros(2,10)
    ]

B = [
     zeros(5,3);
     1/M 0 0;
     zeros(2,3);
     0 1/I 0;
     0 0 1/I;
    ]

C = [
     zeros(3,3);
     mu * eye(3,3);
     zeros(4,3);
    ]


function f(x)
  p = x[1:3]
  dp = x[4:6]
  a = x[7:8]
  da = x[9:10]

  acc_des = (Kp*(-p) + Kdp*(-dp));
  att_des = 1/g .* [0 -1.0 0; 1.0 0 0] * acc_des

  U = [
       (acc_des[3]);
       Ka * (att_des-a) + Kda * (-da)
      ]

  A*x + B*U
end;

t = [0:dt:tmax]

X,T = SDE.em(f, x->C, init, tmax, dt, R)

figure()
temp = X'
plot3D(temp[:,1],temp[:,2],temp[:,3])
