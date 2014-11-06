using PyPlot
using Sundials
import Base.show

tmax = 1e0
dt = 1e-6
init = [0.2,0.2,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0]

I = 1.0
M = 1.0
g = 9.8
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

#show(A)

function f(t, yin, ydot)
  p = yin[1:3]
  dp = yin[4:6]
  a = yin[7:8]
  da = yin[9:10]
  acc_des = (Kp*(-p) + Kdp*(-dp));
  att_des = 1/g .* [0 -1.0 0; 1.0 0 0] * acc_des
  u = [
       (acc_des[3]);
       Ka * (att_des-a) + Kda * (-da)
      ]
  #show(u)
  #print('\n')
  ydot = A*yin + B*u
  #show(ydot)
  print('\n')
end;

t = [0:dt:tmax]

pos = Sundials.cvode(f, init, t)

figure()
#show(pos)
plot3D(pos[:,1],pos[:,2],pos[:,3])
