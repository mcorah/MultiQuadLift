using PyPlot
using Sundials

tmax = 1e2;
dt = 1e-3;
d = 0.1
init = [1.0,1.0,3.0,0.0];

function f(t, yin, ydot)
  x = yin[1];
  y = yin[2];
  z = yin[3];
  dz = yin[4];
  ydot[1] = y;
  ydot[2] = -x;
  ydot[3] = dz;
  ydot[4] = -z - d*dz;
end

t = [0:dt:tmax];

pos = Sundials.cvode(f, init, t);

figure();
plot3D(pos[:,1],pos[:,2],pos[:,3]);
