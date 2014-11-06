tmax = 1e0;
dt = 1e-2;
init = [0.2,0.2,0.1,0.1,0.0,0.0,0.0,0.0,0.0,0.0];

global M
global Kp
global Kdp
global Ka
global Kda
global g
global A
global B
I = 1.0;
M = 1.0;
g = 9.8;
wn_rp = 9;
wn_z = 5;
wn_xy = 2;
Kp=M * diag([wn_xy^2,wn_xy^2,wn_z^2]);
Kdp=M * diag([wn_xy,wn_xy,wn_z]);
Ka=I * diag([wn_rp^2,wn_rp^2]);
Kda=I * diag([wn_rp,wn_rp]);

A = [
     zeros(3,3) eye(3) zeros(3,4);
     zeros(2,6) [0 g; -g 0] zeros(2,2);
     zeros(1,10);
     zeros(2,8) eye(2);
     zeros(2,10)
    ];
B = [
     zeros(5,3);
     1/M 0 0;
     zeros(2,3);
     0 1/I 0;
     0 0 1/I;
    ];



t = [0:dt:tmax];

[T,pos] = ode45(@f, t, init);

figure();
%disp(pos);
plot3(pos(:,1),pos(:,2),pos(:,3));
