function ydot = f(t, yin)
  global M
  global Kp
  global Kdp
  global Ka
  global Kda
  global g
  global A
  global B
  p = yin(1:3,1);
  dp = yin(4:6,1);
  a = yin(7:8,1);
  da = yin(9:10,1);
  acc_des = (Kp*(-p) + Kdp*(-dp));
  att_des = 1/g .* [0 -1 0; 1 0 0] * acc_des;
  u = [
       (acc_des(3));
       [Ka * (att_des-a) + Kda * (-da)]
      ];
  ydot = A*yin + B*u;
  disp(ydot');
end
