

Files provided:

1.- att_truth.mat:

truth attitude in euler angles. vector <4x12041>. sampled at t=0.01seg
format (time;roll;pitch;yaw)

2.- vel_e.mat
gps velocity coverted to (NED). vector <4x121>, sampled at t=1seg 


3.-pos_llh.mat
4.-pos_ecef.mat
gps position in llh and ecef respectively
. sampled at 1 seg. <4x121> vector
fotmat: (time;lat;long;alt)

5.-sensors_clean.mat
6.-sensors_noisy.mat
sensors measuremenst without and with noise, respectively. <7x12041>.
sampled at t=0.01 seg
format: (time; ax,;ay ;az; wx; wy; wz )


