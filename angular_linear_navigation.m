%% Choose clean or noisy data
data = sensors_clean(:,:);
%data = sensors_noisy(:,:);

truth = att_truth(2:4,:);
%% Initialisation 
% Initialise values
earthRate = 7.292115e-5;
initial_lat_rad = pos_llh(2,1);
R_0 = 6378137;
g = 9.79;
initial_height = pos_llh(4,1);

% Get data
accelerations = data(2:4,:);
roll = data(5,:);
pitch = data(6,:);
yaw = data(7,:);

len = length(data(1,:));
deltaT = data(1,2) - data(1,1);

%% Values for Angular Propagation
% Initial angles (phi, theta, psi)
phi = 0; theta = 0; psi = 5.8407;

% Make initial DCM
C = eulerToDCM(phi, theta, psi);

% Make initial q
q = zeros(4,len);
q(:,1) = eulerToQuart(phi, theta, psi);
    
% E will store angle data
E = zeros(3,len);

%% Values for Linear Propagation
% Initial Velocity Vectors
v_n = zeros(3,len);
v_n(:,1) = [vel_e(2,1) vel_e(3,1) vel_e(4,1)]';
v_ndot = zeros(3,len);

% Initialise Position Vectors
% Initial position (lat, lon, height)
initial_lon_rad = pos_llh(3,1);
pos = zeros(3,len);
pos(:,1) = [initial_lat_rad, initial_lon_rad, initial_height]';
lat = zeros(1,len);
lon = zeros(1,len);
h = zeros(1,len);

% Rates
% Transport rate, w_en (will be updated in loop)
w_en1 = v_n(2,1)/(R_0 + initial_height);
w_en2 = -v_n(1,1)/(R_0 + initial_height);
w_en3 = (-v_n(2,1) * tan(initial_lat_rad))/(R_0 + initial_height);
w_en = [w_en1, w_en2, w_en3]';

% Earth rate, w_ie (will be updated in loop)
w_ie = [earthRate*cos(initial_lat_rad), 0, -earthRate*sin(initial_lat_rad)]';

% Plumb bob gravity (will be updated in loop)
vec = [sin(2*initial_lat_rad), 0, 1+cos(2*initial_lat_rad)]';
gvec = [0,0,g]';
g_l = gvec - ((earthRate^2 * (R_0 + initial_height))/2) * vec; %

%% Propagation
for k = 1:len    
    % Angular
    w_ib = [roll(k) pitch(k) yaw(k)]';    
    ratec = C'*(w_ie + w_en);
    w_nb = w_ib - ratec;
    p = [0 w_nb']';
    
    Q = getQ(q(1,k), q(2,k), q(3,k), q(4,k));
    q_dot = 0.5 * Q * p;
    q(:,k+1) = q(:,k) + q_dot * deltaT; 
    
    % Normalise
    q(:,k+1) = q(:,k+1)/norm(q(:,k+1));
    
    C = quartToDCM(q(1,k+1), q(2,k+1), q(3,k+1), q(4,k+1));
    
    % Convert quaternion to Euler and store
    E(:,k) = quartToEuler(q(1,k+1), q(2,k+1), q(3,k+1), q(4,k+1));
    
    % Linear
    % Measured strap-down acceleration
    fb = [accelerations(1,k), accelerations(2,k), accelerations(3,k)]';
    fn = C * fb;
    
    % Acceleration and velocity
    v_ndot(:,k) = fn - cross((2*w_ie + w_en), v_n(:,k)) + g_l;
    v_n(:,k+1) = v_n(:,k) + deltaT*v_ndot(:,k);
    
    % Updating position
    pos_dot = [v_n(1,k+1)/(R_0 + pos(3,k)); 
              (v_n(2,k+1)*sec(pos(1,k)))/(R_0 + pos(3,k)); 
              -v_n(3,k+1)];
    pos(:,k+1) = pos(:,k) + deltaT*pos_dot;
    lat(k) = pos(1,k);
    lon(k) = pos(2,k);
    h(k) = pos(3,k);
    
    % Update rates and plumb bob gravity for next run
    % Transport rate, w_en
    w_en1 = v_n(2,k)/(R_0 + h(k));
    w_en2 = -v_n(1,k)/(R_0 + h(k));
    w_en3 = (-v_n(2,k) * tan(lat(k)))/(R_0 + h(k));
    w_en = [w_en1, w_en2, w_en3]';
    
    % Earth rate, w_ie
    w_ie = [earthRate*cos(lat(k)), 0, -earthRate*sin(lat(k))]';
    
    % Plumb bob gravity
    vec = [sin(2*lat(k)), 0, 1+cos(2*lat(k))]';
    gvec = [0,0,g]';
    g_l = gvec - ((earthRate^2 * (R_0 + h(k)))/2) * vec; %
end

%% Plotting Functions
% Plot 1: Angles
figure
% Plot Roll
subplot(3,1,1)
plot(att_truth(1,:), att_truth(2,:), 'm', att_truth(1,:), E(1,:), 'b');
title('Roll')
legend('Truth Data', 'Estimate', 'Location', 'northeast')
xlabel('Time') % x-axis label
ylabel('Radians') % y-axis label
set(gcf,'color','w');

% Plot Pitch
subplot(3,1,2)
plot(att_truth(1,:), att_truth(3,:), 'g', att_truth(1,:), E(2,:), 'b');
title('Pitch')
legend('Truth Data', 'Estimate', 'Location', 'southeast')
xlabel('Time') % x-axis label
ylabel('Radians') % y-axis label

% Plot Yaw
subplot(3,1,3)
plot(att_truth(1,:), att_truth(4,:), 'r', att_truth(1,:), E(3,:)+2*pi, 'b');
title('Yaw')
legend('Truth Data', 'Estimate', 'Location', 'northwest')
xlabel('Time') % x-axis label
ylabel('Radians') % y-axis label

% Plot 2: Estimation Error
figure
attitudeError = abs(truth-E).^2;
plot(att_truth(1,:), attitudeError(1,:), 'r', att_truth(1,:), attitudeError(2,:), 'g');
title('Estimation Error')
legend('Roll', 'Pitch', 'Location', 'northwest')
xlabel('Time') % x-axis label
ylabel('Estimation Error') % y-axis label
set(gcf,'color','w');

%% Plotting Linear Data
% Plot 3: Position
figure
plot3(pos_llh(2,:),pos_llh(3,:),pos_llh(4,:),pos(1,:),pos(2,:),pos(3,:))
title('Position')
legend('Truth Data', 'Estimate', 'Location', 'northeast')
xlabel('Latitude (radians)')
ylabel('Longitude (radians)')
zlabel('Altitude (metres)')
grid on
set(gcf,'color','w');

% For plotting against true data
time_vec = 1:120;
lat_truth = pos_llh(2,time_vec) * 180/pi;
lon_truth = pos_llh(3,time_vec) * 180/pi;
h_truth = pos_llh(4,time_vec);

figure
% Plot 4: Latitude vs Longitude
plot(lat_truth, lon_truth, 'g', lat_truth, lon(1,100*time_vec) * 180/pi, 'b');
title('Latitude vs Longitude');
legend('Truth Data', 'Estimate', 'Location', 'southwest');
xlabel('Latitude'); ylabel('Longitude');
set(gcf,'color','w');

% Plot 5: All Position Components
figure
% Plot Latitude
subplot(3,1,1)
plot(time_vec, lat_truth, 'g', time_vec, lat(1,100*time_vec) * 180/pi, 'b');
title('Latitude Estimate');
legend('Truth Data', 'Estimate', 'Location', 'southwest');
xlabel('Time'); ylabel('Degrees');

% Plot Longitude
subplot(3,1,2)
plot(time_vec, lon_truth, 'g', time_vec, lon(1,100*time_vec) * 180/pi, 'b');
title('Longitude Estimate')
legend('Truth Data', 'Estimate', 'Location', 'southwest');
xlabel('Time'); ylabel('Degrees');

% Plot Height
subplot(3,1,3)
plot(time_vec, h_truth, 'g', time_vec, h(1,100*time_vec), 'b');
title('Height Estimate')
legend('Truth Data', 'Estimate', 'Location', 'northwest');
xlabel('Time'); ylabel('Metres');
set(gcf,'color','w');

% Plot 6: Estimation Error
figure
positionError = abs(pos_llh(2:4,1:120)-pos(:,100*time_vec)).^2;
plot(time_vec, positionError(1,:), 'r', time_vec, positionError(2,:), 'g');
title('Estimation Error')
legend('Latitude', 'Longitude', 'Location', 'northwest')
xlabel('Time') % x-axis label
ylabel('Estimation Error') % y-axis label
set(gcf,'color','w');

% Plot 7: Estimation Error for Height
figure
plot(time_vec, positionError(3,:), 'r');
title('Estimation Error for Height')
xlabel('Time') % x-axis label
ylabel('Estimation Error') % y-axis label
set(gcf,'color','w');

%% Functions
function C = eulerToDCM(phi, theta, psi)
    C = zeros(3,3);
    C(1,:) = [cos(theta)*cos(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)  sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)];
    C(2,:) = [cos(theta)*sin(psi)  cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)];
    C(3,:) = [-sin(theta)          sin(phi)*cos(theta)                             cos(phi)*cos(theta)];
end

function x = eulerToQuart(phi, theta, psi)
    a = cos(psi/2)*cos(theta/2)*cos(phi/2)+sin(psi/2)*sin(theta/2)*sin(phi/2);
    b = cos(psi/2)*cos(theta/2)*sin(phi/2)-sin(psi/2)*sin(theta/2)*cos(phi/2);
    c = cos(psi/2)*sin(theta/2)*cos(phi/2)+sin(psi/2)*cos(theta/2)*sin(phi/2);
    d = sin(psi/2)*cos(theta/2)*cos(phi/2)-cos(psi/2)*sin(theta/2)*sin(phi/2);
    x = [a b c d]';
end

function Q = getQ(a,b,c,d)
    Q = zeros(4,4);
    Q(1,:) = [a -b -c -d];
    Q(2,:) = [b  a -d  c];
    Q(3,:) = [c  d  a -b];
    Q(4,:) = [d -c  b  a];
end

function C = quartToDCM(a, b, c, d)
    C = zeros(3,3);
    C(1,1) = a^2 + b^2 - c^2 - d^2;
    C(1,2) = 2*(b*c-a*d);
    C(1,3) = 2*(b*d+a*c);
    C(2,1) = 2*(b*c+a*d);
    C(2,2) = a^2-b^2+c^2-d^2;
    C(2,3) = 2*(c*d-a*b);
    C(3,1) = 2*(b*d-a*c);
    C(3,2) = 2*(c*d+a*b);
    C(3,3) = a^2-b^2-c^2+d^2;
end

function E = quartToEuler(a, b, c, d)
    tanYaw = ((2*(b*c+a*d))/(a^2+b^2-c^2-d^2));
    sinPitch = -2*(b*d-a*c);
    tanRoll = ((2*(c*d+a*b))/(a^2-b^2-c^2+d^2));
    
    yaw = atan2((2*(b*c+a*d)),(a^2+b^2-c^2-d^2));
    pitch = asin(sinPitch);
    roll = atan2((2*(c*d+a*b)),(a^2-b^2-c^2+d^2));
    
    E = [roll pitch yaw]';
end