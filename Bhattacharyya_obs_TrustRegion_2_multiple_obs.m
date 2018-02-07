
clc;
clear all;
close all;


%% Drone is supposed to reach (dest_x, dest_y, dest_z) from (0,0,0) in
%% 'n' timesteps with each timestep of duration del_t
%% Version 2 (for higher values of n)
n = 40;      % Number of timesteps
del_t = 1.0; % Time duration for which some command is executed.
tot_time = n*del_t; %% Total time of flight

%% defining drone uncertainty
drone_sig = diag([0.5^2 0.5^2 0.5^2]);

%% Defining obstacle trajectory
%% For now we will assume that obstacle is moving along a straing line. it's starting from (0,0,0) and It's trajectory will be determined based on the destination of the drone

%% Drone Initial and Destination variables
start_x = 0;
start_y = 0;
start_z = 0;
dest_x = 12;
dest_y = 12;
dest_z = 12;

%%%%%%%%%%%%%% Obstacle stuff starts here!!!! %%%%%%%%%%%%%%%%%
%% Here, the drone is moving from (0,0,0) to (dest_x, dest_y, dest_z) in n timesteps
%% The obstacle will reach from start_obs to end_obs in n timesteps 
Ob_start_1=[0, dest_y/2, dest_z];
Ob_end_1=[dest_x, dest_y/2, 0];

Ob_start_2=[-dest_z/3, 0, dest_x/3];                              %%Meeting at 1/3
Ob_end_2=[dest_x + 2*dest_z/3, dest_y, dest_z - 2*dest_x/3];

Ob_start_3=[2*dest_x/3, 0, 2*dest_z/3];                           %%Meeting at 2/3
Ob_end_3=[2*dest_x/3, dest_y, 2*dest_z/3];

%% Calculating velocity of the obstacle in all the directions
vel_obs_1 = (Ob_end_1 - Ob_start_1)/tot_time;
vel_obs_2 = (Ob_end_2 - Ob_start_2)/tot_time;
vel_obs_3 = (Ob_end_3 - Ob_start_3)/tot_time;

%% Defining Obstacle path location at each timestep
for i = 1:n
  Ox_1(i) = Ob_start_1(1) + del_t*i*vel_obs_1(1);
  Oy_1(i) = Ob_start_1(2) + del_t*i*vel_obs_1(2);
  Oz_1(i) = Ob_start_1(3) + del_t*i*vel_obs_1(3);
  
  Ox_2(i) = Ob_start_2(1) + del_t*i*vel_obs_2(1);
  Oy_2(i) = Ob_start_2(2) + del_t*i*vel_obs_2(2);
  Oz_2(i) = Ob_start_2(3) + del_t*i*vel_obs_2(3);
  
  Ox_3(i) = Ob_start_3(1) + del_t*i*vel_obs_3(1);
  Oy_3(i) = Ob_start_3(2) + del_t*i*vel_obs_3(2);
  Oz_3(i) = Ob_start_3(3) + del_t*i*vel_obs_3(3);
end

%% The confidence has to be between [0, 1]
confidence = 0.1;

%% Defining obstalce uncertainty
obs_sig = diag([1.0^2 1.0^2 1.0^2]);

%% Getting mean sigma of the obstacle and the robot
mean_sig = (obs_sig + drone_sig)/2;

%% getting Mahalanobis threshold based on confidence interval
mahalanobisThresh = chi2inv(confidence, 2);

%% Radius of avoidance
R = 2.0; %% This is a Bhattacharyya radius

%% Calculating Eigen_values of the covariance matrix.
%% The eigen_values will give us lengths of semi-axis lengths
sig_obs_x = obs_sig(1,1);
sig_obs_y = obs_sig(2,2);
sig_obs_z = obs_sig(3,3);
eigen_values = eig(obs_sig);

%% Getting eigen_values for X,Y,Z direction
%% The eigen values returned from eig() are sorted in decreasing order
if(sig_obs_x <= sig_obs_y && sig_obs_x <= sig_obs_z)
	eig_x = eigen_values(1,1);

	if(sig_obs_y < sig_obs_z)
		eig_y = eigen_values(2,1);
		eig_z = eigen_values(3,1);
	else
		eig_y = eigen_values(3,1);
		eig_z = eigen_values(2,1);
	end
end

if(sig_obs_y <= sig_obs_z && sig_obs_y <= sig_obs_x)
	eig_y = eigen_values(1,1);

	if(sig_obs_x < sig_obs_z)
		eig_x = eigen_values(2,1);
		eig_z = eigen_values(3,1);
	else
		eig_x = eigen_values(3,1);
		eig_z = eigen_values(2,1);
	end
end

if(sig_obs_z <= sig_obs_x && sig_obs_z <= sig_obs_y)
	eig_z = eigen_values(1,1);

	if(sig_obs_x < sig_obs_y)
		eig_x = eigen_values(2,1);
		eig_y = eigen_values(3,1);
	else
		eig_x = eigen_values(2,1);
		eig_y = eigen_values(3,1);
	end
end


%% Getting lengths of semi major axis based on confidence values

%% Radius of ellipse along X-direction
rad_X = sqrt(mahalanobisThresh*eig_x); 

%% Radius of ellipse along Y-direction
rad_Y = sqrt(mahalanobisThresh*eig_y); 

%% Radius of ellipse along Z-direction
rad_Z = sqrt(mahalanobisThresh*eig_z); 


%%%%%%%%%Obstacle stuff ends here %%%%%%%%%%%%%%%

%% Difference between subsequent velocities
del_Vx = 0.25;
del_Vy = 0.25;
del_Vz = 0.25;

delta_bar = rand;   %%No between 0-1
%delta_bar = 0.2;
trust_del_x = 0 + (delta_bar-0).*rand(n,1); %%No between 0-delta_bar   
trust_del_y = 0 + (delta_bar-0).*rand(n,1);
trust_del_z = 0 + (delta_bar-0).*rand(n,1);
%nu = 0 + (0.25-0).*rand;
guess_vx = 0.2*ones(n,1);
guess_vy = 0.15*ones(n,1);
guess_vz = 0.1*ones(n,1);

l = 0;
iter = 0;
%no_of_iter = 3;

while l >= 0
    
    cvx_begin
    
    variables Vx(n) Vy(n) Vz(n) %% Velocities
    %variables Px(n) Py(n) Pz(n) %% Drone locations at each timestep
    
    for i = 1:n
        Px(i)  =  start_x + sum(Vx(1:i))*del_t;
        Py(i)  =  start_y + sum(Vy(1:i))*del_t;
        Pz(i)  =  start_z + sum(Vz(1:i))*del_t;
    end
    
    %% Cost function
    %%minimise((sum(Vx)*del_t - dest_x)^2 + (sum(Vy)*del_t - dest_y)^2 + (sum(Vz)*del_t - dest_z)^2)
    minimise((Px(n) - dest_x)^2 + (Py(n) - dest_y)^2 + (Pz(n) - dest_z)^2)
    
    %%Obstacle avoindace constrains
    %% original constraint is as follows
    %% (Px(i) - Ox(i))^2/(sig_x^2) + (Py(i) - Oy(i))^2/(sig_y^2) + (Pz(i) - Oz(i))^2/(sig_z^2) >= R^2
    %% Since the above constraint is quadratic in nature, we will linearize it around point x(i), y(i), z(i)
    %% f(Px(i),Py(i),Pz(i)) = (Px(i) - Ox(i))^2 + (Py(i) - Oy(i))^2 + (Pz(i) - Oz(i))^2 linearized around x(i), y(i), z(i) will look like following
    %% (x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2(x(i) - Ox(i))*(Px(i) -  Ox(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i))
    %% Constraint will be (x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2(x(i) - Ox(i))*(Px(i) -  Ox(i)) + 2(y(i) - Oy(i))*(Py(i) -  Oy(i)) + 2(z(i) - Oz(i))*(Pz(i) -  Oz(i)) > R^2
    vx = guess_vx;
    vy = guess_vy;
    vz = guess_vz;
    
    for i = 1:n
        %% Linearized Bhattacharyya constraint in terms of velocity
        actual_lin_1(i) = ((start_x + sum(vx(1:i))*del_t - Ox_1(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_1(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_1(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_1(i) = ((start_x + sum(vx(1:i))*del_t - Ox_1(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_1(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_1(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_1(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_1(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_1(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_1(i) >= R
        
        actual_lin_2(i) = ((start_x + sum(vx(1:i))*del_t - Ox_2(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_2(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_2(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_2(i) = ((start_x + sum(vx(1:i))*del_t - Ox_2(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_2(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_2(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_2(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_2(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_2(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_2(i) >= R
        
        actual_lin_3(i) = ((start_x + sum(vx(1:i))*del_t - Ox_3(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_3(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_3(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_3(i) = ((start_x + sum(vx(1:i))*del_t - Ox_3(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_3(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_3(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_3(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_3(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_3(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_3(i) >= R
    end
    %%%%%%%% Constraints %%%%%%%%
    Vx >= 0;
    Vy >= 0;
    Vz >= 0;
    
    %% The velocity during initialization shouldn't be much higher
    0 <= Vx(1) <= 1;
    0 <= Vy(1) <= 1;
    0 <= Vz(1) <= 1;
    
    (vx - trust_del_x) <= Vx <= (vx + trust_del_x);
    (vy - trust_del_y) <= Vy <= (vy + trust_del_y);
    (vz - trust_del_z) <= Vz <= (vz + trust_del_z);
    
    %% Ensuring that subsequent X vel commands do not have difference more than del_Vx
    0 <= Vx(2:n) - Vx(1:n-1) <= del_Vx;
    
    %% Ensuring that subsequent Y vel commands do not have difference more than del_Vy
    0 <= Vy(2:n) - Vy(1:n-1) <= del_Vy;
    
    %% Ensuring that subsequent Z vel commands do not have difference more than del_Vz
    0 <= Vz(2:n) - Vz(1:n-1) <= del_Vz;
    
    %%%%%%%% Constraints definition over!! %%%%%%%%%
    cvx_end
    
    for i = 1:n
        actual_opt_1(i) = ((start_x + sum(Vx(1:i))*del_t - Ox_1(i))^2/mean_sig(1,1) + (start_y + sum(Vy(1:i))*del_t - Oy_1(i))^2/mean_sig(2,2) + (start_z + sum(Vz(1:i))*del_t - Oz_1(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_lin_1(i) = actual_lin_1(i);
        cost_opt_1(i) = ((start_x + sum(vx(1:i))*del_t - Ox_1(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_1(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_1(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_1(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_1(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_1(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        
        actual_opt_2(i) = ((start_x + sum(Vx(1:i))*del_t - Ox_2(i))^2/mean_sig(1,1) + (start_y + sum(Vy(1:i))*del_t - Oy_2(i))^2/mean_sig(2,2) + (start_z + sum(Vz(1:i))*del_t - Oz_2(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_lin_2(i) = actual_lin_2(i);
        cost_opt_2(i) = ((start_x + sum(vx(1:i))*del_t - Ox_2(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_2(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_2(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_2(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_2(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_2(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        
        actual_opt_3(i) = ((start_x + sum(Vx(1:i))*del_t - Ox_3(i))^2/mean_sig(1,1) + (start_y + sum(Vy(1:i))*del_t - Oy_3(i))^2/mean_sig(2,2) + (start_z + sum(Vz(1:i))*del_t - Oz_3(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_lin_3(i) = actual_lin_3(i);
        cost_opt_3(i) = ((start_x + sum(vx(1:i))*del_t - Ox_3(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy_3(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz_3(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox_3(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy_3(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz_3(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
    end
    for i = 1:n
        ratio_1 = (actual_lin_1(i)-actual_opt_1(i))/(cost_lin_1(i)-cost_opt_1(i));
        ratio_2 = (actual_lin_2(i)-actual_opt_2(i))/(cost_lin_2(i)-cost_opt_2(i));
        ratio_3 = (actual_lin_3(i)-actual_opt_3(i))/(cost_lin_3(i)-cost_opt_3(i));
        if(ratio_1<1/4 && ratio_2<1/4 && ratio_3<1/4)
            trust_del_x(i) = (1/4)*trust_del_x(i);
            trust_del_y(i) = (1/4)*trust_del_y(i);
            trust_del_z(i) = (1/4)*trust_del_y(i);
        else
            if(ratio_1>3/4 && ratio_1>3/4 && ratio_1>3/4)
                trust_del_x(i) = min(2*trust_del_x(i), delta_bar);
                trust_del_y(i) = min(2*trust_del_y(i), delta_bar);
                trust_del_z(i) = min(2*trust_del_z(i), delta_bar);
            else
                trust_del_x(i) = trust_del_x(i);
                trust_del_y(i) = trust_del_y(i);
                trust_del_z(i) = trust_del_y(i);
            end
        end
    end 
    guess_vx = Vx;
    guess_vy = Vy;
    guess_vz = Vz;
    if(l~=0)
        if(cost > prev_cost)
            if((cost - prev_cost) < 0.05*ones(1,n))
                l = -1*l;
            else
                clear Vx Vy Vz Px Py Pz
            end
        else
            if((-cost + prev_cost) < 0.05*ones(1,n))
                l = -1*l;
            else
                clear Vx Vy Vz Px Py Pz
            end
        end
    end
    if (l == 0)
        l = l + 1;
        clear Vx Vy Vz Px Py Pz 
    end
    prev_cost = cost;
    clear cost
    iter = iter + 1;
end %% End of while loop
fprintf('Iterations: %d',iter);

%% Printing obstalce avoidance optimization conditions
 for i = 1:n
        %(x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2*(x(i) - Ox(i))*(Px(i) -  x(i)) + 2*(y(i) - Oy(i))*(Py(i) -  y(i)) + 2*(z(i) - Oz(i))*(Pz(i) -  z(i))
        disp(sprintf('Distances: %f   %f   %f ', sqrt((Px(i) - Ox_1(i))^2 + (Py(i) - Oy_1(i))^2 + (Pz(i) - Oz_1(i))^2), sqrt((Px(i) - Ox_2(i))^2 + (Py(i) - Oy_2(i))^2 + (Pz(i) - Oz_2(i))^2), sqrt((Px(i) - Ox_3(i))^2 + (Py(i) - Oy_3(i))^2 + (Pz(i) - Oz_3(i))^2)))
 end

%% Defining stuff for plotting uncertainities

figure('Name', 'Drone_Path') 

%% Starting the scatter plot 
scatter3(0,0,0,5, [1 1 1], 'filled'); 

%% Giving it a title
%% title(strcat('Mahalanobis threshold = ', num2str(R), ',',  ' Blue line: Obstacle,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Confidence interval = ', num2str(confidence*100), '%.'));

h1 = subplot(1,2,1);
%% setting the max/min limits for the axis
set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);

h2 = subplot(1,2,2);
set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
hold on;
pause(2);

for i = 1:n
	tic	
	if i > 1	
		%scatter3(Px(i),Py(i),Pz(i),100, [1 - i/n, 1/2, i/n], 'filled');
		i
		h1 = subplot(1,2,1);
		scatter3(Px(i),Py(i),Pz(i),100, [0, 0, 0], 'filled');
	    	hold on;
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		scatter3(Ox_1(i), Oy_1(i),Oz_1(i), 100,  [1, 0, 1], 'filled');
		hold on;
	    scatter3(Ox_2(i), Oy_2(i),Oz_2(i), 100,  [1, 0, 1], 'filled');
		hold on;
        scatter3(Ox_3(i), Oy_3(i),Oz_3(i), 100,  [1, 0, 1], 'filled');
		hold on;
		%% Drawing confidence ellipses around the obstacle
		[el_x_1, el_y_1, el_z_1] = ellipsoid(Ox_1(i),  Oy_1(i), Oz_1(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_1, el_y_1, el_z_1, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;

        [el_x_2, el_y_2, el_z_2] = ellipsoid(Ox_2(i),  Oy_2(i), Oz_2(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_2, el_y_2, el_z_2, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.9 0.9]);
		hold on;

        [el_x_3, el_y_3, el_z_3] = ellipsoid(Ox_3(i),  Oy_3(i), Oz_3(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_3, el_y_3, el_z_3, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.1 0.1 0.1]);
		hold on;
        
		plot3(Px((i-1):i), Py((i-1):i), Pz((i-1):i), '-', 'Color', 'r', 'LineWidth', 2.0); 
		hold on;
		plot3(Ox_1((i-1):i), Oy_1((i-1):i), Oz_1((i-1):i), '-', 'Color', 'b', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
        plot3(Ox_2((i-1):i), Oy_2((i-1):i), Oz_2((i-1):i), '-', 'Color', 'y', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
        plot3(Ox_3((i-1):i), Oy_3((i-1):i), Oz_3((i-1):i), '-', 'Color', 'g', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
		title(strcat(' Blue line: Obstacle_1,',' Yellow line: Obstacle_2,',' Green line: Obstacle_3,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Bhattacharyya distance = ', num2str(sqrt(R)), '.'));
		
		h2 = subplot(1,2,2);
		scatter3(Px(i),Py(i),Pz(i), 300, [0, 1, 0.5], 'filled');
	    	hold on;
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		scatter3(Ox_1(i), Oy_1(i),Oz_1(i), 100,  [1, 0, 1], 'filled');
		hold on;
	    scatter3(Ox_2(i), Oy_2(i),Oz_2(i), 100,  [1, 0, 1], 'filled');
		hold on;
        scatter3(Ox_3(i), Oy_3(i),Oz_3(i), 100,  [1, 0, 1], 'filled');
		hold on;
	
		%% Drawing confidence ellipses around the obstacle
		[el_x_1, el_y_1, el_z_1] = ellipsoid(Ox_1(i),  Oy_1(i), Oz_1(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_1, el_y_1, el_z_1, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;

        [el_x_2, el_y_2, el_z_2] = ellipsoid(Ox_2(i),  Oy_2(i), Oz_2(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_2, el_y_2, el_z_2, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.9 0.9]);
		hold on;
        [el_x_3, el_y_3, el_z_3] = ellipsoid(Ox_3(i),  Oy_3(i), Oz_3(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x_3, el_y_3, el_z_3, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.1 0.1 0.1]);
		hold on;

		plot3(Px((i-1):i), Py((i-1):i), Pz((i-1):i), '-', 'Color', 'r', 'LineWidth', 2.0); 
		hold on;
		plot3(Ox_1((i-1):i), Oy_1((i-1):i), Oz_1((i-1):i), '-', 'Color', 'b', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
        plot3(Ox_2((i-1):i), Oy_2((i-1):i), Oz_2((i-1):i), '-', 'Color', 'y', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
        plot3(Ox_3((i-1):i), Oy_3((i-1):i), Oz_3((i-1):i), '-', 'Color', 'g', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox_1, Ox_2, Ox_3]))-1, floor(max([Px, Ox_1, Ox_2, Ox_3]))+1], 'YLim', [floor(min([Py, Oy_1, Oy_2, Ox_3]))-1, floor(max([Py, Oy_1, Oy_2, Ox_3]))+1], 'ZLim',[floor(min([Pz, Oz_1, Oz_2, Ox_3])) - 1, floor(max([Pz, Oz_1, Oz_2, Ox_3]))+1]);
		hold on;
		title('Drone movements!!');
		pause(0.2)
		if (i~= n) 
			delete(h2);		
		end
	end


end
