
clc;
clear all;
close all;


%% Drone is supposed to reach (dest_x, dest_y, dest_z) from (0,0,0) in
%% 'n' timesteps with each timestep of duration del_t
%% Version 2 (for higher values of n)
n = 40;      % Number of timesteps
del_t = 1.0; % Time duration for which some command is executed.
tot_time = n*del_t; %% Total time of flight

%% Defining obstacle trajectory
%% For now we will assume that obstacle is moving along a straing line. it's starting from (0,0,0) and It's trajectory will be determined based on the destination of the drone

%% Drone Destination variables
dest_x = 12;
dest_y = 12;
dest_z = 12;

%%%%%%%%%%%%%% Obstacle stuff starts here!!!! %%%%%%%%%%%%%%%%%
%% Here, the drone is moving from (0,0,0) to (dest_x, dest_y, dest_z) in n timesteps
%% The obstacle will reach from start_obs to end_obs in n timesteps 
%%Ob_start=[0, dest_y/2, dest_z];
%%Ob_end=[dest_x, dest_y/2, 0];
Ob_start=[0,5,0];
Ob_end=[dest_x, dest_y+5, dest_z];

%% Calculating velocity of the obstacle in all the directions
vel_obs = (Ob_end - Ob_start)/tot_time;

%% Defining Obstacle path location at each timestep
for i = 1:n
  Ox(i) = Ob_start(1) + del_t*i*vel_obs(1);
  Oy(i) = Ob_start(2) + del_t*i*vel_obs(2);
  Oz(i) = Ob_start(3) + del_t*i*vel_obs(3);
end

%% The confidence has to be between [0, 1]
confidence = 0.68;

%% defining drone uncertainty
drone_sig = diag([1.5^2 1.5^2 1.5^2]);

%% Defining obstalce uncertainty
obs_sig1 = diag([1.5^2 1.5^2 1.5^2]);

%% Getting mean sigma of the obstacle and the robot
mean_sig = (obs_sig1 + drone_sig)/2;

%% getting Mahalanobis threshold based on confidence interval
mahalanobisThresh = chi2inv(confidence, 2);

%% Radius of avoidance
R = 1.0; %% This is a Bhattacharyya radius

%% Calculating ev_obs1 of the covariance matrix.
%% The ev_obs1 will give us lengths of semi-axis lengths
sig_obs1_x = obs_sig1(1,1);
sig_obs1_y = obs_sig1(2,2);
sig_obs1_z = obs_sig1(3,3);
ev_obs1 = eig(obs_sig1);

sig_drone_x = drone_sig(1,1);
sig_drone_y = drone_sig(2,2);
sig_drone_z = drone_sig(3,3);
ev_drone = eig(drone_sig);
%% Getting ev_obs1 for X,Y,Z direction
%% The eigen values returned from eig() are sorted in decreasing order
if(sig_obs1_x <= sig_obs1_y && sig_obs1_x <= sig_obs1_z)
	eig_obs1_x = ev_obs1(1,1);

	if(sig_obs1_y < sig_obs1_z)
		eig_obs1_y = ev_obs1(2,1);
		eig_obs1_z = ev_obs1(3,1);
	else
		eig_obs1_y = ev_obs1(3,1);
		eig_obs1_z = ev_obs1(2,1);
	end
end

if(sig_obs1_y <= sig_obs1_z && sig_obs1_y <= sig_obs1_x)
	eig_obs1_y = ev_obs1(1,1);

	if(sig_obs1_x < sig_obs1_z)
		eig_obs1_x = ev_obs1(2,1);
		eig_obs1_z = ev_obs1(3,1);
	else
		eig_obs1_x = ev_obs1(3,1);
		eig_obs1_z = ev_obs1(2,1);
	end
end

if(sig_obs1_z <= sig_obs1_x && sig_obs1_z <= sig_obs1_y)
	eig_obs1_z = ev_obs1(1,1);

	if(sig_obs1_x < sig_obs1_y)
		eig_obs1_x = ev_obs1(2,1);
		eig_obs1_y = ev_obs1(3,1);
	else
		eig_obs1_x = ev_obs1(2,1);
		eig_obs1_y = ev_obs1(3,1);
	end
end

%% Drone Eigen Value ordering
if(sig_drone_x <= sig_drone_y && sig_drone_x <= sig_drone_z)
	eig_drone_x = ev_drone(1,1);

	if(sig_drone_y < sig_drone_z)
		eig_drone_y = ev_drone(2,1);
		eig_drone_z = ev_drone(3,1);
	else
		eig_drone_y = ev_drone(3,1);
		eig_drone_z = ev_drone(2,1);
	end
end

if(sig_drone_y <= sig_drone_z && sig_drone_y <= sig_drone_x)
	eig_drone_y = ev_drone(1,1);

	if(sig_drone_x < sig_drone_z)
		eig_drone_x = ev_drone(2,1);
		eig_drone_z = ev_drone(3,1);
	else
		eig_drone_x = ev_drone(3,1);
		eig_drone_z = ev_drone(2,1);
	end
end

if(sig_drone_z <= sig_drone_x && sig_drone_z <= sig_drone_y)
	eig_drone_z = ev_drone(1,1);

	if(sig_drone_x < sig_drone_y)
		eig_drone_x = ev_drone(2,1);
		eig_drone_y = ev_drone(3,1);
	else
		eig_drone_x = ev_drone(2,1);
		eig_drone_y = ev_drone(3,1);
	end
end


%% Getting lengths of semi major axis based on confidence values

%% Radius of ellipse along X-direction
rad1_X = sqrt(mahalanobisThresh*eig_obs1_x); 

%% Radius of ellipse along Y-direction
rad1_Y = sqrt(mahalanobisThresh*eig_obs1_y); 

%% Radius of ellipse along Z-direction
rad1_Z = sqrt(mahalanobisThresh*eig_obs1_z); 

%% Getting lengths of semi major axis for Drone based on confidence values
radD_X = sqrt(mahalanobisThresh*eig_drone_x);
radD_Y = sqrt(mahalanobisThresh*eig_drone_y);
radD_Z = sqrt(mahalanobisThresh*eig_drone_z);
%%%%%%%%%Obstacle stuff ends here %%%%%%%%%%%%%%%

%% Difference between subsequent velocities
del_Vx = 10;
del_Vy = 10;
del_Vz = 10;

l = 0;
no_of_iter = 2;

while l<no_of_iter
    
    cvx_begin
    
    variables Vx(n) Vy(n) Vz(n) %% Velocities
    %variables Px(n) Py(n) Pz(n) %% Drone locations at each timestep
    
    for i = 1:n
        Px(i)  =  sum(Vx(1:i))*del_t;
        Py(i)  =  sum(Vy(1:i))*del_t;
        Pz(i)  =  sum(Vz(1:i))*del_t;
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
    if(l~=0)
    	for i = 1:n
    	   %% Linearized Bhattacharyya constraint	2.9957-95%  4.6052-99%  1.1394-68%
    	  2.9957 >= ((x(i) - Ox(i))^2/mean_sig(1,1) + (y(i) - Oy(i))^2/mean_sig(2,2) + (z(i) - Oz(i))^2/mean_sig(3,3) + 2*(x(i) - Ox(i))*(Px(i) -  x(i))/mean_sig(1,1) + 2*(y(i) - Oy(i))*(Py(i) -  y(i))/mean_sig(2,2) + 2*(z(i) - Oz(i))*(Pz(i) -  z(i))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig1)))/2.0 >= 1.1394
    	   
    	end
    end
    %%%%%%%% Constraints %%%%%%%%
    Vx >= 0;
    Vy >= 0;
    Vz >= 0;
    
    %% The velocity during initialization shouldn't be much higher
    0 <= Vx(1) <= 2.5;
    0 <= Vy(1) <= 2.5;
    0 <= Vz(1) <= 2.5;
    
    %% Ensuring that subsequent X vel commands do not have difference more than del_Vx
    0 <= Vx(2:n) - Vx(1:n-1) <= del_Vx;
    
    %% Ensuring that subsequent Y vel commands do not have difference more than del_Vy
    0 <= Vy(2:n) - Vy(1:n-1) <= del_Vy;
    
    %% Ensuring that subsequent Z vel commands do not have difference more than del_Vz
    0 <= Vz(2:n) - Vz(1:n-1) <= del_Vz;
    
    %%%%%%%% Constraints definition over!! %%%%%%%%%
    cvx_end
 
    %% Updating the linearization points!
    x = Px;
    y = Py;
    z = Pz;
    %cvx_optval
    if(l ~= no_of_iter-1)
        clear Vx Vy Vz Px Py Pz
    end	

    l = l + 1;
    R = R + 1
end %% End of while loop

%% Printing obstalce avoidance optimization conditions
 for i = 1:n
        %(x(i) - Ox(i))^2 + (y(i) - Oy(i))^2 + (z(i) - Oz(i))^2 + 2*(x(i) - Ox(i))*(Px(i) -  x(i)) + 2*(y(i) - Oy(i))*(Py(i) -  y(i)) + 2*(z(i) - Oz(i))*(Pz(i) -  z(i))
        %disp(sqrt((Px(i) - Ox(i))^2 + (Py(i) - Oy(i))^2 + (Pz(i) - Oz(i))^2))
        ((x(i) - Ox(i))^2/mean_sig(1,1) + (y(i) - Oy(i))^2/mean_sig(2,2) + (z(i) - Oz(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig1)))/2.0 
 end

%% Defining stuff for plotting uncertainities

figure('Name', 'Drone_Path') 

%% Starting the scatter plot 
scatter3(0,0,0,5, [1 1 1], 'filled'); 

%% Giving it a title
%% title(strcat('Mahalanobis threshold = ', num2str(R), ',',  ' Blue line: Obstacle,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Confidence interval = ', num2str(confidence*100), '%.'));

h1 = subplot(1,2,1);
%% setting the max/min limits for the axis
set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]);

h2 = subplot(1,2,2);
set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]);
hold on;
pause(1);

for i = 1:n
	tic	
	if i > 1	
		%scatter3(Px(i),Py(i),Pz(i),100, [1 - i/n, 1/2, i/n], 'filled');
		i
		h1 = subplot(1,2,1);
		scatter3(Px(i),Py(i),Pz(i), 100, [0, 1, 0]);
	    	hold on;
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		scatter3(Ox(i), Oy(i),Oz(i), 100,  [1, 0, 1], 'filled');
		hold on;
	
		%% Drawing confidence ellipses around the obstacle1
		[el1_x, el1_y, el1_z] = ellipsoid(Ox(i),  Oy(i), Oz(i), rad1_X, rad1_Y, rad1_Z, 50);
		surf(el1_x, el1_y, el1_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;
        
        %% Drawing confidence ellipses around the Drone
		[elD_x, elD_y, elD_z] = ellipsoid(Px(i),  Py(i), Pz(i), radD_X, radD_Y, radD_Z, 50);
		surf(elD_x, elD_y, elD_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;

		plot3(Px((i-1):i), Py((i-1):i), Pz((i-1):i), '-', 'Color', 'r', 'LineWidth', 2.0); 
		hold on;
		plot3(Ox((i-1):i), Oy((i-1):i), Oz((i-1):i), '-', 'Color', 'b', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]); 
		hold on;
		title(strcat(' Blue line: Obstacle,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Bhattacharyya distance = ', num2str(sqrt(R)), '.'));
		
		h2 = subplot(1,2,2);
		scatter3(Px(i),Py(i),Pz(i), 300, [0, 1, 0.5], 'filled');
	    	hold on;
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		scatter3(Ox(i), Oy(i),Oz(i), 100,  [1, 0, 1], 'filled');
		hold on;
	
		%% Drawing confidence ellipses around the obstacle
		[el1_x, el1_y, el1_z] = ellipsoid(Ox(i),  Oy(i), Oz(i), rad1_X, rad1_Y, rad1_Z, 50);
		surf(el1_x, el1_y, el1_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;
        
        %% Drawing confidence ellipses around the Drone
		[elD_x, elD_y, elD_z] = ellipsoid(Px(i),  Py(i), Pz(i), radD_X, radD_Y, radD_Z, 50);
		surf(elD_x, elD_y, elD_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
		hold on;

		plot3(Px((i-1):i), Py((i-1):i), Pz((i-1):i), '-', 'Color', 'r', 'LineWidth', 2.0); 
		hold on;
		plot3(Ox((i-1):i), Oy((i-1):i), Oz((i-1):i), '-', 'Color', 'b', 'LineWidth', 2.0);
		set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]); 
		hold on;
		title('Drone movements!!');
		pause(0.2)
		if (i~= n) 
			delete(h2);		
		end
    end
end
figure('Name','Start')
set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]);
scatter3(Px(1),Py(1),Pz(1), 100, [0, 1, 0]);
hold on;

%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
scatter3(Ox(1), Oy(1),Oz(1), 100,  [1, 0, 1], 'filled');
hold on;

%% Drawing confidence ellipses around the obstacle1
[el1_x, el1_y, el1_z] = ellipsoid(Ox(1),  Oy(1), Oz(1), rad1_X, rad1_Y, rad1_Z, 50);
surf(el1_x, el1_y, el1_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
hold on;

%% Drawing confidence ellipses around the Drone
[elD_x, elD_y, elD_z] = ellipsoid(Px(1),  Py(1), Pz(1), radD_X, radD_Y, radD_Z, 50);
surf(elD_x, elD_y, elD_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
hold on;
title(strcat(' Blue line: Obstacle,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Bhattacharyya distance = ', num2str(sqrt(R)), '.'));

%%End
figure('Name','End')
set(gca,'XLim', [floor(min([Px, Ox]))-1, floor(max([Px, Ox]))+1], 'YLim', [floor(min([Py, Oy]))-1, floor(max([Py, Oy]))+1], 'ZLim',[floor(min([Pz, Oz])) - 1, floor(max([Pz, Oz]))+1]);
scatter3(Px(n),Py(n),Pz(n), 100, [0, 1, 0]);
hold on;

%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
scatter3(Ox(n), Oy(n),Oz(n), 100,  [1, 0, 1], 'filled');
hold on;

%% Drawing confidence ellipses around the obstacle1
[el1_x, el1_y, el1_z] = ellipsoid(Ox(n),  Oy(n), Oz(n), rad1_X, rad1_Y, rad1_Z, 50);
surf(el1_x, el1_y, el1_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
hold on;

%% Drawing confidence ellipses around the Drone
[elD_x, elD_y, elD_z] = ellipsoid(Px(n),  Py(n), Pz(n), radD_X, radD_Y, radD_Z, 50);
surf(elD_x, elD_y, elD_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
hold on;
title(strcat(' Blue line: Obstacle,', ' Red line: Drone path,', ' Total time = ', num2str(tot_time), ' Bhattacharyya distance = ', num2str(sqrt(R)), '.'));