%%
%% Optimization based code for flying Pelican
%% The optimization satisfied the velocity constraints
%% It tries to minimize the error of cost function
%%

%%
%% First go to terminal and run following
%% cd rotors_ws
%% source devel/setup.bash
%% roslaunch rotors_gazebo mav_hovering_example.launch mav_name:=pelican
%% Only after running that, you should run this code.
%% Don't forget to bring a sphere into the gazebo simulation.
%% Use gazebo top bar to get a sphere, on the left column, 
%% go to models, unit_sphere_0 and turn off the gravity factor!
%%

%%
%% The code integrates position controller and uses a velocity controller.
%% The code assumes that the starting position of drone is (0,0,0)
%% The updated version will be more generic in nature!
%% This code uses Bhattacharyya distance as a metric to avoid Obstalcle
%% Here, both, the obstacle as well as the robot(drone) have gaussian uncertainties
%%

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
Ob_start=[0, dest_y/2, dest_z];
Ob_end=[dest_x, dest_y/2, 0];
%Ob_start=[-3,-3,0];
%Ob_end=[12,12,6];

%% Calculating velocity of the obstacle in all the directions
vel_obs = (Ob_end - Ob_start)/tot_time;

%% Defining Obstacle path location at each timestep
for i = 1:n
  Ox(i) = Ob_start(1) + del_t*i*vel_obs(1);
  Oy(i) = Ob_start(2) + del_t*i*vel_obs(2);
  Oz(i) = Ob_start(3) + del_t*i*vel_obs(3);
  %Ox(i) = dest_x/2;
  %Oy(i) = dest_y/2;
  %Oz(i) = 1;
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

guess_vx = 0.07*ones(n,1);
guess_vy = 0.04*ones(n,1);
guess_vz = 0.06*ones(n,1);
trust_del_x = 0.25*ones(n,1);
trust_del_y = 0.25*ones(n,1);
trust_del_z = 0.25*ones(n,1);

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
        actual_lin(i) = ((start_x + sum(vx(1:i))*del_t - Ox(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost(i) = ((start_x + sum(vx(1:i))*del_t - Ox(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        save_cost = cost;
        cost(i) >= R
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
    %%	Defining path of the obstacle based on first optimization
    %if l ==1
    %	    Ob_middle=[Px(n/2), Py(n/2), Pz(n/2)];
    %	    Ob_end = 2*Ob_middle - Ob_start; %% This will be determined based on first velocity solution of the drone
    %end
    %% Updating the linearization points!
% %     x = Px;
% %     y = Py;
% %     z = Pz;
    %cvx_optval
% %     if(l ~= no_of_iter-1)
% %         clear Vx Vy Vz Px Py Pz
% %     end
    %%R = R + 1
    for i = 1:n
        actual_opt(i) = ((start_x + sum(Vx(1:i))*del_t - Ox(i))^2/mean_sig(1,1) + (start_y + sum(Vy(1:i))*del_t - Oy(i))^2/mean_sig(2,2) + (start_z + sum(Vz(1:i))*del_t - Oz(i))^2/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
        cost_lin(i) = actual_lin(i);
        cost_opt(i) = ((start_x + sum(vx(1:i))*del_t - Ox(i))^2/mean_sig(1,1) + (start_y + sum(vy(1:i))*del_t - Oy(i))^2/mean_sig(2,2) + (start_z + sum(vz(1:i))*del_t - Oz(i))^2/mean_sig(3,3) + 2*(start_x + sum(vx(1:i))*del_t - Ox(i))*(sum(Vx(1:i)) - sum(vx(1:i)))/mean_sig(1,1) + 2*(start_y + sum(vy(1:i))*del_t - Oy(i))*(sum(Vy(1:i)) - sum(vy(1:i)))/mean_sig(2,2) + 2*(start_z + sum(vz(1:i))*del_t - Oz(i))*(sum(Vz(1:i)) - sum(vz(1:i)))/mean_sig(3,3))/8.0 + log(det(mean_sig)/sqrt(det(drone_sig)*det(obs_sig)))/2.0;
    end
    for i = 1:n
        ratio = (cost_opt(i)-cost_lin(i))/(actual_opt(i)-actual_lin(i));
        if(ratio>3/4)
            trust_del_x(i) = trust_del_x(i) + 0.01;
            trust_del_y(i) = trust_del_y(i) + 0.01;
            trust_del_z(i) = trust_del_z(i) + 0.01;
        elseif(ratio<1/4)
            trust_del_x(i) = trust_del_x(i) - 0.01;
            trust_del_y(i) = trust_del_y(i) - 0.01;
            trust_del_z(i) = trust_del_z(i) - 0.01;
        end
    end 
    guess_vx = Vx;
    guess_vy = Vy;
    guess_vz = Vz;
    if(l~=0)
        if(cost > prev_cost)
            if((cost - prev_cost) < 0.01*ones(1,n))
                l = -1*l;
            else
                clear Vx Vy Vz Px Py Pz
            end
        else
            if((-cost + prev_cost) < 0.01*ones(1,n))
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
        disp(sqrt((Px(i) - Ox(i))^2 + (Py(i) - Oy(i))^2 + (Pz(i) - Oz(i))^2))
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
pause(10);

for i = 1:n
	tic	
	if i > 1	
		%scatter3(Px(i),Py(i),Pz(i),100, [1 - i/n, 1/2, i/n], 'filled');
		i
		h1 = subplot(1,2,1);
		scatter3(Px(i),Py(i),Pz(i),100, [0, 0, 0], 'filled');
	    	hold on;
		%scatter3(Ox(i), Oy(i),Oz(i), 100,  [1 - i/n, 1/2, i/n], 'filled');
		scatter3(Ox(i), Oy(i),Oz(i), 100,  [1, 0, 1], 'filled');
		hold on;
	
		%% Drawing confidence ellipses around the obstacle
		[el_x, el_y, el_z] = ellipsoid(Ox(i),  Oy(i), Oz(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x, el_y, el_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
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
		[el_x, el_y, el_z] = ellipsoid(Ox(i),  Oy(i), Oz(i), rad_X, rad_Y, rad_Z, 50);
		surf(el_x, el_y, el_z, 'FaceAlpha', 0.4, 'EdgeColor','none', 'FaceColor', [0.9 0.5 0.1]);
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



