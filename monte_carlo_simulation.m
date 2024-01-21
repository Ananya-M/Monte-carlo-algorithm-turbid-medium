 % BIOE 5648 Biomedical Optics
% Monte Carlo Simulation of photon propagation in turbid media
%
% This code simulates tissue as a 20 x 20 x 20 mm rectangular grid
% with 1mm^3 voxels  that have specific optical properties
%
% The code then simulates the propagation of photons through the tissue
% from a pencil beam incident in the z direction in the middle of the tissue
% and displays the fluence of photons in tissue over 1 ns.

clear; close all; clc;
global epsilon
epsilon = 0.0001;
%% Section 1 Defining the tissue geometry

length_x =20; %??**************************   %tissue geometry in mm
length_y = 20;%??**************************
length_z = 20;%??**************************

dx = 1 ;%??**************************; %voxel size in mm
dy = 1;%??**************************; 
dz = 1;%??**************************;  

dim_x = 20;%??**************************;   %size of tissue matrix in x-direction
dim_y = 20;%??**************************;  
dim_z = 20;%??**************************;  

g = 0.8*ones(dim_x,dim_y,dim_z); %??**************************;   % Anisotropy constant 
n = 1.33*ones(dim_x,dim_y,dim_z); %??**************************;   % Index of refraction 
mu_a = 0.1*ones(dim_x,dim_y,dim_z);%??**************************;   % Array of mu_a values of all 20x20x20 voxels
mu_s = 1.0 *ones(dim_x,dim_y,dim_z);%??**************************;   % Array of all mu_s' values

 
mu_a(:,:,3:5) = 1;%??**************************; Problem 10. changing some of the optical properties
mu_s(:,:,3:5) = 2.5; %??**************************;
g(:,:,3:5) =0; %??**************************;


%% Section 2 Defining initial light simulation conditions

photon_total = 1e5;  % total number of photons for simulation

c = 3e11;%??**************************; % Speed of light in mm/sec 
final_T = 40; %??**************************;  % total simulation time in ps
dt = 1 ;%??**************************;  % 1 ps intervals
dim_t = 40;%??**************************;
final_D = 12; %??**************************; % Threshold for total distance traveled by a photon  (speed of light * T_final in ps)


photon_launch.p = [10.5 10.5 eps(single(1.0))]; %initial location of photon at tissue surface p = [p.x, p.y, p.z]
photon_launch.v = [0 0 1]; %??**************************;  % initial direction of photon upon launch (think of it as a unit velocity vector with no length dimension)
photon_launch.w = 1;%??**************************;   % initial weight of photon
photon_launch.t = 0;%??**************************;   % initial photon time

m_roulette = 50 ;%??**************************;     % threshold for photon termination ("Russian Roulette") weight




%% Section 3 Simulating photon propagation
% refer to the flow chart in the assignment directions for

Fluence = zeros (dim_x, dim_y, dim_z, dim_t);
tic
for i_photon = 1:photon_total
    
    photon = photon_launch;  %set this new photon to have the initial position, direction, weight, time
    
    s = -log(rand());     % calculate dimensionless scattering length
    
%     roulette_rand= -log(rand());
    
    while (s> 0)
        
        
        
        %first lets assess where the photon is (and "when.") Get the current
        % voxel location and the corresponding optical properties
        current_voxel = ceil(photon.p);
        mu_s_step = mu_s(current_voxel(1), current_voxel(2), current_voxel(3));
        mu_a_step = mu_a(current_voxel(1), current_voxel(2), current_voxel(3));
        n_step = n(current_voxel(1), current_voxel(2), current_voxel(3));
        current_t = 1+ floor(photon.t);
        
        
        % call the "one_move_in_cube" function to figure out
        % the distance to advance a photon, and into which neighboring voxel
        
        [d_step, p_next] = one_move_in_cube(photon.p,photon.v);  %calculate minimum distance required to enter a neighboring voxel and new photon position
        
        
        % "Advance" the photon.
        % 1. decrement the remaining scattering length, s
        % 2. calculate "weight loss" from absorption
        % 3. Add the "weight loss" to the Fluence matrix
        %
        
        %??**************************;
        d_incr = min(d_step,s/mu_s_step);%??************************** assign the distance to increment 
        %??**************************;
        
        
        %??**************************;
        if d_incr == d_step %??**************************;
            s = s-d_incr* mu_s_step;%??**************************;    %decrement the length of this scattering event
        else
            s = 0; %??**************************;
        end
        
        delta_w = photon.w *(1-exp(-mu_a_step * d_incr));  %calculate "weight loss" from absorption
        photon.w = photon.w - delta_w;   % decrement the photon weight

    
        
        
        %??**************************;
        % increment the Fluence matrix at the current voxel and time
        %??**************************;
        Fluence(current_voxel(1),current_voxel(2),current_voxel(3),current_t)= Fluence(current_voxel(1),current_voxel(2),current_voxel(3),current_t)+delta_w;
        %??**************************;
        % update the photon position 
        %??**************************;
        photon.p= p_next;
        
        t_incr = d_incr* (n_step/c)*1e12;  % calculate and increment time
        photon.t = photon.t + t_incr;
        
        
        
        %check if we've reached the time limit
        
        if photon.t>= final_T
            break;
        end
        
        %check if the photon has dropped below the threshold weight
        % if so, the photon engages in "Russian Roulette"
        %??**************************;
        % Roulette condition should be entered here
        %??**************************;
        if photon.w<= 1/m_roulette
            if rand < 1/m_roulette
                photon.w = photon.w*m_roulette
            end
        else
             Fluence(current_voxel(1),current_voxel(2),current_voxel(3),current_t)=Fluence(current_voxel(1),current_voxel(2),current_voxel(3),current_t)+photon.w;
        
            
       end   

        

        %check if the photon has escaped the tissue
        
        if (p_next(1) < 0) ||  (p_next(1) > size(Fluence,1))
            break;
        end
        
        if (p_next(2) < 0) ||  (p_next(2) > size(Fluence,2))
            break;
        end
        
        if (p_next(3) < 0) ||  (p_next(3) > size(Fluence,3))
            break;
        end
        
        
        %if the end of the current scattering event has been reached,
        %calculate a new scattering direction and length
        if s<= 0  
            v_new = update_scattering_direction(photon.v, g(current_voxel(1), current_voxel(2), current_voxel(3))) ;
            s = -log(rand());     % calculate new dimensionless scattering length
            photon.v = v_new;
        end
        
        
    end
end
toc

%% Section 4 displaying the data
f = squeeze(Fluence(:,11,:,:));

f_int = zeros (size(f,1,2,3));

h = figure(); set (gca, 'FontWeight', 'Bold');

for index_t = 1:final_T
    
        f_int (:,:,index_t) = rot90(sum(f(:,:,1:index_t),3),3);
    
%     f_int (:,:,index_t) = rot90(f(:,:,index_t),3);
    h;
    imagesc(log10(squeeze(f_int(:,:,index_t)))); axis image;
    xlabel('X (mm)', 'FontWeight', 'Bold'); ylabel ('Depth (mm)', 'FontWeight', 'Bold'); title(['log(Fluence), y = 11; t = ' num2str(index_t) ' ps']);
    colorbar
    drawnow;
    pause (0.25);
    frame = getframe(h);
    im{index_t} = frame2im(frame);
end


filename = 'sample.gif'; % Specify the output file name
for idx = 1:final_T
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.25);
    end
end




%% Section 5

function [dist, htime] = one_move_in_cube(p,v)
% takes the position and movement direction of the photon and calculates
% which of the neighboring voxels to step into
% the shortest (unitless) distances to the voxel walls (in x y or z)  is computed
%

% calculate "time of fly" to hit the wall in each direction


% Move photon to the first intersecting wall inside a 1x1x1mm cube
% Input:
%     p: {x,y,z}current x/y/z position of the photon
%     v: {vx,vy,vz}, direction unitary vector
% Output:
%     htime: store the intersection position (x/y/z) of the first wall
%     id:  id=0: intersect with x-plane; id=1: y-plane; id=3: z-plane


global epsilon;


htime(1)=abs((floor(p(1))+(v(1)>0)-p(1))/v(1)); % Calculate the "time of flight" to the next voxel wall in x, y, z direction
htime(2)=abs((floor(p(2))+(v(2)>0)-p(2))/v(2));
htime(3)=abs((floor(p(3))+(v(3)>0)-p(3))/v(3));



% get the direction with smallest time of fly
[dist, dir]= min(htime);% gives out minimum distance a photon will have to travel to get to voxel wall and the direction plane (dist = htime*norm(v) = htime*1)
% xdirection: dir = 1;  ydirection: dir = 2;
% zdirection dir = 3;

% Once its determine which direction is easiest to move with least
% distance, then it will move the photon to that direction and update its
% x,y,z positions
htime(1) = p(1)+dist*v(1); % update the photon position in x direction
htime(2) = p(2)+dist*v(2); % upate the photon position in y direction
htime(3) = p(3)+dist*v(3); % update the photon position in z direction

xi(1)=round(htime(1))+epsilon*((v(1)>0)-(v(1)<0));
xi(2)=round(htime(2))+epsilon*((v(2)>0)-(v(2)<0));
xi(3)=round(htime(3))+epsilon*((v(3)>0)-(v(3)<0));


htime(dir) = xi(dir);
end



% Function to calculate new scattering direction vector with random # and
% Henyey Greenstein equation

function [v_new] = update_scattering_direction(v, g)

if g ==0;
    cos_theta = 2*rand - 1;
else
    
    dummy_var = (1-g^2) / (1-g+2*g*rand);
    cos_theta = (1+g^2 - dummy_var^2)/ (2*g);  % Henyey Greenstein CDF
end

theta = acos(cos_theta);   % Calculate new Zenith angle

phi = 2 * pi * rand;    %Calculate new Azimuthal Angle


%??**************************;
% compute A, B, and the new values for the 3 v_new vector
B = 1 - v(3)^2  ;%??**************************;
A = sin(theta)/sqrt(B) ;%??**************************;

if v(1)~=0 || v(2) ~= 0
    v(1)= A *(v(1)*v(3)*cos(phi)-v(2)*sin(phi))+v(1)*cos(theta);
    v(2)= A *(v(2)*v(3)*cos(phi)+v(1)*sin(phi))+v(2)*cos(theta);
    v(3)= -A*B*cos(phi)+v(3)*cos(theta);
else
    v(1)=sin(theta)*cos(phi);
    v(2)=sin(theta)*sin(phi);
    v(3)=sign(v(3))*cos(theta);
end
v_new=[v(1) v(2) v(3)];
v_new = v_new/norm(v_new);  %ensure that new direction is a unit vector (length = 1)
end


