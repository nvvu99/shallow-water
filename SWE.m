%% Model
% dh/dt    + d(hu)/dx          + d(hv)             = 0
% d(hu)/dt + d(hu^2+gh^2/2)/dx + d(huv)/dy         = 0
% d(hv)/dt + d(huv)/dx         + d(hv^2+gh^2/2)/dy = 0

clear all
clc

%Constants
g = 9.81; % units of gravity are m/s^2
EXPORT_VIDEO = false;

%Grid length, number and spacing
Lx = 10;
Ly = 10;
nx = 52;
ny = 52;
dx = Lx / (nx - 2);
dy = Ly / (ny - 2);

% set up finite-difference mesh or grid:
[x y] = meshgrid(linspace(0,Lx,nx), linspace(0,Ly,ny));


%% Initial condition. height at t=0.
u = zeros(nx,ny); 
v = zeros(nx,ny);
H = 5;
height = zeros(nx,ny);

%Move the initial column of water around by changing io and jo. 
%w will change the width of the column
io = 38;
jo = 19;
w = 5;

for i = 1:nx 
    for j = 1:ny
        height(i,j) = H + 5 * exp((-((i-io)^2 + (j-jo)^2))/(w^2));
    end
end

% U = (h, hu, hv)
U_current = zeros(nx, ny, 3);
U_current(:,:,1) = height;
U_previous = U_current;
U_next = U_current;


%%  Move through time.  Make plots.
Nsteps = 500;
t = 0; dt = 0.01; 
i = 2:nx-1; left = 1:nx-2; right = 3:nx;
j = 2:ny-1; up   = 1:ny-2; down  = 3:ny;

fig = figure(1);
fig.Position = [0 0 800 600];

surf(x(i,j), y(i,j), U_current(i,j,1))
% zlim([0 15])
% shading interp
% pause;

if EXPORT_VIDEO
    video_writer = VideoWriter('SWE_conservative_form.avi');
    open(video_writer);
end

for n=1:Nsteps
    huv = U_current(:,:,2) .* U_current(:,:,3) ./ U_current(:,:,1);
    huu = U_current(:,:,2).^2 ./ U_current(:,:,1);
    hvv = U_current(:,:,3).^2 ./ U_current(:,:,1);
    ghh = g * U_current(:,:,1).^2;
    % calculate F = (hu, hu^2 + gh^2/2, huv)
    F = cat(3, U_current(:,:,2), huu + ghh/2, huv);
    % calcualte G = (hv, huv, hv^2 + gh^2/2)
    G = cat(3, U_current(:,:,3), huv, hvv + ghh/2);

    U_next(i,j,:) = U_previous(i,j,:) ...
        - dt/dx * (F(i,right,:) - F(i,left,:))...
        - dt/dy * (G(up,j,:) - G(down,j,:));

    % impose boundary conditions on h
    % reflection
    U_next(i,end-1,1) =  U_next(i,end-2,1); U_next(i,2,1) =  U_next(i,3,1);
    U_next(end-1,j,1) =  U_next(end-2,j,1); U_next(2,j,1) =  U_next(3,j,1);
    % on hu
    U_next(i,end-1,2) = -U_next(i,end-2,2); U_next(i,2,2) = -U_next(i,3,2);
    U_next(end-1,j,2) =  U_next(end-2,j,2); U_next(2,j,2) =  U_next(3,j,2);
    % on hv
    U_next(i,end-1,3) =  U_next(i,end-2,3); U_next(i,2,3) =  U_next(i,3,3);
    U_next(end-1,j,3) = -U_next(end-2,j,3); U_next(2,j,3) = -U_next(3,j,3);

    % u_next = U_next(:,:,2)./U_next(:,:,1);   v_next = U_next(:,:,3)./U_next(:,:,1);

    figure(1)

    surf(x(i,j),y(i,j),U_current(i,j,1))

    zlim([0 15])
    shading interp

    if EXPORT_VIDEO
        frame = getframe(gcf);
        writeVideo(video_writer, frame);
    end

    U_previous = U_current;
    U_current = U_next;
end

if EXPORT_VIDEO
    close(video_writer);
end
