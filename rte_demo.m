% rte demo script
%% setup
N = 128; NA = 128;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));
nodes = [cos(theta); sin(theta)];

%% mesh generation
mesh = TriangleMesh();
hull = reshape(nodes, 2 * size(nodes, 2), 1);
idx  = circshift(reshape(repmat(0:size(nodes, 2)-1, 2, 1),...
    2 * size(nodes, 2), 1), [-1,1]);
mesh.set_points_tri(hull);
mesh.set_facets_tri(idx);
mesh = mesh.build_tri(); % not ready to go
mesh = mesh.refine_tri(sprintf('q34.0a%f', 0.00125));  

%% mesh properties
[p,s,t,e, n] = mesh.getData_tri();

%% angular space
angle = linspace(0, 2*pi, NA + 1);
angle = angle(1:end-1);

%% tracing rays, takes time.
tic;rays = tracer(angle, p, t, n);toc;

%% boundary contribution, which should be in solution. 
bc = zeros(NA, size(p,2));
bc(:, 1) = 1.0; % isotropic point source.

tic;u = rays.boundary_transport(p, t, 0.002*ones(size(p,2),1), bc);toc;
fl = sum(u)/NA;
trisurf(t(1:3,:)', p(1,:), p(2,:), fl', 'EdgeColor', 'none');shading interp;view(2);colorbar;
