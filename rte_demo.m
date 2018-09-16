% rte demo script
%% setup
NA = 128;
nodes = [0 0; 1 0; 1 1; 0 1]';

%% mesh generation
mesh = TriangleMesh();
hull = reshape(nodes, 2 * size(nodes, 2), 1);
idx  = circshift(reshape(repmat(0:size(nodes, 2)-1, 2, 1),...
    2 * size(nodes, 2), 1), [-1,1]);
mesh.set_points_tri(hull);
mesh.set_facets_tri(idx);
mesh = mesh.build_tri(); % not ready to go
mesh = mesh.refine_tri(sprintf('q34.0a%f', 0.001));  

%% mesh properties
[p,s,t,e, n] = mesh.getData_tri();

%% angular space
angle = linspace(0, 2*pi, NA + 1) ;
angle = angle(1:end-1);
dtheta = angle(2) - angle(1);

%% tracing rays, takes time.
tic;rays = tracer(angle, p, t, n);toc;

%% boundary contribution, which should be in solution. 
bcs = BC('dirichlet');
bcs.set_constraint('x - 1.0');
bcs.set_constraint('y - 1.0');
bcs.set_constraint('x');
bcs.set_constraint('y');
[bc1, bc2, bc3, bc4] =  bcs.get_boundary(e, p, 4);
bc = unique([bc1 bc2 bc3 bc4]);

boundary_source = zeros(NA, size(p, 2));

for bid = 1:length(bc)
    for aid = 1:NA
        boundary_source(aid, bc(bid)) = 1;
    end
end

N = 500;
[X,Y] =meshgrid(linspace(0, 1, N));
X = X(:);
Y = Y(:);
ph = phantom('Modified Shepp-Logan', N);
F = scatteredInterpolant(X,Y,ph(:));
sigma_t = 0.8 * F(p(1,:), p(2,:))' + 1;

% sigma_t = ones(size(p, 2), 1);


tic; u = rays.boundary_transport(p, t, sigma_t, boundary_source); toc;
boundaryFluence = (sum(u,1)/NA)';

forwardMap = @(X) ( X - (sum( rays.interior_transport(p, t, sigma_t,  repmat(X', NA, 1)),1)/NA)' );

tic; x = gmres(forwardMap, boundaryFluence, 10, 1e-10, 400); toc;

trisurf(t(1:3,:)', p(1,:), p(2,:), x', 'EdgeColor', 'None');shading interp;view(2);colorbar;colormap jet;

