% clear;clc;

function adjoint()

opt = struct('anisotropy', 0.1, 'angle',64, ...
    'nodes', [0 0;1 0;1 1;0 1]', 'minArea', 2e-4);

obj = rte(opt);

f = @(x,y,v) (x);
sigmaS = @(x) (1.5 + 0.6* x(1, :).*x(2, :));
sigmaT = @(x) (1.5 + 0.6*  x(1, :).*x(2,:) + 0.05 );

obj.setBoundaryCondition(f);
obj.setCoefficents(sigmaT, sigmaS);

x1 = obj.ForwardSolve();

% then we compute the adjoint equation.


g = @(x,y ,v) (y);
obj.setBoundaryCondition(g);
obj.setCoefficents(sigmaT, sigmaS);

x2 = obj.ForwardSolve();

% flip directions.
nAngle = obj.nAngle;

x3 = [x2(nAngle/2+1:end, :) ; x2(1:nAngle/2, :)]; % flipped solution (true solution)

% on boundary.
Lbd = size(obj.segms, 2);

s = 0;


for i = 1:Lbd
    cur_edge = obj.segms(:, i);
    
    edge_v = obj.nodes(:, cur_edge(1))   - obj.nodes(:, cur_edge(2));
    
    n = [edge_v(2), -edge_v(1)];
    n = n/norm(n);
    
    l = norm(obj.nodes(:, cur_edge(1))   - obj.nodes(:, cur_edge(2)));
    for j = 1:nAngle
        
        w = [cos((j-1) * obj.dtheta), sin((j-1) * obj.dtheta)];
        
        s = s + 0.5 * x1(j, cur_edge(1)) * x3(j, cur_edge(1)) * (n * w') * obj.dtheta * l;
        s = s + 0.5 * x1(j, cur_edge(2)) * x3(j, cur_edge(2)) * (n * w') * obj.dtheta * l; 
        
    end
    
end

disp(sprintf('error of adjoint %f\n', s));

end


