opt = struct('anisotropy', 0.8);

obj = rte(opt);

f = @(x,y,v) 1;
sigmaS = @(x) (5 + 0.2 * rand(size( x(1, :))));
sigmaT = @(x) (0.2 * rand(size( x(1, :))) + 5 + x(2,:) * 0.2);

obj.setBoundaryCondition(f);
obj.setCoefficents(sigmaT, sigmaS);

x1 = obj.ForwardSolve();

obj.plot(x1);