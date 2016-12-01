global numf numg 



n=1000;

mvec = [3 5 8 12 20]';

for i=1:5
  x = struct('p',ones(n,1)*10);
  lbfgsparams=struct('toler',1.e-4,'maxit',1000,'m',mvec(i));
  [inform,xnew] = LBFGS(@tridia,x,lbfgsparams);
  fprintf('m=%d\t', lbfgsparams.m),
  fprintf('%d\n', inform.iter);
end

  
