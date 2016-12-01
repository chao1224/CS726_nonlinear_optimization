global numf numg numH

x = struct('p',[-1.2; 1]);
trparams = struct('maxit',100,'delta',.01,'hatDelta',10,...
                  'eta',.01,'Delta0',1,'toler',1.0e-6);

numf=0; numg=0;
[inform,xnew] = DoglegTR(@obja,x,trparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %10.6g.\n', trparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('  Ending point: '); fprintf('%10.6g ',xnew.p);
fprintf('\n  Ending function value: %10.6g\n', xnew.f);
fprintf('  No. function evaluations: %d, No. gradient evaluations %d\n',...
    numf, numg);
fprintf('  Norm of ending gradient: %10.6g\n\n\n', norm(xnew.g));

%%%%%%%%%%%%%%

x = struct('p',[-1.2; 1]);
trparams = struct('maxit',100,'delta',.01,'hatDelta',10,...
                  'eta',.01,'Delta0',1,'toler',1.0e-6);

numf=0; numg=0;
[inform,xnew] = DoglegTR(@objb,x,trparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %10.6g.\n', trparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('  Ending point: '); fprintf('%10.6g ',xnew.p);
fprintf('\n  Ending function value: %10.6g\n', xnew.f);
fprintf('  No. function evaluations: %d, No. gradient evaluations %d\n',...
    numf, numg);
fprintf('  Norm of ending gradient: %10.6g\n\n\n', norm(xnew.g));

%%%%%%%%%%%%%%%

x = struct('p',[-1.2; 1]);
trparams = struct('maxit',100,'delta',.01,'hatDelta',10,...
                  'eta',.01,'Delta0',1,'toler',1.0e-6);

numf=0; numg=0;
[inform,xnew] = DoglegTR(@objc,x,trparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %10.6g.\n', trparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('  Ending point: '); fprintf('%10.6g ',xnew.p);
fprintf('\n  Ending function value: %10.6g\n', xnew.f);
fprintf('  No. function evaluations: %d, No. gradient evaluations %d\n',...
    numf, numg);
fprintf('  Norm of ending gradient: %10.6g\n\n\n', norm(xnew.g));

