function spdprotein(node, node_path, outputs_path, solver, max_iter)
% ''' Solves the Relaxad SPD programm for a given distance list'''
fprintf('# Protein: %s || max iterations: %s \n', node, string(max_iter))

% Enale Matlab Log option
log_path = sprintf('%s\\SPD_log.txt', outputs_path);

% Get distance data file from the respective node;
dist_data = readtable(sprintf('%s\\dist.txt', node_path));
disp('# Distance file read: Complete')

% Using the data, creating an acess matrix as [i,j, l_ij, u_ij]
Data = [dist_data(:,1), dist_data(:,2), dist_data(:,9), dist_data(:,10)];

clearvars dist_data

% table ~~> array
Data = table2array(Data);

% Initiate basic variables:
num_dist = length(Data(:, 1));
num_atom = Data(num_dist, 2);
u = Data(:, 1);
v = Data(:, 2);

% Inicialization for guide matrix H (contraints matrix)
H = zeros(num_dist, num_atom^2);
z =  [1, -1, -1, 1];

for i = 1 : num_dist
    key = [u(i), v(i), u(i), v(i)];
    target = [u(i), u(i), v(i), v(i)];
    E = sparse(key, target, z, num_atom, num_atom);
    transposed = E';
    H(i,:) =  transposed(:)';
end

clearvars Z key target E

% Restriction inicialization
lim_inf = Data(:, 3).^2 ;
lim_sup = Data(:, 4).^2 ;
disp('# All restrictions loaded')

% YALMIP configuration

% Matrix G --
G = sdpvar(num_atom, num_atom);

% Restriction/conditions F

F = [G >= 0, G*ones(num_atom) == 0, lim_inf <= H*G(:) <= lim_sup];
disp('# YALMIP config sucessefuly loaded')

% Solving
% Solver options:
if solver == "sedumi"
    ops = sdpsettings('solver', 'sedumi', 'savesolveroutput', 1);
    ops.sedumi.maxiter = 15;
elseif solver == "sdpt3"
    ops = sdpsettings('solver', 'sdpt3', 'savesolveroutput', 1);
    ops.sedumi.maxiter = 15;
else
    disp('# Solver incorrect or not informed, Sedumi will be chosen !')
    ops = sdpsettings('solver', 'sedumi', 'savesolveroutput', 1);
    ops.sedumi.maxiter = 15;
end


fprintf('# Set SDP solver as: %s\n', solver);
disp('# Start Programm:')

try
    tic
    diagnostics = optimize(F,  - trace(G), ops);
    t = toc;
catch
    % Try with other solver
    msg = 'Unfortunately, a memory error probably occurred.';
    disp(msg)
    msg1 = 'Trying with SDPT3 solver...';
    disp(msg1)
    clearvars ops
    ops = sdpsettings('solver', 'sdpt3', 'savesolveroutput', 1);
    ops.sedumi.maxiter = 15;
    try
        tic
        diagnostics = optimize(F,  - trace(G), ops);
        t = toc;
    catch
        disp('# Not possible to continue due to a critical error.')
        disp('# The process will be closed.')
        return
    end
end

fprintf('\n')

fprintf('>> Solver conclusion: %s\n', yalmiperror(diagnostics.problem))

if diagnostics.problem == 9
    % Try with other solver
    msg = 'Unfortunately, sedumi found an unexpected error.';
    disp(msg)
    msg1 = 'Trying with SDPT3 solver...';
    disp(msg1)
    clearvars ops
    ops = sdpsettings('solver', 'sdpt3', 'savesolveroutput', 1);
    try
        tic
        diagnostics = optimize(F,  - trace(G), ops);
        t = toc;
    catch
        disp('# Not possible to continue due to a critical error.')
        disp('# The process will be closed.')
        return
    end
end
disp(diagnostics)
disp(diagnostics.solveroutput)
fprintf('%s\n',' ');

disp('# Get spectral decomposition of solution (solver output)')
Y = double(G);
[V, D] = svd(Y);

% non scaled solution --Solve
non_scaled = sqrt(D) * V';
non_scaled = non_scaled';

% K-scaled solution : K = 3, "projection" onto 3-dimension Euclidian Space
Q = V(:,1:3);
start_conformation = sqrt(D(1:3,1:3)) * Q';
start_conformation = start_conformation';

% Generating output files from obtained solutions
disp('# Generating output files from obtained solutions')

% first --nonscaled solution
file_dir = sprintf('%s\\relax_np.txt', outputs_path);
dlmwrite(file_dir, non_scaled, 'delimiter', ' ');

% second -- K-scaled solution
file_dir = sprintf('%s\\relax_k.txt', outputs_path);
dlmwrite(file_dir, start_conformation, 'delimiter', ' ');

% disable output diary(Log) feature;
disp('# Process completed successfully!')
disp('------------------------------------------------------')

% Save solver -- Yalmip output to file;
file_dir = sprintf('%s\\solver_varargout.txt', outputs_path);
yalmiptime = string(diagnostics.yalmiptime);
solvertime = string(diagnostics.solvertime);
generalelapsedtime = string(t);

info = string(diagnostics.info);
msg3 = sprintf('{"yalmiptime": "%s", "solvertime": "%s", "elapsed time": "%s", "info": "%s"}', yalmiptime, solvertime, generalelapsedtime, info);
diagnostics_out = msg3;
fid = fopen(file_dir, 'w+');
fprintf(fid, diagnostics_out);
fclose(fid);
end