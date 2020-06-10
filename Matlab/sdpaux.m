% This function reads the avaliable proteins for test,
% and the executes the SDP relax algorithm for them.

% Open the file containing the protein names and paths;
fid = fopen('proteins.txt', 'r+');
proteins = textscan(fid, '%s %s %s', 'delimiter', ',');
fclose(fid);
disp('# proteins.txt read complete!')

% style: (Nodes) (Nodes-path) (Tests-path)
Nodes = proteins{1, 1};
Nodes_path = proteins{1, 2};
Tests_path = proteins{1, 3};


n = length(Nodes);
fprintf('# Total number of nodes read: %d \n', n)
for i = 1 : n
    spdprotein(Nodes{i}, Nodes_path{i}, Tests_path{i}, 'sedumi', 10)
end
disp('# The SDP process was successfully completed !')