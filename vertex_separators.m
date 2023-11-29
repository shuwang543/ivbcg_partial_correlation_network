function [separator_list, level_list, parent_list] = vertex_separators(adjmat, x_a, x_b, maxk)

%This MATLAB function implements the algorithm for enumerating minimal
%separator sets of a graph, as described in:
%"Efficient enumeration of all minimal separators in a graph 
% - Hong Shen & Weifa Liang - Theoretical Computer Science 1997"

    % graphmodel input as a binary adjacency matrix, x_a, x_b the indices
    % of the two points a,b being separated.
    % maxk specifies the maximum number of levels to search, to save time
    
    % separator_list: a logical matrix of width equal to number of nodes, and height
    % equal to the total number of minimal separators for a,b
    % level_list: vector indicating the level of a separator
    % parent_list: index of parent separator

    n = size(adjmat,1);
    nodenames = cellstr(string(1:n)');
    
    if maxk>n-3
        maxk=n-3;
    end
    
    %1) compute C_b connected component of b of graph G[V-N(a)]
    %2) compute I(N(a)) the isolateset
        N_a = adjmat(x_a,:); %logical vector of x_a's neighbors
        [isolate_N_a,C_b_names] = isolateset(adjmat,N_a,x_b,nodenames);
    %3) define level 0 of the tree, L_0 = {N(a)-I(N(a))}, set k=0
        separator_list = N_a & ~isolate_N_a; %initialize
        level_list = 0;
        parent_list = 0;
        k = 0;
    %4) while (k \le n-3) AND C_b is not empty
    while (k <= maxk) && (numel(C_b_names)>0)
        tic
        disp(k)
    %   for each separator S \in L_k, 
        L_k = separator_list(level_list==k,:); %each row is a separator
        if k==0
            L_k_parents = separator_list(1,:)~=separator_list(1,:);
        else
            L_k_parents = separator_list(parent_list(level_list==k),:);
        end
        for i = 1:size(L_k,1)
    %       for each node x in separator S not adjacent to b
            parent_separator = L_k_parents(i,:);
            x_S_b = find(L_k(i,:) & ~adjmat(x_b,:));
            for j = 1:numel(x_S_b)
                proceed_id = proceedset(adjmat,parent_separator,x_S_b(j));
                X = L_k(i,:) | proceed_id;
    %           compute C_b of G[V-(S\union N^+(x))]
    %           if C_b is not empty
                    [isolate_id,C_b_names] = isolateset(adjmat,X,x_b,nodenames);
                    %               compute isolateset I(S\union N^+(x))
                    if numel(C_b_names)==0
                        continue
                    else
    %               S' = S\unionN^+(x)) - I(S\union N^+(x)); generating a
    %               potential new separator for the level k+1
                    sep_candid = (X & ~isolate_id);
    %               if S' is not an existing separator in previous levels,
    %               then add S' to the list in L_k+1
                        if any(all(sep_candid==separator_list,2))==0
                        separator_list = [separator_list; sep_candid];
                        level_list = [level_list; k+1];
                        parent_list = [parent_list; find(all(separator_list==L_k(i,:),2))];
                        else
                            continue
                        end
                    end
                    
            end
        
        end
        k = k+1;
        toc
    end
end

function [isolate_id,C_b_names] = isolateset(adjmat,X,x_b,nodenames) %compute the isolated set of a set of nodes
    %give adjmat as adjacency matrix
    %give separator as logical vector indicating nodes
    %give X as a logical vector indicating nodes
    %x_a and x_b are indices
    G_X = graph(adjmat(~X,~X),nodenames(~X));
    bins = conncomp(G_X);
    C_b = bins==bins(strcmp(G_X.Nodes.Name,nodenames{x_b}));
    C_b_names = str2double(string(G_X.Nodes.Name(C_b)));
    
    non_adj_nodes = all(adjmat(:,C_b_names)==0,2);
    isolate_id = non_adj_nodes' & X;
end

function proceed_id = proceedset(adjmat, S_parent, x) %compute the proceeding set (N^+(x)) of a node x
%S_parent is a logical vector indicating nodes in a parent Separator

preced_x = S_parent & adjmat(x,:);
proceed_id = adjmat(x,:) & ~preced_x;

end

