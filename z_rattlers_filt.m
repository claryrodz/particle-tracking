% Author: Clary Rodriguez-Cruz (clarydrc@seas.upenn.edu)

% Z-Zc calculation from expt. data
% input lub without ANY filters, goodenough = 2
% d = dimensionality
% max,min xyz values should be 1-2 diameters from the edges

function [z_zc,ztst,track_nonrattlers,neighbors] = z_rattlers_filt(lub,L,d,sim,constantN)

xyzrt = lub;
unqtns = unique(xyzrt(:,5));
tn = numel(unique(xyzrt(:,5)));
drops_tot = numel(unique(xyzrt(:,6)));
z_zc = zeros(tn,1);
ztst = zeros(tn,1);
Zc = 2*d;
track_nonrattlers = [];
neighbors = cell(tn,drops_tot);

% Loop over all t
for c = 1:tn
    xyzrtn = xyzrt(xyzrt(:,5)==unqtns(c),:);
    X = xyzrtn(:,2);
    Y = xyzrtn(:,1);
    Z = xyzrtn(:,3);
    R = xyzrtn(:,4);
    T = xyzrtn(:,5);
    ID = xyzrtn(:,6);
    ndrops = length(ID);
    zi = zeros(ndrops,1);

    % Loop over all N
    for i = 1:ndrops
        amat = ones(ndrops,1)*i;
        bmat = (1:ndrops)';
        pairs = [amat,bmat]; pairs(i,:) = []; % All possible pairs for N_i
        if sim==1 % Simulation data (periodic boundaries)
            rdist = sqrt((X(pairs(:,1))-X(pairs(:,2))-L*round((X(pairs(:,1))-X(pairs(:,2)))./L)).^2+...
            (Y(pairs(:,1))-Y(pairs(:,2))-L*round((Y(pairs(:,1))-Y(pairs(:,2)))./L)).^2+...
            (Z(pairs(:,1))-Z(pairs(:,2))-L*round((Z(pairs(:,1))-Z(pairs(:,2)))./L)).^2);
        elseif sim==0 % Expt. data
            rdist = sqrt((X(pairs(:,1))-X(pairs(:,2))).^2+...
            (Y(pairs(:,1))-Y(pairs(:,2))).^2+...
            (Z(pairs(:,1))-Z(pairs(:,2))).^2);
        end
        radsum = R(pairs(:,1))+R(pairs(:,2));
        fidx = find(rdist<=radsum); % Find bubbles in contact with N_i
        zi(i) = numel(fidx);
        if constantN == 1
            neighbors{c,i} = ID(pairs(fidx,2));
        end
    end
    
    z_zc(c) = mean(zi(zi>=Zc))-Zc;
    ztst(c) = sum(zi<Zc); % track # of rattlers at each t
    find_rattlers = zi<Zc; 
    rattlers = ID(find_rattlers); % get rattler IDs, Z < Z_tol
    [~,~,keep] = setxor(rattlers,ID); % get IDs that are not rattlers
    track_nonrattlers = [track_nonrattlers;xyzrtn(keep,:)]; % keep xyzrt of non-rattlers
end

end