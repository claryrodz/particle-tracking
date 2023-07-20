% Author: Clary Rodriguez-Cruz (clarydrc@seas.upenn.edu)

% Function calculates the mean squared displacement of all particles using
% the overcount method (counting all displacements at each dt)
% track = [x y z r t id]
% tau_spacing = 0 for linear,  = m for log-spaced where m is the # of points
% L = box size (accounts for periodic boundaries)

function [lagtime,msd] = msd_overcount(track,tau_spacing,L)

X = track(:,2);
Y = track(:,1);
Z = track(:,3);
T = track(:,5);
ID = track(:,6);
ts = unique(T);
dT = ts(2)-ts(1);
tn = numel(unique(T));
unqdropsID = unique(ID);
numdrops = numel(unqdropsID);
tauspaced = 1:((round(tn/dT))-2);
if tau_spacing ~= 0
    tauspaced = unique(round(logspace(log10(1),log10(tn/dT-3),tau_spacing)));
end
tnL=length(tauspaced);

msd = zeros(tnL,1);
delxmean = zeros(tnL,1);
delymean = zeros(tnL,1);
delzmean = zeros(tnL,1);
delx2mean = zeros(tnL,1);
dely2mean = zeros(tnL,1);
delz2mean = zeros(tnL,1);

for tau = 1:length(tauspaced)
    sumdelx = 0; sumdely = 0; sumdelz = 0; sumdelx2 = 0; sumdely2 = 0; sumdelz2 = 0; N = 0; Neff = 0; tot=[]; radsum=0;
    for j = 1:numdrops
        drop_idx = find(ID == unqdropsID(j));  
        if length(drop_idx) > tauspaced(tau)+1
        x_0 = X(drop_idx);
        y_0 = Y(drop_idx);
        z_0 = Z(drop_idx);
        t_0 = T(drop_idx);
        delx = x_0(tauspaced(tau)+1:end)-x_0(1:end-tauspaced(tau))-L*round((x_0(tauspaced(tau)+1:end)-x_0(1:end-tauspaced(tau)))./L);
        dely = y_0(tauspaced(tau)+1:end)-y_0(1:end-tauspaced(tau))-L*round((y_0(tauspaced(tau)+1:end)-y_0(1:end-tauspaced(tau)))./L);
        delz = z_0(tauspaced(tau)+1:end)-z_0(1:end-tauspaced(tau))-L*round((z_0(tauspaced(tau)+1:end)-z_0(1:end-tauspaced(tau)))./L);
        delx2 = delx.^2;
        dely2 = dely.^2; 
        delz2 = delz.^2;
        delt = t_0(tauspaced(tau)+1:end)-t_0(1:end-tauspaced(tau));
        tfind = delt == tauspaced(tau)*dT;
        delx=delx(tfind);dely=dely(tfind);delz=delz(tfind);sumr=sumr(tfind);
        delx2=delx2(tfind);dely2=dely2(tfind);delz2=delz2(tfind);
        sumdelx = sumdelx+sum(delx);
        sumdelx2 = sumdelx2+sum(delx2);
        sumdely = sumdely+sum(dely);
        sumdely2 = sumdely2+sum(dely2);
        sumdelz = sumdelz+sum(delz);
        sumdelz2 = sumdelz2+sum(delz2);
        N = N+sum(tfind);
        end
    end
    delxmean(tau) = sumdelx/N;
    delx2mean(tau) = sumdelx2/N;
    delymean(tau) = sumdely/N;
    dely2mean(tau) = sumdely2/N;
    delzmean(tau) = sumdelz/N;
    delz2mean(tau) = sumdelz2/N;
    msd(tau) = (delx2mean(tau)-delxmean(tau)^2)+(dely2mean(tau)-delymean(tau)^2)+(delz2mean(tau)-delzmean(tau)^2);
end

lagtime = tauspaced'*dT;

end