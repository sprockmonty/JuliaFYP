function [a,e,i,rightAsc,omega] = earth2orbital(r,v)
%earth2orbital Transforms an earth centred reference frame into orbital
%parameters
%   Inputs:
%       r - an array containing the position coordinates in the earth
%       centred frame r = [x,y,z]
%       v - an array containing the velocity vector in the earth centred
%       frame v = [vx, vy, vz]
%   Outputs:
%       a        - semi-major axis
%       e        - eccentricity
%       i        - inclination (in degrees)
%       rightAsc - right ascension of the ascending node (in degrees)
%       omega    - argument of perigee (in degrees)

% No error checking occurs to speed up code assuming user has imput correct
% array lengths
    mu = 3.986012e14;
    h = cross(r,v);
    ecc = cross(v,h) / mu - r / norm(r);
    n = cross([0,0,1]',h) ;
    i = acos(h(3) / norm(h));
    e = norm(ecc);
    if n(2) >= 0
        rightAsc = acos(n(1) / norm(n));
    else
        rightAsc =  2*pi - acos(n(1) / norm(n));
    end
    if ecc(3) >= 0
        omega = acos(dot(n,ecc) / (norm(n)*norm(ecc)));
    else
        omega = 2*pi - acos(dot(n,ecc) / (norm(n)*norm(ecc)));
    end
    rightAsc = rad2deg(rightAsc);
    omega = rad2deg(omega);
    i = rad2deg(i);
    a = 1 / (2/norm(r) - norm(v)^2 / mu);
    if isnan(i)
        [a,e,i,rightAsc,omega] = earth2orbital(r,v+eps); % add small eps if singularity to remove nan
    end
end

