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
%       i        - inclination
%       rightAsc - right ascension of the ascending node
%       omega    - argument of perigee
    mu = 3.986012e14;
    h = cross(r,v);
    n = cross([0,0,1], h);
    ecc = ((norm(v)^2 - mu / norm(r)) * r - (dot(r,v)) * v ) / mu;
    e = norm(ecc);
    E = norm(v)^2 /2 - mu / norm(r);
    if abs(e-1) > eps
        a = -mu/(2*E);
        p = a*(1-e^2);
    else
        a = inf;
        p = norm(h)^2 / mu;
    end
    i = acosd(h(3) / norm(h));
    rightAsc = acosd(n(1) / norm(n));
    omega = acosd(dot(n,ecc) / (norm(n)*e));
    if n(2) < 0
        rightAsc = 360 - rightAsc;
    end
    
    if ecc(3) < 0
        omega = 360 - omega;
    end
    
end

