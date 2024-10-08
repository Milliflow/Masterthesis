% This codes gives a polar representation by using a bisection algorithm; it's quite fast.
n = 1000; % resolution
q = [1 1]; % twist
r = zeros(1,n); % stores the radius values
positive = false;
for phi = 1:n
    r(phi) = 0.69*pi; % start value, choosen near the stable change at 0.68*pi in the 1D-case
    for j = 1:30 % refinement steps
        rgeo = r(phi)*[sin(pi*phi/(2*n)) cos(pi*phi/(2*n))]; % cartesian representation
        positive = false;
        for a = -10:10 % This and the following line give the range of eigenfunctions which are analyzed
            for b = 0:10
                k = [a b];
                c1 = C1(q,k,rgeo);
                if c1 > 0
                    positive = true;
                end
            end
        end
        % bisection: if positive "go down", if negative "go up"
        if positive
            r(phi) = r(phi)-pi*2^(-j-1);
        else
            r(phi) = r(phi)+pi*2^(-j-1);
        end
    end
end
r = r/(2*pi); % choosen for comparison with the values in the 1D-case in "The size of the sync basin" (Wiley)
plot(r);
xticklabels({'0','','','','','\pi/4','','','','','\pi/2'});
title('Plot of boundary of the stable region in polar coordinates');
xlabel('\phi \in [0,\pi/2]');
ylabel('r/(2\pi)');

function c1 = C1(q,k,r)
% Calculates the eigenvalue
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function v = V(r,k)
% V Calculates the fourier coefficient for the given values
    v = 2;
    for i=1:length(k)
        if k(i) == 0
            v_i = r(i)/pi;
        else
            v_i = sin(k(i)*r(i))/(pi*k(i));
        end
        v = v*v_i;
    end
end
