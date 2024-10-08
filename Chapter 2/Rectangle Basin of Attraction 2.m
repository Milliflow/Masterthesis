% What happens here? This is a n1 x n2 discretization of the 2D-Kuramoto model, with an interaction domain of lambda x 2pi; the basin of attraction
% is calculated for each twist, stored in the matrix c and displayed during runtime. As predicted by theory, no twists in y-direction occur. 
% This program is unfortunately very slow. An explicit euler approximation is used.
n1 = 15; % Nodes in x-direction
n2 = 15; % Nodes in y-direction; if set to 1 the normal cyclic Kuramoto model should emerge
x = zeros(n1,n2);
dt = 0.005; % size of timestep
T = 30000; % Number of timesteps: full runtime = dt*T
Test = 10;
c = zeros(5,5,Test); % stores the 2D-twist in the matrix where an (q,p)-twist is displayed at position (q+1,p+1) in the matrix
for lambda = 1:Test % for t = 0 all rings will be independent -> random twists -> programm not suited for this
    for Stat = 1:100 % Number of runs
        x = 2*pi*(rand(n1,n2)-0.5); % Initialisation
        x = dyKuramoto(x,lambda,dt,T); % Simulation
        q = zeros([0,n2]); % twist in x-direction (may vary)
        p = zeros([0,n1]); % twist in y-direction (should be 0)
        for j = 1:n2
            q(j) = mod(pi+x(1,j)-x(n1,j),2*pi)-pi;
            for i = 1:n1-1
                q(j) = q(j)+mod(pi+x(i+1,j)-x(i,j),2*pi)-pi; % Detection
            end
        end
        for j = 1:n1
            p(j) = mod(pi+x(j,1)-x(j,n2),2*pi)-pi;
            for i = 1:n2-1
                p(j) = p(j)+mod(pi+x(j,i+1)-x(j,i),2*pi)-pi; % Detection
            end
        end
        q = q./(2*pi);
        p = p./(2*pi);
        for j1 = 0:4
            qtrue = true;
            for i = 1:n2
                if abs(abs(q(i))-j1)>0.4
                    qtrue = false;
                end
            end
            if qtrue
                for j2 = 0:4
                    ptrue = true;
                    for i = 1:n1
                        if abs(abs(p(i))-j2)>0.4
                            ptrue = false;
                        end
                    end
                    if ptrue
                        c(j1+1,j2+1,lambda) = c(j1+1,j2+1,lambda)+1;
                    end
                end
            end
        end
        disp(c(:,:,lambda)./Stat); % Display of current basin of attraction calculations
    end
end
disp(c);

function y = dyKuramoto(x,lambda,dt,T)
    sz = size(x);
    n = sz(1);
    k = sz(2);
    yalt = [x; x];
    yneu = yalt;
    for t = 1:T % Loop over time
        for i = 1:n % Loop over nodes in x-direction
            for j = 1:lambda % Loop over interaction domain in x-direction
                for l = 1:k % Loop over nodes in y-direction
                    for m = 1:k % Loop over interaction domain in y-direction
                        yneu(i,l) = yneu(i,l)+dt*(sin(yalt(i+j,m)-yalt(i,l))+sin(yalt(i+n-j,m)-yalt(i,l)));
                    end
                end
            end                                                 
        end
        for i = 1:k
            yalt(:,i) = [yneu(1:n,i); yneu(1:n,i)];
        end
    end
    y = yneu(1:n,1);
    for i = 2:k
        y = [y yneu(1:n,i)];
    end
end
