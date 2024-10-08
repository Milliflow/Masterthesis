% Basin of Attraction of model with rectangle interaction domain
% where the rectangle is 3x3 nodes large; the ODE is solved by
% an explicit euler approximation

n1 = 20; % resolution of x-axis
n2 = 20; % resolution of y-axis
x = zeros(n1,n2);
dt = 0.005; % step-width
T = 300000; % number of steps
c = zeros(5,5); % stores the end configuration
for Stat = 1:100 % Number of runs
    x = 2*pi*(rand(n1,n2)-0.5); % random start configuration
    x = dyKuramoto(x,dt,T); % Simulation
    xend = x;
    q = zeros([0,n2]); % checks twists in x-direction
    p = zeros([0,n1]); % checks twists in y-direction
    for j = 1:n2
        q(j) = mod(pi+x(1,j)-x(n1,j),2*pi)-pi;
        for i = 1:n1-1
            q(j) = q(j)+mod(pi+x(i+1,j)-x(i,j),2*pi)-pi;%Detektion
        end
    end
    for j = 1:n1
        p(j) = mod(pi+x(j,1)-x(j,n2),2*pi)-pi;
        for i = 1:n2-1
            p(j) = p(j)+mod(pi+x(j,i+1)-x(j,i),2*pi)-pi;%Detektion
        end
    end
    q = q./(2*pi);
    p = p./(2*pi);
    for j1 = 0:4
        qtrue = true; % checks if all twists in x-direction are equal
        for i = 1:n2
            if abs(abs(q(i))-j1)>0.4
                qtrue = false;
            end
        end
        if qtrue
            for j2 = 0:4
                ptrue = true; % checks if all twists in y-direction are equal
                for i = 1:n1
                    if abs(abs(p(i))-j2)>0.4
                        ptrue = false;           
                    end
                end
                if ptrue
                    c(j1+1,j2+1) = c(j1+1,j2+1)+1;
                end
            end
        end
    end
    disp(c./Stat);
end
disp(c);

function y = dyKuramoto(x,dt,T)
    sz = size(x);
    n1 = sz(1);
    n2 = sz(2);
    yalt = [x x; x x];
    yneu = yalt;
    for t = 1:T
        for i = 1:n1
            for j = 1:n2
                yneu(i,j) = yneu(i,j)+dt*(sin(yalt(i+1,j)-yalt(i,j))+sin(yalt(i+n1-1,j)-yalt(i,j))+sin(yalt(i,j+1)-yalt(i,j))+sin(yalt(i,j+n1-1)-yalt(i,j))+sin(yalt(i+1,j+1)-yalt(i,j))+sin(yalt(i+n1-1,j+n1-1)-yalt(i,j))+sin(yalt(i+n1-1,j+1)-yalt(i,j))+sin(yalt(i+1,j+n1-1)-yalt(i,j)));
            end                                                
        end
    yalt = [yneu(1:n1,1:n2) yneu(1:n1,1:n2); yneu(1:n1,1:n2) yneu(1:n1,1:n2)];
    end
    y = yneu(1:n1,1:n2);
end
