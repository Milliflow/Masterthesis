q = 1; % twist
r = [pi/3 pi/3]; % the first number is r_{1,1}, the second one r_{2,2}
a = 0; % stores the highest eigenvalue for the q-twist
b = 0; % stores the highest eigenvalue for the 0-twist
for i = 1:400 % combined with the next code line this translates to the intervall [0,4]
    x = i/100;
    k = 1;
    a(i) = max(real(Tra(q,k,x,r)/2+sqrt(Tra(q,k,x,r)^2/4-Dete(q,k,x,r))),real(Tra(q,k,x,r)/2-sqrt(Tra(q,k,x,r)^2/4-Dete(q,k,x,r))));
    b(i) = max(real(Tra(0,k,x,r)/2+sqrt(Tra(0,k,x,r)^2/4-Dete(0,k,x,r))),real(Tra(0,k,x,r)/2-sqrt(Tra(0,k,x,r)^2/4-Dete(0,k,x,r))));
    for k = 2:5 % iterates through several k's to find the highest eigenvalue
        achallenger = max(real(Tra(q,k,x,r)/2+sqrt(Tra(q,k,x,r)^2/4-Dete(q,k,x,r))),real(Tra(q,k,x,r)/2-sqrt(Tra(q,k,x,r)^2/4-Dete(q,k,x,r))));
        bchallenger = max(real(Tra(0,k,x,r)/2+sqrt(Tra(0,k,x,r)^2/4-Dete(0,k,x,r))),real(Tra(0,k,x,r)/2-sqrt(Tra(0,k,x,r)^2/4-Dete(0,k,x,r))));
        if a(i) < achallenger
            a(i) = achallenger;
        end
        if b(i) < bchallenger
            b(i) = bchallenger;
        end
    end
end
plot(a./b); % plots the quotient; for the other graphs, plot -1./a and -1./b
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt/100)
xlabel("r_{1,2}");
ylabel("quotient");
title("Double Ring Network Characteristic Return Time quotient");


function tr = Tra(q,k,x,r)
% Calculates the trace
    tr = C1(q,k,r(1))+C1(q,k,r(2))-2*V(x,q);
end

function det = Dete(q,k,x,r)
% Calculates the determinant
    det = (C1(q,k,r(1))-V(x,q))*(C1(q,k,r(2))-V(x,q))-((V(x,q-k)+V(x,q+k))^2)/4;
end

function c1 = C1(q,k,r)
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
