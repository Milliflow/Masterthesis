q = 1; % twist
r = pi/10; % should translate to the discretized model with 20 nodes
for i = 1:1000 % combined with the next code line this translates to the intervall [0,10]
    x = i/100;
    k = 1;
    a(i) = C1(q,k,r)+2*x*(W(r,q,q-k)+W(r,q,q+k)-2*W(r,q,q)); % stores the highest eigenvalue for the q-twist
    b(i) = C1(q,k,r)+2*x*(W(r,0,-k)+W(r,0,k)-2*W(r,0,0)); % stores the highest eigenvalue for the 0-twist
    for k = 2:5 % iterates through several k's to find the highest eigenvalue
        achallenger = C1(q,k,r)+2*x*(W(r,q,q-k)+W(r,q,q+k)-2*W(r,q,q));
        bchallenger = C1(q,k,r)+2*x*(W(r,0,-k)+W(r,0,k)-2*W(r,0,0));
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
xlabel("\lambda");
ylabel("quotient");
title("Higher Order 2 Characteristic Return Time quotient");

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function v = W(r,q,k)
    if q==0 && k == 0
        v = 3*r^2;
    elseif q == 0 && k~=0
        v = 3*pi*r*V(r/2,k)*cos(k*r/2)+2*sin(k*r/2)*(2*sin(k*r/2)-k*r*cos(k*r/2))/(k^2);
    elseif k == 0 && q~=0
        v = 3*pi*r*V(r/2,q)*cos(q*r/2)+2*sin(q*r/2)*(2*sin(q*r/2)-q*r*cos(q*r/2))/(q^2);
    else
        v = pi/k*V(r/2,q)*(sin((k-q/2)*r)+sin((k+q/2)*r));
        v = v+pi/k*V(r/2,q-k)*sin((k+q)*r/2)+pi/q*V(r/2,q+k)*sin((q-k)*r/2);
    end
    v = v/pi;
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
