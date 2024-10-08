% blue is the stable, yellow the instable region; to get the colours from
% the thesis for q = [1 1] and q = [2 3] type colormap(cmap); please note that
% this colormap is not yet suited for the occurence of supercritical
% bifurcations, in the case other q are choosen

q = [1 1]; % twist
m = 1000; % resolution
A = zeros(m,m); % stores eigenvalues above 0
n = 10; % number of analyzed eigenfunctions
B = ones(m,m); % stores bifurcation values
for i = 1:m
    for j = 1:m
        for k1 = 0:n
            for k2 = -n:n
                k = [k1 k2];
                r = [pi*i/m pi*j/m];
                z = C1(q,k,r); % respective eigenvalue
                if z>0 && A(i,j) == 0
                    A(i,j) = 1;
                    B(i,j) = C5(q,k,r)-C2(q,k,r)*C3(q,2*k,k,r)/C1(q,2*k,r); % bifurcation formula
                end
            end
        end
    end
end

% post-processing
for i = 1:m
    for j = 1:m
        if A(i,j)~0
            if i>5 && j>5
                if A(i-5,j-5)==0 && B(i,j)>0
                    B(i,j) = 2;
                elseif A(i-5,j-5)==0
                    B(i,j) = 4;
                else
                    B(i,j) = 3;
                end
            elseif i<6 && j>5
                if A(i,j-5)==0 && B(i,j)>0
                    B(i,j) = 2;
                elseif A(i,j-5)==0
                    B(i,j) = 4;
                else
                    B(i,j) = 3;
                end
            else
                if A(i-5,j)==0 && B(i,j)>0
                    B(i,j) = 2;
                elseif A(i-5,j)==0
                    B(i,j) = 4;
                else
                    B(i,j) = 3;
                end
            end
        end
    end
end
A = A';
B = B';
cmap = zeros(10, 3);
cmap = [0, 1, 0; ...   % Green for 1 -> stable region
  1, 1, 0; ...       % Yellow for 2 -> bifurcation region
  1, 0, 0];       % Red for 3 -> instable region
figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[0,pi],B);
xlabel("r_1");
ylabel("r_2");
ax.YDir = "normal";
title('Bifurcation Portrait');
C = cmap();
L = size(C,1);
Gs = round(interp1(linspace(min(B(:)),max(B(:)),L),1:L,B));
H = reshape(C(Bs,:),[size(Bs) 3]);
subplot(1,2,2);
imshow(B, cmap);

function c1 = C1(q,k,r)
    c1 = (2*pi)^length(r)*[V(r,q+k)+V(r,q-k)-2*V(r,q)]/4;
end

function c2 = C2(q,k,r)
    c2 = (2*pi)^length(r)*(-V(r,q-2*k)+2*V(r,q-k)-2*V(r,q+k)+V(r,q+2*k))/8;
end

function c3 = C3(q,k,j,r)
    c3 = (2*pi)^length(r)*(-V(r,q-k)+V(r,q-k+j)+V(r,q-j)+V(r,q+k)-V(r,q+k-j)-V(r,q+j))/8;
end

function c5 = C5(q,k,r)
    c5 = (2*pi)^length(r)*(V(r,q-2*k)-4*V(r,q-k)+6*V(r,q)-4*V(r,q+k)+V(r,q+2*k))/16;
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
