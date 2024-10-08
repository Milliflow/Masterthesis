m = 1000;
% spherecoeff has a range of 1:m*1:21*1:21 which corresponds to
% eigenfunctions from -10 to 10

% First step: Calculating the fourier coefficients
for a=1:m
    for i=1:21
        for j=1:21
            r = pi*a/m;
            k = [i j]-[11 11]; % (0 0) in analysis = (11 11) in the array
            if k ~ [0 0];
                spherecoeff(a,i,j) = integral(@(x) sin(r*k(1)*sin(x)+r*k(2)*cos(x))./(k(1)*sin(x)+k(2)*cos(x)),0,2*pi);
            else
                spherecoeff(a,i,j) = pi*r*r;
            end
        end
    end
end

% Second step: Calculating the stable zone and the bifurcation at the boundary
A = zeros(m,m);
B = A;
n = 4;
q = [1 1]; % twist; please be beware of out of range errors if you set this value higher up
qfA = q + [11 11]; % =qforArray; gives the twist in relation to the array (remember that the eigenfunction  [0 0] equals [11 11] there)
for i=1:m
    for j=1:i-1
        for k1=-n:n
            for k2=0:n                
                k = [k1 k2];
                r = [i j];
                % Please notice that in the following for z, c1, c2, c3 and c5 we not only have to calculate the values for a circle, but also subtract the values for the
                % smaller circle, hence exactly double as many terms as usual occur, the first ones with r(1), the second ones with r(2)
                z = spherecoeff(r(1),qfA(1)+k(1),qfA(2)+k(2))+spherecoeff(r(1),qfA(1)-k(1),qfA(2)-k(2))-2*spherecoeff(r(1),qfA(1),qfA(2))-spherecoeff(r(2),qfA(1)+k(1),qfA(2)+k(2))-spherecoeff(r(2),qfA(1)-k(1),qfA(2)-k(2))+2*spherecoeff(r(2),qfA(1),qfA(2));
                if z>0
                    A(i,j) = 1;
                    c1 = (spherecoeff(r(1),qfA(1)+2*k(1),qfA(2)+2*k(2))+spherecoeff(r(1),qfA(1)-2*k(1),qfA(2)-2*k(2))-2*spherecoeff(r(1),qfA(1),qfA(2)))/(2*pi)-(spherecoeff(r(2),qfA(1)+2*k(1),qfA(2)+2*k(2))+spherecoeff(r(2),qfA(1)-2*k(1),qfA(2)-2*k(2))-2*spherecoeff(r(2),qfA(1),qfA(2)))/(2*pi);
                    c2 = (-spherecoeff(r(1),qfA(1)-2*k(1),qfA(2)-2*k(2))+2*spherecoeff(r(1),qfA(1)-k(1),qfA(2)-k(2))-2*spherecoeff(r(1),qfA(1)+k(1),qfA(2)+k(2))+spherecoeff(r(1),qfA(1)+2*k(1),qfA(2)+2*k(2)))/(4*pi)-(-spherecoeff(r(2),qfA(1)-2*k(1),qfA(2)-2*k(2))+2*spherecoeff(r(2),qfA(1)-k(1),qfA(2)-k(2))-2*spherecoeff(r(2),qfA(1)+k(1),qfA(2)+k(2))+spherecoeff(r(2),qfA(1)+2*k(1),qfA(2)+2*k(2)))/(4*pi);
                    c3 = (-spherecoeff(r(1),qfA(1)-2*k(1),qfA(2)-2*k(2))+spherecoeff(r(1),qfA(1)-2*k(1)+k(1),qfA(2)-2*k(2)+k(2))+spherecoeff(r(1),qfA(1)-k(1),qfA(2)-k(2))+spherecoeff(r(1),qfA(1)+2*k(1),qfA(2)+2*k(2))-spherecoeff(r(1),qfA(1)+2*k(1)-k(1),qfA(2)+2*k(2)-k(2))-spherecoeff(r(1),qfA(1)+k(1),qfA(2)+k(2)))/(4*pi)-(-spherecoeff(r(2),qfA(1)-2*k(1),qfA(2)-2*k(2))+spherecoeff(r(2),qfA(1)-2*k(1)+k(1),qfA(2)-2*k(2)+k(2))+spherecoeff(r(2),qfA(1)-k(1),qfA(2)-k(2))+spherecoeff(r(2),qfA(1)+2*k(1),qfA(2)+2*k(2))-spherecoeff(r(2),qfA(1)+2*k(1)-k(1),qfA(2)+2*k(2)-k(2))-spherecoeff(r(2),qfA(1)+k(1),qfA(2)+k(2)))/(4*pi);
                    c5 = (spherecoeff(r(1),qfA(1)-2*k(1),qfA(2)-2*k(2))-4*spherecoeff(r(1),qfA(1)-k(1),qfA(2)-k(2))+6*spherecoeff(r(1),qfA(1),qfA(2))-4*spherecoeff(r(1),qfA(1)+k(1),qfA(2)+k(2))+spherecoeff(r(1),qfA(1)+2*k(1),qfA(2)+2*k(2)))/(8*pi)-(spherecoeff(r(2),qfA(1)-2*k(1),qfA(2)-2*k(2))-4*spherecoeff(r(2),qfA(1)-k(1),qfA(2)-k(2))+6*spherecoeff(r(2),qfA(1),qfA(2))-4*spherecoeff(r(2),qfA(1)+k(1),qfA(2)+k(2))+spherecoeff(r(2),qfA(1)+2*k(1),qfA(2)+2*k(2)))/(8*pi);
                   B(i,j) = c5-c2*c3/c1;
                end
            end
        end
    end
end

% post-processing
for i = 1:m
    for j = 1:m
        if(i<=j)
            A(i,j)=1; % Gives unattainable domain the value 1
        end
    end
end

for i = 1:m
    for j = 1:i-1
        if A(i,j)~0
            if i>5 && j>5
                if A(i-5,j-5)==0 && B(i,j)<0
                    B(i,j)=2;
                elseif A(i-5,j-5)==0
                    B(i,j)=4;
                else
                    B(i,j)=0;
                end
            elseif j<6 && i>6
                if A(i-10,j)==0 && B(i,j)<0
                    B(i,j)=2;
                elseif A(i-10,j)==0
                    B(i,j)=4;
                else
                    B(i,j)=0;
                end
            end
        end
    end
end
for i = 1:m
    for j = 1:i-1 % loop over the attendable region
        if B(i,j)==2
            B(i,j)=2;
        elseif B(i,j)==4
            B(i,j)=4;
        elseif A(i,j)==0
            B(i,j)=1;
        else
            B(i,j)=3;
        end
    end
end
figure('pos',[100 100 1000 800]);
ax = subplot(1,1,1);
imagesc([0,pi],[0,pi],B);
xlabel("r_1");
ylabel("r_2");
ax.YDir = "normal";
title('Bifurcation Portrait')
cmap = zeros(10, 4);
cmap = [0, 0, 0; ...
  0, 1, 0; ...   % Blue for 1
  0, 1, 1; ...       % Black for 2
  1, 0, 0];       % Green for 3
C = cmap;
L = size(C,1);
Gs = round(interp1(linspace(min(A(:)),max(A(:)),L),1:L,A));
H = reshape(C(As,:),[size(As) 3]);
subplot(1,2,2);
