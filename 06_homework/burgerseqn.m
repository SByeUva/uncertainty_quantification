%%
a=-1; b=1; alpha=1; beta=-1; % left and right boundary points and their values
n=100; % number of spatial grid points for discretizing the PDE
x=linspace(a,b,n); % create spatial grid for solving PDE
M=8;
tt=linspace(0,M,M+1);

k=rescale(cos(tt*pi/M),0.05,0.35);

textprintout = false;
ws = [];
for i=1:M+1
    kappa = k(i);

    w = bsv ( a, b, alpha, beta, kappa, n, textprintout );
    ws = [ws, w];
end
figure(1); plot(k,ws(60,:),'ob'); grid;


%%
dz=0.001;
z=linspace(0.05,0.35,n);
zx=size(z);
fa=zeros(n, zx(2));
for h=1:n
    for j=1:length(z)
        for i=1:M+1
            qi=lagrangebasisfunc(z(j),k,i);
            fa(h,j)=fa(h,j)+ws(h,i)*qi;
        end
    end
end

figure(2); plot(x,fa,'-r'); grid;
%%
function y=lagrangebasisfunc(x,xi,i)

    % evaluate i-th Lagrange basis function at x, given nodes xi
    % y = l_i(x) = \product(j \neq i) (x-x_j)/(x_i-x_j)
    % according to the definition in Xiu

    y=1;

    for j=1:length(xi)
        if i~=j
            y = y*(x-xi(j))/(xi(i)-xi(j));
        end
    end

end

