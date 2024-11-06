function [Hessian, eigen_values, eigen_vectors] = HessYale(x, y, Dn, N, L)
% L = Ly (periodic in y)
% x = xy(1:N);
% y = xy(N+1:end);

% ID the wall particles
left_wall_list = (x<Dn/2);
right_wall_list = (x>Lx-Dn/2);

d2Ddr1dr2 = zeros(2*N, 2*N);

for n = 1:N-1
    ix_n = n;
    iy_n = n + N;
    for m = n+1:N
        ix_m = m;
        iy_m = m + N;
        Dnm = Dn(n) + Dn(m);
        
        dx = x(m) - x(n);
        % dx = dx - round(dx / L) * L;
        if abs(dx) < Dnm
            dy = y(m) - y(n);
            dy = dy - round(dy / L) * L;
            d = sqrt(dx^2 + dy^2);
            if d < Dnm
                dd = 1 - Dnm / d;
                Dnm_sq = Dnm^2;
                dx = dx / d;
                dy = dy / d;
                dx_sq = dx^2;
                dy_sq = dy^2;
                d2Ddx2 = (dx_sq + dd * dy_sq) / Dnm_sq / 2;
                d2Ddy2 = (dy_sq + dd * dx_sq) / Dnm_sq / 2;
                d2Ddxdy = dx * dy * (1 - dd) / Dnm_sq;

                d2Ddr1dr2(ix_n, ix_n) = d2Ddr1dr2(ix_n, ix_n) + d2Ddx2;
                d2Ddr1dr2(ix_n, iy_n) = d2Ddr1dr2(ix_n, iy_n) + d2Ddxdy;
                d2Ddr1dr2(iy_n, iy_n) = d2Ddr1dr2(iy_n, iy_n) + d2Ddy2;

                d2Ddr1dr2(ix_m, ix_m) = d2Ddr1dr2(ix_m, ix_m) + d2Ddx2;
                d2Ddr1dr2(ix_m, iy_m) = d2Ddr1dr2(ix_m, iy_m) + d2Ddxdy;
                d2Ddr1dr2(iy_m, iy_m) = d2Ddr1dr2(iy_m, iy_m) + d2Ddy2;

                d2Ddr1dr2(ix_n, ix_m) = d2Ddr1dr2(ix_n, ix_m) - 2 * d2Ddx2;
                d2Ddr1dr2(ix_n, iy_m) = d2Ddr1dr2(ix_n, iy_m) - d2Ddxdy;
                d2Ddr1dr2(ix_m, iy_n) = d2Ddr1dr2(ix_m, iy_n) - d2Ddxdy;
                d2Ddr1dr2(iy_n, iy_m) = d2Ddr1dr2(iy_n, iy_m) - 2 * d2Ddy2;
            end
        end
    end                    
end
%%
Hessian = d2Ddr1dr2 + d2Ddr1dr2';

[eigen_vectors, eigen_values] = eig(Hessian); % [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.