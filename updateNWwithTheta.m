function [nw,dAdt,dBdt] = updateNWwithTheta(nw,theta,u0)

theta = theta(:).';

num_lines = nw.num_lines;
non_zib_index_array = nw.non_zib_index_array;
elems_comb = nw.elems_comb;
downstream_matrix = nw.downstream_matrix;
downstream_matrix_p = nw.downstream_matrix_p;

r = theta(1:num_lines);
x = theta(num_lines+1:2*num_lines);

ne = length(elems_comb);
nl = num_lines;

nw.i1 = [nw.comb_idx1(:);nl+nw.comb_idx1(:)];
nw.i2 = [nw.comb_idx2(:);nl+nw.comb_idx2(:)];

big_r_mat_coeff = nw.big_r_mat_coeff;
% big_r_mat = zeros(ne,nl);
% big_r_mat(:) = big_r_mat_coeff.' * r(:);
% big_x_mat = zeros(ne,nl);
% big_x_mat(:) = big_r_mat_coeff.' * x(:);
% %nw.big_r_mat = big_r_mat;
% %nw.big_x_mat = big_x_mat;
br = big_r_mat_coeff(:,1:ne).' * r(:);
bx = big_r_mat_coeff(:,1:ne).' * x(:);

V_node_idx = nw.V_node_idx;
nn = length(non_zib_index_array);
nv = length(V_node_idx);

sens_mat_r_coeff = nw.sens_mat_r_coeff;
R_sens_mat = zeros(nn,nv);
R_sens_mat(:) = sens_mat_r_coeff.' * r(:);
X_sens_mat = zeros(nn,nv);
X_sens_mat(:) = sens_mat_r_coeff.' * x(:);
%nw.R_sens_mat = R_sens_mat;
%nw.X_sens_mat = X_sens_mat;

z_mat_coeff = nw.z_mat_coeff;
z = r.^2 + x.^2;
z_mat_loss = zeros(ne,nv);
z_mat_loss(:) = z_mat_coeff.' * z(:);
%nw.z_mat_loss = z_mat_loss;

rhat_inner_coeff = nw.rhat_inner_coeff;
r_mat_inner = zeros(ne,nl);
r_mat_inner(:) = rhat_inner_coeff.' * r(:);
x_mat_inner = zeros(ne,nl);
x_mat_inner(:) = rhat_inner_coeff.' * x(:);
%nw.r_mat_inner = r_mat_inner;
%nw.x_mat_inner = x_mat_inner;

downstream_matrix_of_V = downstream_matrix.';
r_mat_out = r .* downstream_matrix_of_V;
x_mat_out = x .* downstream_matrix_of_V;
r_mat_loss = r_mat_out * r_mat_inner.';
x_mat_loss = x_mat_out * x_mat_inner.';
rx_mat_loss = r_mat_loss + x_mat_loss;
%nw.rx_mat_loss = rx_mat_loss;

%A1 = eye(2*nl);
%B1 = zeros(2*nl,2*ne);

%A2 = blkdiag(downstream_matrix_p,downstream_matrix_p);
%B2 = [big_r_mat.' , big_r_mat.';big_x_mat.' big_x_mat.']/u0;

A3 = -2*[R_sens_mat.' X_sens_mat.'];
B3 = -[z_mat_loss.'+2*rx_mat_loss , z_mat_loss.'+2*rx_mat_loss]/u0;

zz = zeros(1,nl);
dd = downstream_matrix_p(1,:);

%nw.A = [A1;A2([1 nl+1],:);A3];
nw.A = [[dd zz;zz dd];A3];
%nw.B = [B1;B2([1 nl+1],:);B3];
nw.B = [[br.' br.';bx.' bx.'];B3];
%nw.y0 = [zeros(2*nl+2,1);u0*ones(size(B3,1),1)];
nw.y0 = [0;0;u0*ones(size(B3,1),1)];

if nargout > 1
    dbrdt = zeros(1,ne,2*nl);
    dbxdt = zeros(1,ne,2*nl);
    for i=1:nl
        dbrdt(:,:,i) = big_r_mat_coeff(i,1:ne);
        dbxdt(:,:,nl+i) = big_r_mat_coeff(i,1:ne);
    end

    dd = zeros(nn,nv);
    dR_sens_mat_dt = zeros(nn,nv,nl*2);
    dX_sens_mat_dt = zeros(nn,nv,nl*2);
    for i=1:nl
        dd(:) = sens_mat_r_coeff(i,:).';
        dR_sens_mat_dt(:,:,i) = dd;
        dX_sens_mat_dt(:,:,nl+i) = dd;
    end

    dd = zeros(ne,nv);
    dz_mat_loss_dt = zeros(ne,nv,2*nl);
    for i=1:nl
        dd(:) = z_mat_coeff(i,:).'; 
        dz_mat_loss_dt(:,:,i) = (2 * r(i)) * dd;
        dz_mat_loss_dt(:,:,nl+i) = (2 * x(i)) * dd;
    end

    dd = zeros(ne,nl);
    drx_mat_loss = zeros(nv,ne,2*nl);
    for i=1:nl
        dd(:) = rhat_inner_coeff(i,:).';
        ee = zeros(1,nl);
        ee(i) = 1;
        d1 = ee .* downstream_matrix_of_V;
        drx_mat_loss(:,:,i) = d1 * r_mat_inner.' + r_mat_out * dd.';
        drx_mat_loss(:,:,nl+i) = d1 * x_mat_inner.' + x_mat_out * dd.';
    end

    dAdt = zeros(size(nw.A,1),size(nw.A,2),2*nl);
    dBdt = zeros(size(nw.B,1),size(nw.B,2),2*nl);
    for i=1:2*nl
        dAdt(3:end,:,i) = -2*[dR_sens_mat_dt(:,:,i).' dX_sens_mat_dt(:,:,i).'];
        d1  = [dbrdt(:,:,i) dbrdt(:,:,i);dbxdt(:,:,i) dbxdt(:,:,i)];
        dB3 = -[dz_mat_loss_dt(:,:,i).'+2*drx_mat_loss(:,:,i) , dz_mat_loss_dt(:,:,i).'+2*drx_mat_loss(:,:,i)]/u0;
        dBdt(:,:,i) = [d1;dB3];
    end
end