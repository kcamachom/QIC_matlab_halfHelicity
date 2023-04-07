%%% Generates input vmec files for eps_eff scans
%%%%
nfp = 2;

%Z1 = 0.05:0.05:0.8;
%Z1 = 0.2;
Z3 = 0.01:0.01:0.1;
R4 = 0:0.005:0.025; %-0.03:0.005:0.03;
Z4_param = 0:0.1:1.0;


% % %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % %         % %%%%%  1st order zeros kappa at Bmin %%%%%%%%% 
% % %         % %%%%%  2nd order zeros kappa at Bmax %%%%%%%%% 
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
count = 0;
for i=1:length(Z3)%length(B0_1)
    %Z2 = 0:0.005:Z1(i)/10;
    for j=1:length(R4)
    for k = 1:length(Z4_param)
    count = count+1;
    %base_name = 'omni_Scans_Z2_',Z2;
    R0c= [1.0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
    Z0s = [0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
    R0c(4) = R4(j);
    Z0s(3) = Z3(i);
    Z0s(4) = Z4_param(k)*R0c(4);
    R0c(2) = -((R0c(4) + 9*nfp*nfp* R0c(4))/(1 + nfp*nfp));
    R0c(3) = -(1/(1 + 4*nfp*nfp));
    Z0s(2) = -((4*Z0s(3) + 8*nfp*nfp*Z0s(3) + 6*Z0s(4)+ 27*nfp*nfp*Z0s(4))/(2 + nfp^2));
    R0s = zeros(size(R0c));
    Z0c = zeros(size(R0c));

    [varphi,phi,iota,alpha_RMS,curvature, torsion, n_for_alpha, sigma_max, elongation] = mGarrenBoozerOmnigenity_scans_halfHelicity(nfp,R0c,R0s,Z0c,Z0s);
%    [varphi,phi,iota,alpha_RMS,curvature, torsion, n_for_alpha, sigma_max, elongation] = mGarrenBoozerOmnigenity_1layerIt_Frenet_FieldPeriods_alpha(nfp,R0c,R0s,Z0c,Z0s,B0_1(i),B0_2(j));
    Scans(count).R4 = R0c(4);
    Scans(count).R2 = R0c(2);
    Scans(count).R3 = R0c(3);
    Scans(count).Z3 = Z0s(3);
    Scans(count).Z4 = Z0s(4);    
    Scans(count).Z1 = Z0s(1);
    Scans(count).Z2 = Z0s(2);

%     Scans(count).B0_1 = B0_1(i);
%     Scans(count).B0_2 = B0_2(j);
    Scans(count).phi = phi;
    Scans(count).varphi = varphi;
    Scans(count).iota = iota;
    Scans(count).alpha_RMS = alpha_RMS;
    Scans(count).curvature = curvature;
    Scans(count).torsion = torsion;
    Scans(count).n_for_alpha = n_for_alpha;
    Scans(count).sigma_max = sigma_max;
    Scans(count).elongation = elongation;
    end
    end
end
%save('Scans_N5_Z1_0.0_0.05_0.8_Z2_0_0.005_Z1.mat', 'Scans')
