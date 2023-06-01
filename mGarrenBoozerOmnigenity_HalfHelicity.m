%function mGarrenBoozerOmnigenity_1layerIt_Frenet_FieldPeriods_alpha(R0c1, R0c2,R0c3,Z0s1,Z0s2,Z0s3,Z0s4,Z0s5)
%function mGarrenBoozerOmnigenity_1layerIt_Frenet_FieldPeriods_alpha(RZvalues)
function mGarrenBoozerOmnigenity_HalfHelicity()

% See "20190701-01 Omnigenity construction with single layer iteration.docx"

% N_phi = Number of grid points in phi
N_phi_array = 201;
%N_phi_array = [31];
%N_phi_array = [7, 15];

N_theta = 32;%256; 


nfp = 5; % For now only nfp=1 is allowed.

% Shape of the magnetic axis:
% R0c = [1.0,0,-0.2,0,0];
% Z0s = [0,0,0.35,0,0.0];
%Z0s = [0,0.5,0.0,0.05,0.0];
%Z0s = [0,0.0,0.15,0,0.0];
%Z0s = [0,0.0,0.07,0.0,0.0];
R0s = [0,0,0,0,0,0,0,0,0,0];
Z0c = [0,0,0,0,0,0,0,0,0,0];
%R0c = [1.0,0,-1/101,0,0];
%R0c = [5.05,0,-0.05,0,0];
% R0c = [2.04,0,-0.12,0,0];
% %R0c = [1.0,0,-1/17,0,0];
% %R0c = [4.16,0,-0.064,0,0];
% Z0s = [0,0.8,0.01,0,0.0];
% 
% Z0s = [0,0.0,0.05,0.0,0,0,0,0,0,0];

%Z0s = [0,0.4,-0.03,-0.005,0,0,0,0,0,0];
%Z0s = [0,-0.1,-0.1,0.03,0,0,0,0,0,0];
%Z0s = [0,1.0e-01, -5.0e-02, -5.0e-03,0,0,0,0,0,0]; %%% 2FP optimization
%new input

%    R0c = [1,0,0,0,0,0,0,0,0,0];
%    R0c(3) = -1/(1+(4*(nfp^2)));
%     
%            
%    %R0c(2) = -0.2;        
%    %R0c(3) = -1/((4*(nfp^2)) +1);
%    %R0c(4) = -(R0c(2)+ (R0c(2)*(nfp^2)) )/((9*(nfp^2)) +1);    
%    %Z0s(3) = -0.3;
%    %Z0s(3) = 0.1;
% % % % %R0c = 2.04.*R0c;
% % % % %Z0s = [0,0,-0.085,0,0,0,0,0,0,0];
%     Z0s = [0,0.8,0.01,0,0,0,0,0,0,0]./2.04;
% Z0s = [0,0.72,0.01,0,0,0,0,0,0,0];
%Z0s = [0,0.4,-0.1,0,0,0,0,0,0,0];
%    R0c = [1.0,0,-0.05,0,-0.15/65,0,0,0,0,0];
%    Z0s = [0,0,0.05,0,0,0,0,0,0,0];
%  % %Z0s = [0,0.4,0.01,0,0.0];
%     Z0s = [0,0.0,0.07,0,0.0];

%    R0c = [1.0,0,-0.01,0,0.01/401];
%    Z0s = [0,0.1,0,0,0.0];

%%% half-helicity curves

R0c= [1.0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
Z0s = [0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
Z0s(3) = 0.007;%0.0065; % 0.003;
R0c(4) = -0.001;
Z0s(4) = 0.0026;%0.0027;
R0c(2) = -((R0c(4) + 9*nfp*nfp* R0c(4))/(1 + nfp*nfp));
R0c(3) = -(1/(1 + 4*nfp*nfp));
Z0s(2) = -((4*Z0s(3) + 8*nfp*nfp*Z0s(3) + 6*Z0s(4)+ 27*nfp*nfp*Z0s(4))/(2 + nfp^2));

R0s = zeros(size(R0c));
Z0c = zeros(size(R0c));


        
%%%%%%%%%%%% 2 Parameters R %%%%%%%%%%
% % % R0c(1) = 1;
% % % R0c(3) = -1/(1+(4*(nfp^2)));
% % %         Z0s(2) = Axis_Shape(1);
% % %         Z0s(3) = Axis_Shape(2);
% % %         Z0s(4) = Axis_Shape(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B0_1 = 0.12;
B0_2 = -0.002;

B0_as_a_function_of_varphi = @(varphi) 1 + B0_1 * cos(varphi) + B0_2*cos(varphi*2);

% Fourier coefficients of d(varphi):
%ds = [0.7, 0.2, 0.35];
%ds = [0.77, 0.255, 0.4]; %% for nfp =2
%ds = [0.715,0.34,0.49]; %for nfp=4, example
%%%%ds =  [0.41,0.12,0.16];
%%ds =  [0.5,0.1,0.0];
%ds = [1.08,0.26,0.46];
% Fraction of each toroidal period devoted to the left buffer region and to
% the right buffer region. The two buffer regions together have a relative
% width of 2*delta, so delta must be <0.5.
delta = 0.1;

% n_for_alpha is an integer that describes how much alpha increases per toroidal period, divided by (2*pi):
%n_for_alpha = 1; 

iota_initial = 0.7;

alpha_type = 2; 
% approach taken to calculate alpha
% 1 = buffer region approach
% 2 = non-zone approach, requires defining k 

k_for_alpha = 2;
alpha_shift = 5*pi/2;


%finite_r_option='linear';
finite_r_option='nonlinear';

extra_condition_option = 1;
% 1 = force sigma=0 at the phi specified by interpolateTo0
% 2 = force sigma=0 at the phi specified by interpolateToPi

sigma_initial = 0;


% Maximum number of Newton iterations to attempt:
N_iterations = 80;


N_line_search = 6;

sign_G = 1;

I2_over_BBar = 0;
%verify_Jacobian = true;
verify_Jacobian = false;

tolerance = 1e-30;

% r is used only for drawing the flux surface shapes at the end.
aspect_ratio = 12;

r = 1.0 / aspect_ratio;

vmec_filename_base = 'input.omni_1layerFrenet';

figure_offset = 40;
more 

% -------------------------------------------------------
% End of input parameters.
% -------------------------------------------------------

assert(delta <= 0.5)
%assert(n_for_alpha == round(n_for_alpha))

figure(1 + figure_offset)

clf
set(gcf,'Units','in','Position',[0,0.5,17.7639,9.2778])
numRows = 4;
numCols = 8;

iota_vs_phi = zeros(size(N_phi_array));

for which_resolution = 1:numel(N_phi_array)
    N_phi = N_phi_array(which_resolution);
    if mod(N_phi,2)==0 % Ensure N_phi is odd.
        N_phi = N_phi + 1;
    end
    fprintf('------- Beginning solve for N_phi = %d -------\n',N_phi)
    tic
    
    scheme = 20;
    [phi, ~, D, ~] = sfincs_uniformDiffMatrices(N_phi, 0, 2*pi/abs(nfp), scheme);  %%% change Katia    
    %    [phi, ~, D, ~] = sfincs_uniformDiffMatrices(N_phi, 0, 2*pi, scheme);
    assignin('base','D',D)
    dphi = phi(2) - phi(1);
    phi = phi + dphi/3;
    assignin('base','phi_1layer_Frenet',phi)

%     phi_original = phi;
%     %%% Change Katia 
%     
%     phi_Frenet = phi;
%     assignin('base','phiF_1layer_Frenet',phi_Frenet)
     phi_Complete = phi;
     if nfp>1
         for jkl = 2:nfp
             phi_Complete = [phi_Complete;phi+((2*pi*(jkl-1))/nfp)];        
         end
     else
         phi_Complete = phi;
 
     end    
    
%     interpolateTo0 = m20180624_02_FourierInterpolationMatrix(N_phi, -(1/3)*2*pi/(nfp*N_phi));%%% change Katia
%     interpolateToPi = m20180624_02_FourierInterpolationMatrix(N_phi, (pi-(1/3)*2*pi/(N_phi))./nfp);%%% change Katia

    interpolateTo0 = m20180624_02_FourierInterpolationMatrix(N_phi, -(1/3)*2*pi/N_phi);
    interpolateToPi = m20180624_02_FourierInterpolationMatrix(N_phi, pi-(1/3)*2*pi/N_phi);

    phi_extended = linspace(0,2*pi,N_phi*abs(nfp)+1); %%%Katia change
    phi_extended(end) = [];
    
    R0_extended = zeros(size(phi_extended));
    Z0_extended = zeros(size(phi_extended));
    for imn = 1:numel(R0c)
        angle = nfp * (imn-1) * phi_extended;
        sinangle = sin(angle);
        cosangle = cos(angle);
        R0_extended = R0_extended + R0c(imn)*cosangle + R0s(imn)*sinangle;
        Z0_extended = Z0_extended + Z0c(imn)*cosangle + Z0s(imn)*sinangle;
    end
    
    X0_extended = R0_extended .* cos(phi_extended);
    Y0_extended = R0_extended .* sin(phi_extended);
    
    if which_resolution==numel(N_phi_array)
        subplot(numRows,numCols,1)
        plot3(X0_extended, Y0_extended, Z0_extended, '.-')
        daspect([1,1,1])
        zoom(1.7)
        set(gca,'clipping','off')
        axis vis3d
        title('Magnetic axis')
        rotate3d on
    end
    
    
    %phi = phi_Frenet; %change Katia, to make sure Frenet quantities are calculated in the whole angle phi
    
    R0 = zeros(size(phi));
    Z0 = zeros(size(phi));
    R0p = zeros(size(phi));
    Z0p = zeros(size(phi));
    R0pp = zeros(size(phi));
    Z0pp = zeros(size(phi));
    R0ppp = zeros(size(phi));
    Z0ppp = zeros(size(phi));
     
    
    for imn = 1:numel(R0c)
        n = nfp*(imn-1); %%% change Katia
        angle = n * phi;
        sinangle = sin(angle);
        cosangle = cos(angle);
        R0 = R0 + R0c(imn)*cosangle + R0s(imn)*sinangle;
        Z0 = Z0 + Z0c(imn)*cosangle + Z0s(imn)*sinangle;
        R0p = R0p - n*R0c(imn)*sinangle + n*R0s(imn)*cosangle;
        Z0p = Z0p - n*Z0c(imn)*sinangle + n*Z0s(imn)*cosangle;
        R0pp = R0pp - n*n*R0c(imn)*cosangle - n*n*R0s(imn)*sinangle;
        Z0pp = Z0pp - n*n*Z0c(imn)*cosangle - n*n*Z0s(imn)*sinangle;
        R0ppp = R0ppp + n*n*n*R0c(imn)*sinangle - n*n*n*R0s(imn)*cosangle;
        Z0ppp = Z0ppp + n*n*n*Z0c(imn)*sinangle - n*n*n*Z0s(imn)*cosangle;
    end
    
    
    d_l_d_phi = sqrt(R0 .* R0 + R0p .* R0p + Z0p .* Z0p);
    d2_l_d_phi2 = (R0 .* R0p + R0p .* R0pp + Z0p .* Z0pp) ./ d_l_d_phi;
    %B0_over_abs_G = 1 / mean(d_l_d_phi);
    
    d_r_d_phi_cylindrical = [R0p, R0, Z0p];
    d2_r_d_phi2_cylindrical = [R0pp-R0, 2*R0p, Z0pp];
    d3_r_d_phi3_cylindrical = [R0ppp-3*R0p, 3*R0pp-R0, Z0ppp];
    
    tangent_cylindrical = [...
        d_r_d_phi_cylindrical(:,1) ./ d_l_d_phi, ...
        d_r_d_phi_cylindrical(:,2) ./ d_l_d_phi, ...
        d_r_d_phi_cylindrical(:,3) ./ d_l_d_phi];
    
    d_tangent_d_l_cylindrical = [...
        (-d_r_d_phi_cylindrical(:,1) .* d2_l_d_phi2 ./ d_l_d_phi + d2_r_d_phi2_cylindrical(:,1)) ./ (d_l_d_phi .* d_l_d_phi), ...
        (-d_r_d_phi_cylindrical(:,2) .* d2_l_d_phi2 ./ d_l_d_phi + d2_r_d_phi2_cylindrical(:,2)) ./ (d_l_d_phi .* d_l_d_phi), ...
        (-d_r_d_phi_cylindrical(:,3) .* d2_l_d_phi2 ./ d_l_d_phi + d2_r_d_phi2_cylindrical(:,3)) ./ (d_l_d_phi .* d_l_d_phi)];
    
    curvature = sqrt(...
        d_tangent_d_l_cylindrical(:,1) .* d_tangent_d_l_cylindrical(:,1) + ...
        d_tangent_d_l_cylindrical(:,2) .* d_tangent_d_l_cylindrical(:,2) + ...
        d_tangent_d_l_cylindrical(:,3) .* d_tangent_d_l_cylindrical(:,3));
    
    curvature_alt = sqrt( ...
        (d_r_d_phi_cylindrical(:,2) .* d2_r_d_phi2_cylindrical(:,3) - d_r_d_phi_cylindrical(:,3) .* d2_r_d_phi2_cylindrical(:,2)) .^ 2 + ...
        (d_r_d_phi_cylindrical(:,3) .* d2_r_d_phi2_cylindrical(:,1) - d_r_d_phi_cylindrical(:,1) .* d2_r_d_phi2_cylindrical(:,3)) .^ 2 +...
        (d_r_d_phi_cylindrical(:,1) .* d2_r_d_phi2_cylindrical(:,2) - d_r_d_phi_cylindrical(:,2) .* d2_r_d_phi2_cylindrical(:,1)) .^ 2 ) ...
        ./ (d_l_d_phi .^ 3);
    
    kappa_cylindrical = [...
        d_tangent_d_l_cylindrical(:,1) ./ curvature, ...
        d_tangent_d_l_cylindrical(:,2) ./ curvature, ...
        d_tangent_d_l_cylindrical(:,3) ./ curvature];
    
    binormal_cylindrical = [...
        tangent_cylindrical(:,2) .* kappa_cylindrical(:,3) - tangent_cylindrical(:,3) .* kappa_cylindrical(:,2), ...
        tangent_cylindrical(:,3) .* kappa_cylindrical(:,1) - tangent_cylindrical(:,1) .* kappa_cylindrical(:,3), ...
        tangent_cylindrical(:,1) .* kappa_cylindrical(:,2) - tangent_cylindrical(:,2) .* kappa_cylindrical(:,1)];
    
    % This script uses the sign convention for torsion used in wikipedia and mathworld.wolfram.com/Torsion.html but
    % opposite to Garren & Boozer's sign convention.
    torsion_numerator = (0 ...
        + d_r_d_phi_cylindrical(:,1) .* (d2_r_d_phi2_cylindrical(:,2) .* d3_r_d_phi3_cylindrical(:,3) - d2_r_d_phi2_cylindrical(:,3) .* d3_r_d_phi3_cylindrical(:,2)) ...
        + d_r_d_phi_cylindrical(:,2) .* (d2_r_d_phi2_cylindrical(:,3) .* d3_r_d_phi3_cylindrical(:,1) - d2_r_d_phi2_cylindrical(:,1) .* d3_r_d_phi3_cylindrical(:,3)) ...
        + d_r_d_phi_cylindrical(:,3) .* (d2_r_d_phi2_cylindrical(:,1) .* d3_r_d_phi3_cylindrical(:,2) - d2_r_d_phi2_cylindrical(:,2) .* d3_r_d_phi3_cylindrical(:,1)));
    
    torsion_denominator = 0 ...
        + (d_r_d_phi_cylindrical(:,2) .* d2_r_d_phi2_cylindrical(:,3) - d_r_d_phi_cylindrical(:,3) .* d2_r_d_phi2_cylindrical(:,2)) .^ 2 ...
        + (d_r_d_phi_cylindrical(:,3) .* d2_r_d_phi2_cylindrical(:,1) - d_r_d_phi_cylindrical(:,1) .* d2_r_d_phi2_cylindrical(:,3)) .^ 2 ...
        + (d_r_d_phi_cylindrical(:,1) .* d2_r_d_phi2_cylindrical(:,2) - d_r_d_phi_cylindrical(:,2) .* d2_r_d_phi2_cylindrical(:,1)) .^ 2;
    
    torsion = torsion_numerator ./ torsion_denominator;
    
    int_torsion= trapz(phi,torsion);
    %%%%%%% Check for curvature
%     
%     
%     numerador = (1 - 2.4*cos(2*phi) + 2.6*(sin(2*phi)).^2 + 2.33*(cos(2*phi).^2) ...
%     - 1.46*(cos(2*phi).^3) + 0.416*(sin(2*phi).^4) + 0.53*(cos(2*phi).^4) ...
%     + 0.9904*(sin(2*phi).^2).*(cos(2*phi).^2) - 0.768*(sin(2*phi).^2).*cos(2*phi));
% 
%     denominador = (0.16*(sin(2*phi).^2) + 1 + 0.53*(cos(2*phi).^2)- 0.4*cos(2*phi));
% 
% 
%     curvature2 = (numerador.^0.5)./(denominador.^1.5);
%     
%     numerador2 = ( -2.1*cos(2*phi) + 0.84*(cos(2*phi).^2) ...
%     + 1.26*(cos(2*phi).^3)  + 3.36*(sin(2*phi)).^2   ...
%     + 1.344*(sin(2*phi).^2).*(cos(2*phi)) );
% 
%     torsion2 = numerador2./(numerador);
    
    figure(999) 
    plot(phi,curvature)
    hold on 
    plot(phi, curvature_alt, '-.')
   
    plot(phi,torsion,'-.')
    hold on
    %plot(phi, torsion2)
    legend('curvature','curvature alt','torsion')
    
    
% %     figure(999)
% %     plot(phi,d_l_d_phi)
% %     legend('dl/dphi')
% %     hold on 
    
    
    % Now convert the 3 vectors from cylindrical to Cartesian coordinates, for plotting:
    % We could have computed the binormal directly in Cartesian coordinates,
    % but I wanted to do it first in cylindrical coordinates since we may want
    % to project results from the Garren-Boozer calculation into the (R,Z)
    % plane.
    
    cosphi = cos(phi);
    sinphi = sin(phi);
    
    tangent_Cartesian = [...
        tangent_cylindrical(:,1) .* cosphi - tangent_cylindrical(:,2) .* sinphi, ...
        tangent_cylindrical(:,1) .* sinphi + tangent_cylindrical(:,2) .* cosphi, ...
        tangent_cylindrical(:,3)];
    
    kappa_Cartesian = [...
        kappa_cylindrical(:,1) .* cosphi - kappa_cylindrical(:,2) .* sinphi, ...
        kappa_cylindrical(:,1) .* sinphi + kappa_cylindrical(:,2) .* cosphi, ...
        kappa_cylindrical(:,3)];
    
    binormal_Cartesian = [...
        binormal_cylindrical(:,1) .* cosphi - binormal_cylindrical(:,2) .* sinphi, ...
        binormal_cylindrical(:,1) .* sinphi + binormal_cylindrical(:,2) .* cosphi, ...
        binormal_cylindrical(:,3)];
    
    phi_tripled  = [ phi-(2*pi/nfp);  phi;  phi+(2*pi/nfp)]; %%% Katia change
    indices_for_sign_flip = (phi > pi/nfp);
     R0_tripled = repmat(R0,3,1);
     Z0_tripled = repmat(Z0,3,1);
    kappa_cylindrical_signed = kappa_cylindrical;
    binormal_cylindrical_signed = binormal_cylindrical;
    kappa_cylindrical_signed(indices_for_sign_flip,:) = -kappa_cylindrical_signed(indices_for_sign_flip,:);
    binormal_cylindrical_signed(indices_for_sign_flip,:) = -binormal_cylindrical_signed(indices_for_sign_flip,:);
    kappa_cylindrical_signed_tripled = repmat(-kappa_cylindrical_signed,3,1);
    binormal_cylindrical_signed_tripled = repmat(-binormal_cylindrical_signed,3,1);
    
    indices_for_sign_flip_tripled = [N_phi+1:2*N_phi];
    binormal_cylindrical_signed_tripled(indices_for_sign_flip_tripled,:) = -binormal_cylindrical_signed_tripled(indices_for_sign_flip_tripled,:);
    kappa_cylindrical_signed_tripled(indices_for_sign_flip_tripled,:) = -kappa_cylindrical_signed_tripled(indices_for_sign_flip_tripled,:);
    
    curvature_signed = curvature;
    curvature_signed(indices_for_sign_flip) = -curvature(indices_for_sign_flip);
    curvature = curvature_signed;
    
    cosphi = cos(phi_tripled);
    sinphi = sin(phi_tripled);
    
    kappa_Cartesian_signed_tripled = [...
        kappa_cylindrical_signed_tripled(:,1) .* cosphi - kappa_cylindrical_signed_tripled(:,2) .* sinphi, ...
        kappa_cylindrical_signed_tripled(:,1) .* sinphi + kappa_cylindrical_signed_tripled(:,2) .* cosphi, ...
        kappa_cylindrical_signed_tripled(:,3)];
    
    binormal_Cartesian_signed_tripled = [...
        binormal_cylindrical_signed_tripled(:,1) .* cosphi - binormal_cylindrical_signed_tripled(:,2) .* sinphi, ...
        binormal_cylindrical_signed_tripled(:,1) .* sinphi + binormal_cylindrical_signed_tripled(:,2) .* cosphi, ...
        binormal_cylindrical_signed_tripled(:,3)];
    
    
%     figure(7)
%     plot(kappa_Cartesian_signed_tripled)
%     hold on 
%     figure(8)
%     plot(binormal_Cartesian_signed_tripled)
%     tangent_Cartesian(1,:)
%     kappa_Cartesian(1,:)
    
%%%%%%%%Helicity Calculation 2 %%%%%%%


d_kappa_d_phi = (-curvature_signed.*tangent_cylindrical + torsion.* binormal_cylindrical_signed);
integrand = kappa_cylindrical_signed(:,3).*d_kappa_d_phi(:,1) - kappa_cylindrical_signed(:,1).*d_kappa_d_phi(:,3);
helicity = cumtrapz(phi,integrand);
helicity = helicity(end)/(2*pi);


%%%%%% Helicity calculation %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quadrant = zeros(size(phi));
    for i=1:size(phi)
        if kappa_cylindrical_signed(i,1)>=0
            if kappa_cylindrical_signed(i,3)>=0
                quadrant(i) = 1;
            else 
                quadrant(i) = 4;
            end
        else 
            if kappa_cylindrical_signed(i,3)>=0
                quadrant(i) = 2;
            else
                quadrant(i) = 4;
            end
        end
    end
    
    crosses = 0;
    for j = 1:size(phi)-1
        if quadrant(j) == 4 && quadrant(j+1) == 1
            crosses = crosses +1;
        elseif quadrant(j) == 1 && quadrant(j+1) == 4
            crosses = crosses -1;
        else
            crosses = crosses + quadrant(j+1) -quadrant(j);
        end
    end
    
    if quadrant(length(phi)) == 4 && quadrant(1) == 1
            crosses = crosses +1;
        elseif quadrant(length(phi)) == 1 && quadrant(1) == 4
            crosses = crosses -1;
        else
            crosses = crosses + quadrant(1) -quadrant(length(phi));
    end
    
    
 n_for_alpha = -nfp*crosses/4;
 n_for_alpha = 0.5;

 figure(357)
    plot(kappa_cylindrical_signed(:,1),kappa_cylindrical_signed(:,3))

 figure(358)
    plot(kappa_cylindrical(:,1),kappa_cylindrical(:,3))
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Curvature fit 

[A_n,B_n,~] = Fseries(phi,curvature_signed,2,true,'sin');

[A_tau,B_tau,~] = Fseries(phi,torsion,10,true,'cos');


    % -------------------------------------------------------------
    % Compute the Boozer toroidal angle varphi on the phi grid.

    fprintf('Computing varphi(phi)...\n')
    
    % The operator D is singular, so we add the equation varphi(phi=0)=0 to the first row.
    %figure(1); clf   
    
    matrix = D;
    matrix(1,:) = matrix(1,:) + interpolateTo0(1,:);
    nu = zeros(N_phi,1); % initial guess
    for j = 1:20
        varphi = phi + nu; 
        B0 = B0_as_a_function_of_varphi(nfp*varphi);
        abs_G0 = nfp*sum(B0 .* d_l_d_phi) * dphi /(2*pi);
        rhs = -1 + d_l_d_phi .* B0 / abs_G0;
        last_nu = nu;
        nu = matrix \ rhs;
        %plot(phi,nu,'.-','displayname',['iteration ',num2str(j)]); hold on
        norm_change = sqrt(sum((nu-last_nu).^2)/N_phi);
        fprintf('  Iteration %3d: |change to nu| = %g\n',j,norm_change)
        if norm_change < 1e-15
            break
        end
    end
    legend('location','eastoutside')
    
    varphi = phi + nu;  
    B0 = B0_as_a_function_of_varphi(nfp*varphi);
    abs_G0 = nfp*sum(B0 .* d_l_d_phi) * dphi / (2*pi); 
    G0 = sign_G * abs_G0;
    B0_over_abs_G = B0 ./ abs_G0;
    assignin('base','B0_1layer_Frenet',B0)
    assignin('base','varphi_1layer_Frenet',varphi)
    fprintf('Done computing varphi(phi).\n')
    
    % -------------------------------------------------------------    
    % Now that we know varphi(phi), compute d on the phi grid:
     d = zeros(size(phi));
     
     d_k = 0.28;
     dk_cos = -0.065;
     dk_sin = 0.04;
    
    for j= 1:size(phi)
        if phi(j) < pi/nfp
            d(j) = curvature(j).*sqrt(d_k./B0(j))+ curvature(j).*dk_cos.*cos(nfp*varphi(j)) + curvature(j).*dk_sin.*sin(nfp*varphi(j))   ;
        elseif phi(j) > pi/nfp
           d(j) = curvature(j).*sqrt(d_k./B0(j)) + curvature(j).*dk_cos.*cos(nfp*varphi(j)) - curvature(j).*dk_sin.*sin(nfp*varphi(j)) ; 
        end
    end
   
%    for j= 1:size(phi)
%         if phi(j) < pi/nfp
%             d(j) = curvature(j).*0.73;%sqrt(0.3./B0(j));
%         elseif phi(j) > pi/nfp
%            d(j) = -curvature(j).*0.73;%sqrt(0.3./B0(j));%.*0.73; 
%         end
%     end
%     

% d = curvature_signed.*0.73; 

%     for j = 1:numel(ds)
%         d = d + ds(j) * sin(nfp*j*varphi);
%     end
    
    %%%%%%%plot stuff%%
    
    figure(1235)
    %plot(phi,curvature,'--')
    hold on 
    plot(phi,binormal_Cartesian,'--')
    %binormal_Cartesian(1,:)
    %plot(phi,kappa_cylindrical,'-.')


    % -------------------------------------------------------------
    % Compute the part of alpha(varphi) that is independent of iota and the
    % part proportional to iota:
    % ----------------------------------------------


      
    %%% different alpha expressions 
    
    if alpha_type == 2
        
        % -------------------------------------------------------------
        % Compute alpha with the non-zone approach
        % -------------------------------------------------------------
        
        phiNpi_factor = ((varphi*nfp/pi)-1).^((2*k_for_alpha)+1);
        alpha_iota = varphi - (pi/nfp) - (pi/nfp)*phiNpi_factor;
        alpha_notIota = (pi/2) + pi*n_for_alpha*(2+phiNpi_factor);
       
        d_alpha_iota_d_varphi = 1-(2*k_for_alpha+1)*((((varphi*(nfp/pi))-1)).^(2*k_for_alpha));
        d_alpha_notIota_d_varphi = (2*k_for_alpha+1)*n_for_alpha*nfp*(((varphi*(nfp/pi))-1).^(2*k_for_alpha));

        alpha_omnigenous_iota =  varphi - (pi/nfp);
        alpha_omnigenous_notIota = pi*(2*n_for_alpha+0.5);
        
    else
        % -------------------------------------------------------------
        % Compute alpha with the buffer regions method
        % -------------------------------------------------------------
        % Set up some quantities we will need.
        delta2pi = 2*pi*delta/nfp;
        indices_left = find(varphi < delta2pi);
        indices_right = find(varphi > ((2*pi/nfp)-delta2pi));

        % ----------------------------------------------
        % Handle region II:

        %alpha = iota * (varphi - pi):
        alpha_iota = varphi - (pi/nfp);
        alpha_notIota = zeros(size(varphi));

        % ----------------------------------------------
        % Handle the contributions proportional to iota.

        % Find alpha at the boundaries between regions I and II, and between
        % regions II and III:
        %alpha_I_II = interp1(phi_grid, alpha, delta);
        alpha_I_II = (2*pi*delta - pi)/nfp; 

        % Set alpha in region I:
        %alpha_0 = 0.5 * (alpha_I_II + alpha_II_III) - pi * n_for_alpha;
        alpha_0 = 0;

        %{
    alpha_1 = -(2*alpha_0 - 2*alpha_I_II + iota*delta) / delta;
    alpha_3 = 2*(alpha_0 - alpha_I_II + iota*delta) / (delta * delta * delta) ;
    alpha_4 = -(alpha_0 - alpha_I_II + iota*delta) / (delta * delta * delta * delta) ;
        %}

        alpha_1 = -(2*alpha_0 - 2*alpha_I_II + 1*delta2pi) / delta2pi;
        alpha_3 = 2*(alpha_0 - alpha_I_II + 1*delta2pi) / (delta2pi * delta2pi * delta2pi) ;
        alpha_4 = -(alpha_0 - alpha_I_II + 1*delta2pi) / (delta2pi * delta2pi * delta2pi * delta2pi) ;
        alpha_iota(indices_left) = alpha_0 + alpha_1 .* varphi(indices_left) + alpha_3 .* (varphi(indices_left) .^3) + alpha_4 .* (varphi(indices_left) .^4);
        alpha_iota(indices_right) = -(alpha_0 + alpha_1 .* ((2*pi/nfp)-varphi(indices_right)) + alpha_3 .* (((2*pi/nfp)-varphi(indices_right)) .^3) + alpha_4 .* (((2*pi/nfp)-varphi(indices_right)) .^4));


        % ----------------------------------------------
        % Handle the contributions independent of iota.

        % Find alpha at the boundaries between regions I and II, and between
        % regions II and III:
        %alpha_I_II = interp1(phi_grid, alpha, delta);
        alpha_I_II = 0;

        % Set alpha in region I:
        %alpha_0 = 0.5 * (alpha_I_II + alpha_II_III) - pi * n_for_alpha;
        alpha_0 =  - pi * n_for_alpha/nfp;

        %{
    alpha_1 = -(2*alpha_0 - 2*alpha_I_II + iota*delta) / delta;
    alpha_3 = 2*(alpha_0 - alpha_I_II + iota*delta) / (delta * delta * delta) ;
    alpha_4 = -(alpha_0 - alpha_I_II + iota*delta) / (delta * delta * delta * delta) ;
        %}
        alpha_1 = -(2*alpha_0 - 2*alpha_I_II + 0*delta2pi) / delta2pi;
        alpha_3 = 2*(alpha_0 - alpha_I_II + 0*delta2pi) / (delta2pi * delta2pi * delta2pi) ;
        alpha_4 = -(alpha_0 - alpha_I_II + 0*delta2pi) / (delta2pi * delta2pi * delta2pi * delta2pi) ;
        alpha_notIota(indices_left) = alpha_0 + alpha_1 .* varphi(indices_left) + alpha_3 .* ((varphi(indices_left)) .^3) + alpha_4 .* ((varphi(indices_left)) .^4);
        alpha_notIota(indices_right) = -(alpha_0 + alpha_1 .* ( (2*pi/nfp)-varphi(indices_right) ) + alpha_3 .* ( ( (2*pi/nfp)-varphi(indices_right) ).^3) + alpha_4 .* (((2*pi/nfp)-varphi(indices_right)).^4));
        alpha_notIota = alpha_notIota + alpha_shift;

    end    
        
    % -------------------------------------------------------------    
    %d_d_varphi = diag(1 ./ (B0_over_abs_G .* d_l_d_phi)) * D;
    
    d_d_varphi = diag(abs_G0 ./ (B0 .* d_l_d_phi)) * D;

   
%    d_alpha_iota_d_varphi = d_d_varphi * alpha_iota;

%     d_alpha_notIota_d_varphi = n_for_alpha + d_d_varphi * (alpha_notIota - varphi*n_for_alpha); % We have to treat the secular part separately here since D assumes periodicity!
    
%      d_alpha_notIota_d_varphi =  d_d_varphi * (alpha_notIota - varphi); % We have to treat the secular part separately here since D assumes periodicity!
    

    
    figure(4567)
    plot(varphi*nfp,d_alpha_iota_d_varphi,'-.')
    hold on 
    plot(varphi*nfp,d_alpha_notIota_d_varphi,'-.')
    
    figure(1 + figure_offset)
    plotNum = 2;

    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, B0, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('B_0')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, d, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('d')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, d_l_d_phi, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('d l / d \phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, curvature, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    %plot(phi, curvature2, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    xlabel('\phi')
    title('\kappa')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, (d ./ curvature) .^ 2 , '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('d^2 / \kappa^2')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, torsion, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('\tau')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, varphi, '.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('\phi')
    title('$\varphi$','interpreter','latex')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, alpha_iota, '.-')
    title('\alpha iota at iota initial')
    xlabel('\phi')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, alpha_notIota, '.-')
    title('\alpha notIota at iota initial')
    xlabel('\phi')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi, alpha_notIota - varphi * n_for_alpha , '.-')
    title('\alpha notIota  - varphi * n_for_alpha at iota initial')
    xlabel('\phi')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    alpha = alpha_iota * iota_initial + alpha_notIota;
    plot(phi, alpha, '.-')
    title('\alpha total at iota initial')
    xlabel('\phi')
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    B1c_over_B0 = d .* cos(alpha);
    B1s_over_B0 = d .* sin(alpha);
    plot(phi, B1c_over_B0, '.-', 'DisplayName','B1c')
    hold on
    plot(phi, B1s_over_B0, '.-', 'DisplayName','B1s')
    legend show
    hold off
    title('B_1 at iota initial')
    
    % Indices for the rows of the state vector and columns of the matrix:
    indices_sigma  = (1:N_phi) + 0*N_phi;
    index_iota = 1*N_phi + 1;
    
    % Indices for the equations, i.e. rows of the matrix:
    indices_ODE    = (1:N_phi) + 0*N_phi;
    index_theta_constraint = 1*N_phi + 1;
    
    matrix_size = 1*N_phi + 1;
    
    state = zeros(matrix_size,1);
    
    % Initial guess:
    state(indices_sigma) = zeros(N_phi,1);
    state(index_iota) = iota_initial;
    
    % Initialize some quantities so they have global scope:
    sigma = 0;
    iota = iota_initial;
    B1c_over_B0 = 0;
    B1s_over_B0 = 0;
    alpha = 0;
    
    if verify_Jacobian
        state = linspace(0.4,1.7,matrix_size)';
        state0 = state;
        extract_quantities_from_state_vector();
        matrix = compute_matrix();
        finite_difference_Jacobian = zeros(matrix_size);
        
        for j = 1:matrix_size
            delta_finite_difference = 1e-6;
            
            state = state0;
            state(j) = state(j) + delta_finite_difference;
            extract_quantities_from_state_vector();
            residual_plus = compute_residual();
            
            state = state0;
            state(j) = state(j) - delta_finite_difference;
            extract_quantities_from_state_vector();
            residual_minus = compute_residual();
            
            finite_difference_Jacobian(:,j) = (residual_plus - residual_minus) / (2*delta_finite_difference);
        end
        
        th = 3e-9;
        
        figure(10 + figure_offset)
        clf
        
        subplot(1,3,1)
        spy(abs(matrix)>th);
        title('Analytic Jacobian')
        
        subplot(1,3,2)
        spy(abs(finite_difference_Jacobian)>th);
        title('Finite difference Jacobian')
        
        subplot(1,3,3)
        spy(abs(matrix - finite_difference_Jacobian)>th);
        title('Difference')
        
        assignin('base','matrix',matrix)
        assignin('base','finite_difference_Jacobian',finite_difference_Jacobian)
        
        fprintf('max(max(abs(matrix - finite_difference_Jacobian))) = %g\n', max(max(abs(matrix - finite_difference_Jacobian))))
        return
    end
    
    residual_norms = zeros(N_iterations+1,1);
    iotas = zeros(N_iterations+1,1);
    
    
    extract_quantities_from_state_vector();
    residual = compute_residual();
    residual_norm = sqrt(sum(residual.*residual));
    residual_norms(1) = residual_norm;
    iotas(1) = iota;
    fprintf('Initial residual norm: %g  iota: %.15g\n',residual_norms(1), iota)
    
    for iteration = 2:(N_iterations+1)
        last_residual_norm = residual_norm;
        if residual_norm < tolerance
            break
        end
        
        matrix = compute_matrix();
        assignin('base','matrix',matrix)
        
        state0 = state;
        fprintf('Newton iteration %d.\n',iteration-1)
        fprintf('Condition number of Jacobian: %g\n',cond(matrix))
        step_direction = - matrix \ residual;
        
        step_scale = 1.0;
        for j_line_search = 1:N_line_search
            state = state0 + step_scale * step_direction;
            extract_quantities_from_state_vector();
            residual = compute_residual();
            residual_norm = sqrt(sum(residual.*residual));
            fprintf('  Line search step %d. Residual norm: %g  iota: %.15g\n',j_line_search,residual_norm, iota)
            if residual_norm <= last_residual_norm
                break
            end
            step_scale = step_scale / 2;
        end
        
        residual_norms(iteration) = residual_norm;
        iotas(iteration) = iota;
    end
    
    fprintf('Seconds for solve: %g\n',toc)
    iota_vs_phi(which_resolution) = iota;
    


    X1s = B1s_over_B0 ./ curvature;
    X1c = B1c_over_B0 ./ curvature;
    
    Y1s = sign_G * curvature .* ( B1c_over_B0 + B1s_over_B0 .* sigma) ./ ((B1c_over_B0.*B1c_over_B0 + B1s_over_B0.*B1s_over_B0).*B0);
    Y1c = sign_G * curvature .* (-B1s_over_B0 + B1c_over_B0 .* sigma) ./ ((B1c_over_B0.*B1c_over_B0 + B1s_over_B0.*B1s_over_B0).*B0);
   
   
    
    figure(41)
    
    
    R1c = (-binormal_cylindrical_signed(:,3) .* X1c + kappa_cylindrical_signed(:,3) .* Y1c) .* (d_l_d_phi) ./ R0;
    R1s = (-binormal_cylindrical_signed(:,3) .* X1s + kappa_cylindrical_signed(:,3) .* Y1s) .* (d_l_d_phi) ./ R0;
    Z1c = ( binormal_cylindrical_signed(:,1) .* X1c - kappa_cylindrical_signed(:,1) .* Y1c) .* (d_l_d_phi) ./ R0;
    Z1s = ( binormal_cylindrical_signed(:,1) .* X1s - kappa_cylindrical_signed(:,1) .* Y1s) .* (d_l_d_phi) ./ R0;
    
    X1c_signed = X1c;
    X1s_signed = X1s;
    Y1c_signed = Y1c;
    Y1s_signed = Y1s;
    X1c_signed(indices_for_sign_flip) = -X1c_signed(indices_for_sign_flip);
    X1s_signed(indices_for_sign_flip) = -X1s_signed(indices_for_sign_flip);
    Y1c_signed(indices_for_sign_flip) = -Y1c_signed(indices_for_sign_flip);
    Y1s_signed(indices_for_sign_flip) = -Y1s_signed(indices_for_sign_flip);

        angle_untwist = -n_for_alpha * nfp * varphi;
        sinangle = sin(angle_untwist);
        cosangle = cos(angle_untwist);
        X1s_untwisted = X1s;%(X1s .*   cosangle  + X1c .* sinangle);%* self.sign_curvature_change
        X1c_untwisted = X1c;%(X1s .* (-sinangle) + X1c .* cosangle);% * self.sign_curvature_change
        Y1s_untwisted = Y1s;%(Y1s .*   cosangle  + Y1c .* sinangle);% * self.sign_curvature_change
        Y1c_untwisted = Y1c;%(Y1s .* (-sinangle) + Y1c .* cosangle);% * self.sign_curvature_change
% 
%     figure(1001)
%     plot(phi,curvature)
%     figure(1002)
%     plot(phi,binormal_cylindrical_signed)
%     figure(1003)
%     plot(phi,kappa_cylindrical_signed)
%     figure(1004)
%     plot(phi,binormal_cylindrical)
%     figure(1005)
%     plot(phi,kappa_cylindrical)
%     figure(41)
    
    assignin('base','R1c_1layer_Frenet',R1c)
    assignin('base','R1s_1layer_Frenet',R1s)
    assignin('base','Z1c_1layer_Frenet',Z1c)
    assignin('base','Z1s_1layer_Frenet',Z1s)
    
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    semilogy(residual_norms(1:iteration),'o-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    xlabel('Iteration')
    title('Residual norm')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,X1c,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('X1c')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,X1s,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('X1s')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Y1c,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Y1c')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Y1s,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Y1s')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end

    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,X1c_signed,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('X1c signed')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,X1s_signed,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('X1s signed')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Y1c_signed,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Y1c signed')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Y1s_signed,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Y1s signed')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,sigma,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('sigma')
    xlabel('\phi')

    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,R1c,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('R1c')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,R1s,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('R1s')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Z1c,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Z1c')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,Z1s,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('Z1s')
    xlabel('\phi')
    if which_resolution == numel(N_phi_array)
        legend show
        set(legend,'Interpreter','none')
    end
    
    % See my note 20180329-03 for derivation of the formula below for
    % elongation:
    p = R1s.*R1s + R1c.*R1c + Z1s.*Z1s + Z1c.*Z1c;
    q = R1s.*Z1c - R1c.*Z1s;
%     p = X1s.*X1s + X1c.*X1c + Y1s.*Y1s + Y1c.*Y1c;
%     q = X1s.*Y1c - X1c.*Y1s;
    elongation = 2*abs(q)./(p - sqrt(p.*p-4*q.*q));

    subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
    plot(phi,elongation,'.-', 'DisplayName',sprintf('N_phi=%d',N_phi))
    hold on
    title('elongation in RZ plane')
    xlabel('\phi')

end

subplot(numRows,numCols,plotNum); plotNum = plotNum + 1;
plot(N_phi_array,iota_vs_phi,'o-')
xlabel('N_phi','Interpreter','none')
title('Iota')
xlim([0,Inf])

    function extract_quantities_from_state_vector()
        % Extract individual quantities from the state vector:
        sigma = state(indices_sigma);
        iota = state(index_iota);
        alpha = alpha_iota * iota + alpha_notIota; %%%alpha change
        B1c_over_B0 = d .* cos(alpha);
        B1s_over_B0 = d .* sin(alpha);
    end

   function residual = compute_residual()
        residual = zeros(matrix_size,1);
        
        residual(indices_ODE) = d_d_varphi * sigma ...
            + (iota - iota * d_alpha_iota_d_varphi - d_alpha_notIota_d_varphi) .* ( (d .^ 4) .* B0 .* B0 ./ (curvature .^ 4) + 1 + sigma .* sigma) ...
            - 2 * (I2_over_BBar - torsion) .* G0 .* d .* d ./ (curvature .* curvature);
        
        switch extra_condition_option
            case 1
                residual(index_theta_constraint) = interpolateTo0 * sigma;
            case 2
                residual(index_theta_constraint) = interpolateToPi * sigma;
            otherwise
                error('Invalid extra_condition_option')
        end
    end

    function matrix = compute_matrix()
        matrix = zeros(matrix_size);
        
        %matrix(indices_ODE,index_iota) =  B1Squared_over_kappaSquared .* B1Squared_over_kappaSquared + 1 + sigma .* sigma;
        matrix(indices_ODE,index_iota) =  (1 - d_alpha_iota_d_varphi) .* ( (d .^ 4) .* B0 .* B0 ./ (curvature .^ 4) + 1 + sigma .* sigma);
        
        %matrix(indices_ODE,indices_sigma) = d_d_varphi + diag((iota - iota_shift) * 2 * sigma);
        matrix(indices_ODE,indices_sigma) = d_d_varphi + diag((iota - iota * d_alpha_iota_d_varphi - d_alpha_notIota_d_varphi) .* 2 .* sigma);
        
        
        switch extra_condition_option
            case 1
                matrix(index_theta_constraint, indices_sigma) = interpolateTo0;
            case 2
                matrix(index_theta_constraint, indices_sigma) = interpolateToPi;
            otherwise
                error('Invalid extra_condition_option')
        end
    end


% **********************************************
%return
% **********************************************

fprintf('maxval(curvature): %.15g\n',max(curvature))

d_phi = phi(2) - phi(1);
axis_length = sum(d_l_d_phi) * d_phi * nfp;
mean_R = sum(R0.* d_l_d_phi) * d_phi * nfp / axis_length;
mean_Z = sum(Z0.* d_l_d_phi) * d_phi * nfp / axis_length;
standard_deviation_of_R = sqrt(sum((R0 - mean_R) .^ 2 .* d_l_d_phi) * d_phi * nfp / axis_length);
standard_deviation_of_Z = sqrt(sum((Z0 - mean_Z) .^ 2 .* d_l_d_phi) * d_phi * nfp / axis_length);
fprintf('mean_R: %.15g\n',mean_R)
fprintf('mean_Z: %.15g\n',mean_Z)
fprintf('standard_deviation_of_R: %.15g\n',standard_deviation_of_R)
fprintf('standard_deviation_of_Z: %.15g\n',standard_deviation_of_Z)

% % % figure(700)
% % % hold on
% % % theta = linspace(0,2*pi,100);
% % % extract_quantities_from_state_vector();
% % % %colors = jet(N_phi);
% % % for j = 1:100:N_phi
% % %     plot(R0(j),Z0(j),'x','Color','m')
% % %     hold on
% % %     plot(R0(j) + r * (R1c(j)*cos(theta) + R1s(j)*sin(theta)), Z0(j) + r * (Z1c(j)*cos(theta) + Z1s(j)*sin(theta)),'.-','Color','m')
% % % end

figure(2345)
plot(phi,elongation)
hold on 

figure(2346)
plot(phi, (d ./ curvature) .^ 2 , '.-', 'DisplayName',sprintf('dk=%d',d_k))
    hold on     
 plot(phi, (d ./ curvature) , '.-', 'DisplayName',sprintf('dk=%d',d_k))
    title('d^2 / \kappa^2')

figure(7 + figure_offset)
hold on
theta = linspace(0,2*pi,100);
extract_quantities_from_state_vector();
colors = jet(N_phi);
for j = 1:N_phi
    plot(R0(j),Z0(j),'x','Color',colors(j,:))
    hold on
    plot(R0(j) + r * (R1c(j)*cos(theta) + R1s(j)*sin(theta)), Z0(j) + r * (Z1c(j)*cos(theta) + Z1s(j)*sin(theta)),'.-','Color',colors(j,:))
end
xlabel('R')
ylabel('Z')
title('Flux surface shapes for the solution')
axis equal


    function residual = phi_residual(phi0)     
        % phi0 is the value of phi corresponding to zeta0.
        % I.e., at given zeta0, phi0 is the cylindrical azimuthal angle on
        % the axis.
        %phi0 = interp1(zeta_tripled, phi_tripled, zeta0, interp1_method, 'extrap');
        
        % Evaluate Cartesian components of r0(zeta):
        R0_at_phi0 = interp1(phi_tripled, R0_tripled, phi0, interp1_method, 'extrap');
        Cartesian_x_at_phi0 = R0_at_phi0 * cos(phi0);
        Cartesian_y_at_phi0 = R0_at_phi0 * sin(phi0);
        
        % Evaluate Cartesian components of the normal and binormal vectors
        % at the specified zeta0:
        normal_x = interp1(phi_tripled, kappa_Cartesian_signed_tripled(:,1), phi0, interp1_method, 'extrap');
        normal_y = interp1(phi_tripled, kappa_Cartesian_signed_tripled(:,2), phi0, interp1_method, 'extrap');
        binormal_x = interp1(phi_tripled, binormal_Cartesian_signed_tripled(:,1), phi0, interp1_method, 'extrap');
        binormal_y = interp1(phi_tripled, binormal_Cartesian_signed_tripled(:,2), phi0, interp1_method, 'extrap');
        
        %normal_R = interp1(phi, kappa_cylindrical(:,1), phi0, interp1_method, 'extrap');
        %normal_phi = interp1(phi, kappa_cylindrical(:,2), phi0, interp1_method, 'extrap');
        %normal_Z = interp1(phi, kappa_cylindrical(:,3), phi0, interp1_method, 'extrap');
        %binormal_R = interp1(phi, binormal_cylindrical(:,1), phi0, interp1_method, 'extrap');
        %binormal_phi = interp1(phi, binormal_cylindrical(:,2), phi0, interp1_method, 'extrap');
        %binormal_Z = interp1(phi, binormal_cylindrical(:,3), phi0, interp1_method, 'extrap');
        
        % We solved for X1 and Y1 on a grid in phi0, since the d/dzeta
        % operator was transformed to a d/dphi0 operator in the ODE.
        X1_at_phi0 = interp1(phi_tripled, X1, phi0, interp1_method, 'extrap');
        Y1_at_phi0 = interp1(phi_tripled, Y1, phi0, interp1_method, 'extrap');
        total_x = Cartesian_x_at_phi0 + r * (X1_at_phi0 * normal_x + Y1_at_phi0 * binormal_x);
        total_y = Cartesian_y_at_phi0 + r * (X1_at_phi0 * normal_y + Y1_at_phi0 * binormal_y);
        
        %residual = phi0 + r * (X1_at_zeta0 * normal_phi + Y1_at_zeta0 * binormal_phi) - phi_1D(j_phi);
        residual = atan2(total_y, total_x) - phi_1D(j_phi);
        % We expect the residual to be less than pi in absolute value, so if it is not, the reason must be the branch cut:
        if (residual >  pi) 
            residual = residual - 2*pi;
        end
        if (residual < -pi) 
            residual = residual + 2*pi;
        end

    end


    function Rxyz = Rxyz(phi0)     
        % phi0 is the value of phi corresponding to zeta0.
        % I.e., at given zeta0, phi0 is the cylindrical azimuthal angle on
        % the axis.
        %phi0 = interp1(zeta_tripled, phi_tripled, zeta0, interp1_method, 'extrap');
        
        % Evaluate Cartesian components of r0(zeta):
        R0_at_phi0 = interp1(phi_tripled, R0_tripled, phi0, interp1_method, 'extrap');
        Z0_at_phi0 = interp1(phi_tripled, Z0_tripled, phi0, interp1_method, 'extrap');
        Cartesian_x_at_phi0 = R0_at_phi0 * cos(phi0);
        Cartesian_y_at_phi0 = R0_at_phi0 * sin(phi0);
        
        % Evaluate Cartesian components of the normal and binormal vectors
        % at the specified zeta0:
        normal_x = interp1(phi_tripled, kappa_Cartesian_signed_tripled(:,1), phi0, interp1_method, 'extrap');
        normal_y = interp1(phi_tripled, kappa_Cartesian_signed_tripled(:,2), phi0, interp1_method, 'extrap');
        normal_z = interp1(phi_tripled, kappa_Cartesian_signed_tripled(:,3), phi0, interp1_method, 'extrap');
        binormal_x = interp1(phi_tripled, binormal_Cartesian_signed_tripled(:,1), phi0, interp1_method, 'extrap');
        binormal_y = interp1(phi_tripled, binormal_Cartesian_signed_tripled(:,2), phi0, interp1_method, 'extrap');
        binormal_z = interp1(phi_tripled, binormal_Cartesian_signed_tripled(:,3), phi0, interp1_method, 'extrap');
        
        % We solved for X1 and Y1 on a grid in phi0, since the d/dzeta
        % operator was transformed to a d/dphi0 operator in the ODE.
        X1_at_phi0 = interp1(phi_tripled, X1, phi0, interp1_method, 'extrap');
        Y1_at_phi0 = interp1(phi_tripled, Y1, phi0, interp1_method, 'extrap');
        total_x = Cartesian_x_at_phi0 + r * (X1_at_phi0 * normal_x + Y1_at_phi0 * binormal_x);
        total_y = Cartesian_y_at_phi0 + r * (X1_at_phi0 * normal_y + Y1_at_phi0 * binormal_y);
        total_z =          Z0_at_phi0 + r * (X1_at_phi0 * normal_z + Y1_at_phi0 * binormal_z);
        
        %residual = phi0 + r * (X1_at_zeta0 * normal_phi + Y1_at_zeta0 * binormal_phi) - phi_1D(j_phi);
        Rxyz = [sqrt(total_x*total_x + total_y*total_y), total_x, total_y, total_z];
    end


% -------------------------------------------------------------

%ntor = 15;
ntor = min([100,floor(N_phi / 2)]);

if strcmp(finite_r_option,'linear')
    rbc1 = zeros(2*ntor+1,1);
    rbs1 = zeros(2*ntor+1,1);
    zbc1 = zeros(2*ntor+1,1);
    zbs1 = zeros(2*ntor+1,1);
    rbc1(ntor+1) = mean(R1c);
    rbs1(ntor+1) = mean(R1s);
    zbc1(ntor+1) = mean(Z1c);
    zbs1(ntor+1) = mean(Z1s);
%     for n = 1:ntor
%         % RBC
%         half_sum = mean(R1c .* cos(nfp*n*phi));
%         half_difference = mean(R1s .* sin(nfp*n*phi));
%         rbc1(ntor+1+n) = half_sum + half_difference;
%         rbc1(ntor+1-n) = half_sum - half_difference;
%         
%         % ZBC
%         half_sum = mean(Z1c .* cos(nfp*n*phi));
%         half_difference = mean(Z1s .* sin(nfp*n*phi));
%         zbc1(ntor+1+n) = half_sum + half_difference;
%         zbc1(ntor+1-n) = half_sum - half_difference;
%         
%         % RBS
%         half_sum = mean(R1s .* cos(nfp*n*phi));
%         half_difference = -mean(R1c .* sin(nfp*n*phi));
%         rbs1(ntor+1+n) = half_sum + half_difference;
%         rbs1(ntor+1-n) = half_sum - half_difference;
%         
%         % ZBS
%         half_sum = mean(Z1s .* cos(nfp*n*phi));
%         half_difference = -mean(Z1c .* sin(nfp*n*phi));
%         zbs1(ntor+1+n) = half_sum + half_difference;
%         zbs1(ntor+1-n) = half_sum - half_difference;
%     end

%%% change Katia 
    for n = 1:ntor
        % RBC
        half_sum = mean(R1c .* cos(nfp*n*phi));
        half_difference = mean(R1s .* sin(nfp*n*phi));
        rbc1(ntor+1+n) = half_sum + half_difference;
        rbc1(ntor+1-n) = half_sum - half_difference;
        
        % ZBC
        half_sum = mean(Z1c .* cos(n*nfp*phi));
        half_difference = mean(Z1s .* sin(n*nfp*phi));
        zbc1(ntor+1+n) = half_sum + half_difference;
        zbc1(ntor+1-n) = half_sum - half_difference;
        
        % RBS
        half_sum = mean(R1s .* cos(nfp*n*phi));
        half_difference = -mean(R1c .* sin(nfp*n*phi));
        rbs1(ntor+1+n) = half_sum + half_difference;
        rbs1(ntor+1-n) = half_sum - half_difference;
        
        % ZBS
        half_sum = mean(Z1s .* cos(nfp*n*phi));
        half_difference = -mean(Z1c .* sin(nfp*n*phi));
        zbs1(ntor+1+n) = half_sum + half_difference;
        zbs1(ntor+1-n) = half_sum - half_difference;
    end
    
%%%%%%%    
    
    rbc1 = rbc1 * r;
    rbs1 = rbs1 * r;
    zbc1 = zbc1 * r;
    zbs1 = zbs1 * r;
    assignin('base','rbc1',rbc1)
    assignin('base','rbs1',rbs1)
    assignin('base','zbc1',zbc1)
    assignin('base','zbs1',zbs1)
    
    % Test that I got the transformation right:
    for j = 1:N_phi
        surface_R = zeros(size(theta)) + R0(j);
        surface_Z = zeros(size(theta)) + Z0(j);
        m=1;
        for n = (-ntor):ntor
            angle = m*theta - n*phi(j); %change Katia
            sinangle = sin(angle);
            cosangle = cos(angle);
            surface_R = surface_R + rbc1(n+ntor+1)*cosangle + rbs1(n+ntor+1)*sinangle;
            surface_Z = surface_Z + zbc1(n+ntor+1)*cosangle + zbs1(n+ntor+1)*sinangle;
        end
        plot(surface_R,surface_Z,'o','Color',colors(j,:))
        plot(surface_R(1),surface_Z(1),'s','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',10)
    end
else
    % Handle the finite_r_option = 'nonlinear' case
    
    % New version of the transformation from Frenent to cylindrical
    % coordinates, which includes higher-order corrections in an attempt to
    % get smoother shapes, particularly for QHS.
    %N_phi = 30;
    theta_1D = linspace(0, 2*pi, N_theta+1);
    theta_1D(end) = [];
    %phi_1D = linspace(0, 2*pi/nfp, N_phi+1);
    %phi_1D(end) = [];
    phi_1D = phi;
    [phi_2D, theta_2D] = meshgrid(phi_1D, theta_1D);
    
    R_2D = zeros(size(phi_2D));
    x_2D = zeros(size(phi_2D));
    y_2D = zeros(size(phi_2D));
    z_2D = zeros(size(phi_2D));
    interp1_method = 'spline';
    
    %zeta_tripled = [zeta-2*pi/nfp; zeta; zeta+2*pi/nfp];
%     X1c_tripled = repmat(X1c_signed,3,1);
%     X1s_tripled = repmat(X1s_signed,3,1);
%     Y1c_tripled = repmat(Y1c_signed,3,1);
%     Y1s_tripled = repmat(Y1s_signed,3,1);

    X1c_tripled = repmat(-X1c_untwisted,3,1);
    X1s_tripled = repmat(-X1s_untwisted,3,1);
    Y1c_tripled = repmat(-Y1c_untwisted,3,1);
    Y1s_tripled = repmat(-Y1s_untwisted,3,1);

    
    index_for_sign_change =  [N_phi+1:2*N_phi];%indices_for_sign_flip + N_phi;
    
    X1c_tripled(index_for_sign_change) = -X1c_tripled(index_for_sign_change);
    X1s_tripled(index_for_sign_change) = -X1s_tripled(index_for_sign_change);
    Y1c_tripled(index_for_sign_change) = -Y1c_tripled(index_for_sign_change);
    Y1s_tripled(index_for_sign_change) = -Y1s_tripled(index_for_sign_change);
    
    figure(6)
    plot(X1c_tripled)
    hold on 
    plot(Y1c_tripled)
    plot(Y1s_tripled)
    plot(X1s_tripled)
    
    fprintf('Beginning finite_r_option nonlinear root-finding.\n')
    for j_theta = 1:N_theta
        costheta = cos(theta_1D(j_theta));
        sintheta = sin(theta_1D(j_theta));
        X1 = X1c_tripled * costheta + X1s_tripled * sintheta;
        Y1 = Y1c_tripled * costheta + Y1s_tripled * sintheta;
        

        for j_phi = 1:N_phi
            % Solve for the phi0 such that r0 + X n + Y b has the desired
            % phi:
            phi0_guess = phi_1D(j_phi);
            x0 = [phi0_guess - 1, phi0_guess + 1];
            phi0_solution = fzero(@phi_residual, x0);
            %[R_2D(j_theta,j_phi), x_2D(j_theta,j_phi), y_2D(j_theta,j_phi), z_2D(j_theta,j_phi)] = Rxyz(zeta_solution);
            this_Rxyz = Rxyz(phi0_solution);
            R_2D(j_theta,j_phi) = this_Rxyz(1);
            x_2D(j_theta,j_phi) = this_Rxyz(2);
            y_2D(j_theta,j_phi) = this_Rxyz(3);
            z_2D(j_theta,j_phi) = this_Rxyz(4);
        end
        
        
        
    end
    fprintf('Done with finite_r_option nonlinear root-finding.\n')
    
    %mpol = floor(N_theta/2 - 1);
    
    %ntor = floor(N_phi/3 - 1);
    mpol = floor(N_theta/3 - 1);
    
    rbc = zeros(2*ntor+1, mpol+1);
    rbs = zeros(2*ntor+1, mpol+1);
    zbc = zeros(2*ntor+1, mpol+1);
    zbs = zeros(2*ntor+1, mpol+1);
    
    % Handle the m=0 n=0 modes:
    rbc(ntor+1, 1) = mean(mean(R_2D));
    zbc(ntor+1, 1) = mean(mean(z_2D));
    
    % Handle the other modes:
    %{
d_theta = 2*pi / N_theta;
d_zeta = nfp*2*pi / N_phi;
normalizer = 2 * d_theta * d*zeta / (4*pi*pi);
    %}
    normalizer = 2 /(N_phi * N_theta);
    for m = 0:mpol
        if m==0
            n_min = 1;
        else
            n_min = -ntor;
        end
        for n = n_min:ntor
            angle = m * theta_2D - n*nfp*phi_2D; %Katia change
            sinangle = sin(angle);
            cosangle = cos(angle);
            rbc(ntor+1+n, m+1) = sum(sum(R_2D .* cosangle)) * normalizer;
            rbs(ntor+1+n, m+1) = sum(sum(R_2D .* sinangle)) * normalizer;
            zbc(ntor+1+n, m+1) = sum(sum(z_2D .* cosangle)) * normalizer;
            zbs(ntor+1+n, m+1) = sum(sum(z_2D .* sinangle)) * normalizer;
        end
    end

    % Test that I got the transformation right:
    for j = 1:N_phi
        surface_R = zeros(size(theta_1D));
        surface_Z = zeros(size(theta_1D));
        for m = 0:mpol
            for n = (-ntor):ntor
                angle = m*theta_1D - n*nfp*phi(j);
                sinangle = sin(angle);
                cosangle = cos(angle);
                surface_R = surface_R + rbc(n+ntor+1, m+1)*cosangle + rbs(n+ntor+1, m+1)*sinangle;
                surface_Z = surface_Z + zbc(n+ntor+1, m+1)*cosangle + zbs(n+ntor+1, m+1)*sinangle;
            end
        end
        plot(surface_R,surface_Z,'o-','Color',colors(j,:))
    end
end


figure(987)

plot(varphi,B0)
hold on 
plot(varphi, 1+0.1*cos(varphi))

figure(989) 

theta = linspace(0, 2*pi, 100);
theta1 = theta;
%theta = -(theta1 + pi);
B1 = (B0.*B1c_over_B0)*cos(theta) + (B0.*B1s_over_B0)*sin(theta);
%size(varphi)
%size(theta)
%size(B1)

colormap jet

contourf(varphi,theta1,B1',23)
colorbar
title(['B_{1}, A = ',num2str(aspect_ratio)])
xlabel('Boozer \phi')
ylabel('Boozer \theta')


figure(990)

B_tot = B0 + r*B1;
contour(varphi,theta1,B_tot',23,':','LineWidth',1.5)
colorbar
title(['Total B, A =',num2str(aspect_ratio)])
xlabel('Boozer \phi')
ylabel('Boozer \theta')

figure(991) 

%theta = linspace(0, 2*pi, 601);
% B1 = B0.*B1c_over_B0*cos(theta) + B0.*B1s_over_B0*sin(theta);
B1 = r*B1;

contour(varphi,theta1,B1',23,':','LineWidth',1.5)
colorbar
title(['B_{1}, A = ',num2str(aspect_ratio)])
xlabel('Boozer \phi')
ylabel('Boozer \theta')

B_filename = ['B_B1_Objective_',num2str(nfp),'FP_A',num2str(aspect_ratio),'_2.mat'];
B0_O = B0;
save(B_filename,'varphi','B0_O','B1s_over_B0','B1c_over_B0','r');

%theta_testVMEC = 0;

 theta_testVMEC = [0,pi/7,pi/6,3*pi/2];
% %
% R1c_VMEC = zeros(length(phi_Complete),length(theta_testVMEC));
% 
% 
% for jj= 1:length(theta_testVMEC)
% 
%     theta_test = theta_testVMEC(jj);
%     
%     for j = 1:numel(R0c)
%         R1c_VMEC(:,jj) = R1c_VMEC(:,jj) + R0c(j)*cos(-(j-1)*phi_Complete*nfp);
%     end
% 
%     for n=(-ntor):ntor
%         R1c_VMEC(:,jj) = R1c_VMEC(:,jj) + rbc1(ntor+1+n)*cos(theta_test-n*nfp*phi_Complete);
%     end
% 
% %     for m = 0:10
% %             if m==0
% %                 n_min = 0;
% %             else
% %                 n_min = -ntor;
% %             end
% %             for n = n_min:ntor
% %                 R1c_VMEC(:,jj) = R1c_VMEC(:,jj) + rbc(ntor+1+m,m+1) *cos(m*theta_test - n*nfp*phi_Complete);
% %             end
% %     end
% 
% 
% end 
% 
% 
figure(531)
% plot(phi_Complete, R1c_VMEC(:,3))
% hold on 
plot(phi,R0+r*(R1c.*cos(theta_testVMEC(2))+R1s.*sin(theta_testVMEC(2))))

figure(532)
plot(phi,Z0+r*(Z1c*cos(theta_testVMEC(2))+Z1s*sin(theta_testVMEC(2))))


    figure(4565)
    
    plot(varphi,alpha_notIota+alpha_iota*iota)
    hold on 
    alpha_omnigenous = alpha_omnigenous_iota*iota + alpha_omnigenous_notIota;
    plot(varphi,alpha_omnigenous)
    
    alpha_RMS = sqrt(mean(alpha-alpha_omnigenous).^2);

    if strcmp(finite_r_option,'nonlinear')
        rbc1 = rbc(:,1);
        rbs1 = rbs(:,1);
        zbs1 = zbs(:,1);
        zbc1 = zbc(:,1);
    end
    

% NOTE: For now, this next part only works for stellarator symmetry!

vmec_filename = [vmec_filename_base,'_finite_r_',finite_r_option,'_A',num2str(aspect_ratio),'_delta',num2str(delta)];
if alpha_type == 1
    alpha_type_name = 'Buffer region approach, see d';
else 
    alpha_type_name = 'Non-zone approach, see k';
end
lasym = (any(R0s ~= 0) || any(Z0c ~= 0));
fid = fopen(vmec_filename,'w');
fprintf(fid,'! Generated by %s.m\n',mfilename);
fprintf(fid,['! on ',datestr(now,'yyyymmdd_HH_MM_SS'),'.\n']);
fprintf(fid,'! N_phi = %d\n',N_phi);
fprintf(fid,'! alpha_type = %s\n',alpha_type_name);
fprintf(fid,'! delta = %g\n',delta);
fprintf(fid,'! k = %g\n',k_for_alpha);
fprintf(fid,'! B0_1 = %g\n',B0_1);
fprintf(fid,'! aspect_ratio = %g\n\n\n',aspect_ratio);
fprintf(fid,'\n');
fprintf(fid,'&INDATA\n');
if (lasym)
    fprintf(fid,'LASYM = T\n');
else
    fprintf(fid,'LASYM = F\n');
end
fprintf(fid,'DELT = 0.9\n');
fprintf(fid,'TCON0 = 2.\n');
fprintf(fid,'NFP = %d\n',nfp);
fprintf(fid,'NCURR = 1\n');
fprintf(fid,'MPOL = 9\n');
fprintf(fid,'NTOR = 80\n');
fprintf(fid,'NS_ARRAY    =    13    25    51   101\n');
fprintf(fid,'NITER_ARRAY =  1000  1500  2000 20000\n');
% % fprintf(fid,'FTOL_ARRAY  = 1e-30 1e-30 1e-30 1e-16\n');
fprintf(fid,'FTOL_ARRAY  = 1e-17 1e-17 1e-15 1e-13\n');
fprintf(fid,'NSTEP = 200\n');
fprintf(fid,'GAMMA = 0\n');
fprintf(fid,'SPRES_PED = 1.0\n');
% % % fprintf(fid,'PHIEDGE = 1\n');
fprintf(fid,'PHIEDGE = %22.15e\n',pi*r*r);
fprintf(fid,'CURTOR = 0\n');
fprintf(fid,'AM = 0\n');
fprintf(fid,'AC = 0\n');
fprintf(fid,'AI = 0\n');
fprintf(fid,'RAXIS_CC = ');
for j = 1:numel(R0c)
    fprintf(fid,'%22.15e ',R0c(j));
end
fprintf(fid,'\nZAXIS_CS = ');
for j = 1:numel(R0c)
    fprintf(fid,'%22.15e ',Z0s(j)); %%%SignChange % Why is the minus sign here?
end
fprintf(fid,'\n');
if (lasym)
    fprintf(fid,'RAXIS_CS = ');
    for j = 1:numel(R0c)
        fprintf(fid,'%22.15e ',R0s(j));
    end
    fprintf(fid,'\nZAXIS_CC = ');
    for j = 1:numel(R0c)
        fprintf(fid,'%22.15e ',Z0c(j)); %%%SignChange
    end
    fprintf(fid,'\n');
end
if strcmp(finite_r_option,'linear')
    for j=1:numel(R0c)
        fprintf(fid,'RBC(%3d,0) = %22.15e   ZBS(%3d,0) = %22.15e\n',j-1,R0c(j),j-1,Z0s(j)); %%%SignChange % The minus sign is here for Z0s because vmec modes are ~sin(mu-nv), which for m=0, is sin(-nv)=-sin(nv).
    end
    for n=(-ntor):ntor
        fprintf(fid,'RBC(%3d,1) = %22.15e   ZBS(%3d,1) = %22.15e\n',n,rbc1(ntor+1+n),n,-zbs1(ntor+1+n));%%%SignChange
    end
    if (lasym)
        fprintf(fid,'\n');
        %fprintf(fid,'RBS(  0,0) = 0.0   ZBC(  0,0) = 0.0\n');
        %fprintf(fid,'RBS(  1,0) = %22.15e   ZBC(  1,0) = %22.15e\n',R0s,Z0c); % Should R0s have a minus sign here like Z0s did?
        for j=1:numel(R0c)
            fprintf(fid,'RBS(%3d,0) = %22.15e   ZBC(%3d,0) = %22.15e\n',j-1,R0s(j),j-1,Z0c(j)); %%%SignChange % Should R0s have a minus sign here like Z0s did?
        end
        for n=(-ntor):ntor
            fprintf(fid,'RBS(%3d,1) = %22.15e   ZBC(%3d,1) = %22.15e\n',n,-rbs1(ntor+1+n),n,zbc1(ntor+1+n));%%%SignChange
        end
    end
else
    for m = 0:mpol
        if m==0
            n_min = 0;
        else
            n_min = -ntor;
        end
        for n = n_min:ntor
            fprintf(fid,'RBC(%3d,%3d) = %22.15e   ZBS(%3d,%3d) = %22.15e\n',n,m,rbc(ntor+1+n,m+1),n,m,-zbs(ntor+1+n,m+1));%%%SignChange
            if (lasym)
                fprintf(fid,'RBS(%3d,%3d) = %22.15e   ZBC(%3d,%3d) = %22.15e\n',n,m,-rbs(ntor+1+n,m+1),n,m,zbc(ntor+1+n,m+1));%%%SignChange
            end
        end
    end
end
fprintf(fid,'\n');
fprintf(fid,'/\n');

fprintf(fid,'\n');
fprintf(fid,'&optimum\n');
fprintf(fid,'NOPTIMIZERS=1\n');
fprintf(fid,'NFUNC_MAX = 1\n');
fprintf(fid,'EQUIL_TYPE = "VMEC2000"\n');
fprintf(fid,'OPT_TYPE   = "one_iter"\n');
fprintf(fid,'FTOL =  1.00000000000000E-06\n');
fprintf(fid,'XTOL =  1.00000000000000E-30\n');
fprintf(fid,'GTOL =  1.00000000000000E-30\n');
fprintf(fid,'FACTOR =   100.0\n');
fprintf(fid,'EPSFCN =   1.0E-05\n');
fprintf(fid,'LKEEP_MINS = T\n');
fprintf(fid,'\n');
fprintf(fid,'MBOZ = 35\n');
fprintf(fid,'NBOZ = 160\n');
fprintf(fid,'TARGET_NEO(2:10) = 9*0.0  SIGMA_NEO(2:10) = 9*1.0\n');
fprintf(fid,'TARGET_NEO(21) = 0.0  SIGMA_NEO(21) = 1.0\n');
fprintf(fid,'TARGET_NEO(31) = 0.0  SIGMA_NEO(31) = 1.0\n');
fprintf(fid,'TARGET_NEO(41) = 0.0  SIGMA_NEO(41) = 1.0\n');
fprintf(fid,'TARGET_NEO(51) = 0.0  SIGMA_NEO(51) = 1.0\n');
fprintf(fid,'TARGET_NEO(61) = 0.0  SIGMA_NEO(61) = 1.0\n');
fprintf(fid,'TARGET_NEO(71) = 0.0  SIGMA_NEO(71) = 1.0\n');
fprintf(fid,'TARGET_NEO(81) = 0.0  SIGMA_NEO(81) = 1.0\n');
fprintf(fid,'TARGET_NEO(91) = 0.0  SIGMA_NEO(91) = 1.0\n');
fprintf(fid,'TARGET_NEO(101) = 0.0  SIGMA_NEO(101) = 1.0\n');
fprintf(fid,'/\n\n');
fprintf(fid,'&NEO_IN\n');
fprintf(fid,'THETA_N = 100\n');
fprintf(fid,'PHI_N = 100\n');
fprintf(fid,'MAX_M_MODE = 0\n');
fprintf(fid,'MAX_N_MODE = 0\n');
fprintf(fid,'NPART = 75\n');
fprintf(fid,'MULTRA = 1\n');
fprintf(fid,'ACC_req = 0.01\n');
fprintf(fid,'NO_BINS = 100\n');
fprintf(fid,'NSTEP_PER = 150\n');
fprintf(fid,'NSTEP_MIN = 500\n');
fprintf(fid,'NSTEP_MAX = 200000\n');
fprintf(fid,'CALC_NSTEP_MAX = 0\n');
fprintf(fid,'EOUT_SWI = 1\n');
fprintf(fid,'LAB_SWI = 0\n');
fprintf(fid,'INP_SWI = 0\n');
fprintf(fid,'REF_SWI = 2\n');
fprintf(fid,'WRITE_PROGRESS = 1\n');
fprintf(fid,'WRITE_OUTPUT_FILES = 0\n');
fprintf(fid,'SPLINE_TEST = 0\n');
fprintf(fid,'WRITE_INTEGRATE = 0\n');
fprintf(fid,'WRITE_DIAGNOSTIC = 0\n');
fprintf(fid,'CALC_CUR = 0\n');
fprintf(fid,'NPART_CUR = 200\n');
fprintf(fid,'/\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fclose(fid);


%****************************************************
% 3D Figure of surface, adding X n + Y b to the axis
%****************************************************

if strcmp(finite_r_option,'linear')
    return
end

figure(51 + figure_offset)
clf

% Close the surface in theta:
x_2D_closed = [x_2D; x_2D(1,:)];
y_2D_closed = [y_2D; y_2D(1,:)];
z_2D_closed = [z_2D; z_2D(1,:)];

assignin('base','x_2D_m',x_2D)
assignin('base','y_2D_m',y_2D)
assignin('base','z_2D_m',z_2D)
assignin('base','R_2D_m',(x_2D.^2 + y_2D.^2))

g = 0.8;
surf(x_2D_closed, y_2D_closed, z_2D_closed, 'edgecolor','r','FaceColor',[g,g,g],'DisplayName','const \phi')
daspect([1,1,1])
axis vis3d
set(gca,'clipping','off')
%light
zoom(1.2)


%N_theta = 20;
%N_phi = 30;
theta_1D = linspace(0, 2*pi, N_theta+1);
%theta_1D = linspace(0, 2*pi, N_theta+1);
theta_1D(end) = [];
%phi_1D = linspace(0, 2*pi/nfp, N_phi+1);
%phi_1D(end) = [];
phi_1D = phi;
[phi_2D, theta_2D] = meshgrid(phi_1D, theta_1D);

x_2D = zeros(size(phi_2D));
y_2D = zeros(size(phi_2D));
z_2D = zeros(size(phi_2D));

x0 = R0(1:N_phi) .* cos(phi); %%% changes Katia
y0 = R0(1:N_phi) .* sin(phi);

for j_theta = 1:N_theta
    costheta = cos(theta_1D(j_theta));
    sintheta = sin(theta_1D(j_theta));
    X1 = X1c_signed * costheta + X1s_signed * sintheta;
    Y1 = Y1c_signed * costheta + Y1s_signed * sintheta;
    
    x_2D(j_theta,:) = x0 + r * (X1 .* kappa_Cartesian(:,1) + Y1 .* binormal_Cartesian(:,1));
    y_2D(j_theta,:) = y0 + r * (X1 .* kappa_Cartesian(:,2) + Y1 .* binormal_Cartesian(:,2));
    z_2D(j_theta,:) = Z0 + r * (X1 .* kappa_Cartesian(:,3) + Y1 .* binormal_Cartesian(:,3));    
end

hold on
surf(x_2D, y_2D, z_2D, 'edgecolor','k','FaceColor','none','DisplayName','Original Frenet')

legend show

end