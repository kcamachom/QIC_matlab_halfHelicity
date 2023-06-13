%%%%%%%%%%
%%%%%%%%%% This code scans d coefficients to find curves with 
%%%%%%%%%% half helicity and low elongation 
%%%%%%%%%%

global nfp 
global max_elongation
global min_elongation
counter = 0;

nfp = 5;

dk_c = [0,0,0];
dk_s = [0,0,0];

for i= -0.05:0.1:0.05
    for j = -0.05:0.1:0.05
        for l = -0.05:0.1:0.05
           
        
        x0(1) = i; %dk_c(1)
        x0(2) = j; %dk_c(2)
        x0(3) = l; %dk_c(3) 
        %x0(4) = 0;
         
        options = optimset('TolX',1e-5,'TolFun',5.0e-2);
        [x,fval] = fminsearch(@NAEfunction,x0,options);
        
        if fval<= 5
            %if abs(Tg(end,1)-Tg(1,1)) <= 3e-3
            %if sqrt(sum((Tg(end,:)-Tg(1,:)).^2)) <= 1e-2
                counter = counter + 1;
                goodValues(counter,1) = x(1);
                goodValues(counter,2:3) = x(2:3);
                goodValues(counter,4) = 123;
                goodValues(counter,5) = fval;
                goodValues(counter,6) = max_elongation;
                goodValues(counter,7) = min_elongation;              
            %end
            
        end
        end
    end
end

function objective_function = NAEfunction(data)
%%%%%%%%% Frenet Frame Solver
global nfp
global max_elongation
global min_elongation 

dk = 0.4;%0.28;
dk_c(1) = data(1);
dk_c(2) = data(2);
dk_c(3) = data(3);

R0c= [1.0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
Z0s = [0,   0,  0,  0,  0,  0, 0, 0,0,0,0,0,0];
R0c(4) = -0.001;
Z0s(3) = 0.007;
Z0s(4) = 0.0026;
R0c(2) = -((R0c(4) + 9*nfp*nfp*R0c(4))/(1 + nfp*nfp));
R0c(3) = -(1/(1 + 4*nfp*nfp));
Z0s(2) = -((4*Z0s(3) + 8*nfp*nfp*Z0s(3) + 6*Z0s(4)+ 27*nfp*nfp*Z0s(4))/(2 + nfp^2));
R0s = zeros(size(R0c));
Z0c = zeros(size(R0c));
 

[varphi,phi,iota,alpha_RMS,curvature, torsion, n_for_alpha, sigma_max, elongation] = mGarrenBoozerOmnigenity_minimizer_halfHelicity(nfp,R0c,R0s,Z0c,Z0s,dk,dk_c);

    objective_function = (max(elongation) - min(elongation)) ;     
     %objective_function = (max(curvature)-3).^2;
    max_elongation = max(elongation);
    min_elongation = min(elongation);
    
%      figure(234)
%      plot(phi,elongation)
%      hold on 
%     plot(phi,torsion)

end


