function Fisz_est  = ddhf_est_Iso_Kronecker_soft_unk(x, AV, zero_levels, thr, theta, BseOnd)

  
% %    Finds the wavelet-Fisz estimate using some given (BseOnd) wavelet.
% %    
% %    takes:
% %    x           - data
% %    AV          - a pre-estimate of the fdistribution of pixel values. Can be an
% %                  average overt time of the dynamic sequence.
% %    zero.levels - number of finest-scale levels at which detail coefficients are set to zero
% %    thr         - type of threshold: 1 for universal, 2 set by user
% %    theta       - multiplier for the (universal) threshold (use 1 if not sure)
% %    zero.levels - number of finest-scale levels at which detail coefficients are set to zero
% %    BseOnd      - wavelet type. 

% %    returns:
% %    the estimate "Fisz_est"
  
% % Younes, 12/04/2017



  n=size(x,1);
  s = size(x,2);
  J = log2(s);



    if nargin<6
     BseOnd = 'db6';
    end

    
     if nargin<5
     theta = 1;
     end
    
     if nargin<4
     thr = 1;
     end
     
     if nargin<3
     zero_levels = 1;
     end

 
     
     
     Fisz_est =zeros(size(x));
     
     
     
    % Stationary Isotropic Wavelet transform

    
    
A_all = zeros(n,s,s,J);
H_all = zeros(n,s,s,J);
D_all = zeros(n,s,s,J);
V_all = zeros(n,s,s,J);
    

 for i=1:n
    
    
    [A,H,V,D]  = swt2(squeeze(x(i,:,:)),J,BseOnd);
    
   
    
    
    % local means
    loc_mean = zeros(size(A,1),size(A,2));
    loc_mean_list = zeros(size(A));
    %Fisz_normal = zeros(size(A));
    
    j=.5*[1:J];
    j = 2.^j;
    
    
    
    h=h_est_isotone_v3(squeeze(x(i,:,:)),1,AV);
    
    for k=1:J
        
        B=A(:,:,k);
        
        loc_mean = squeeze(B/j(k));
        
        C=loc_mean(:);
        
        fun = function_from_vector(h.u, h.dot, C');
        
        Fisz(:,:,k) = reshape(fun,s,s);
        
        loc_mean_list(:,:,k) = loc_mean;
        
 
    end

        
    
   %% Fisz normalization
   
   
   POS_TEST = loc_mean_list > 0;
   
   
    Fisz_normal = sqrt(Fisz(POS_TEST));
    
    H(POS_TEST) = H(POS_TEST)./Fisz_normal;
    D(POS_TEST) = D(POS_TEST)./Fisz_normal;
    V(POS_TEST) = V(POS_TEST)./Fisz_normal;
    
    
    
    
    % Set to zero fine scales
% %     
    if (zero_levels); 	
        
        H(:,:,1:zero_levels) = 0; 
        D(:,:,1:zero_levels) = 0; 
        V(:,:,1:zero_levels) = 0; 
    
    end


    
       
    
    Fisz_normal_all{i}=Fisz_normal;
    A_all(i,:,:,:)=A;
    D_all(i,:,:,:)=D;
    H_all(i,:,:,:)=H;
    V_all(i,:,:,:)=V;
    
    POS_TEST_all{i}=POS_TEST;
    
    clear Fisz_normal
 end
    
    
    
    Hyperbolic_H = swt(H_all,5,'haar');
    Hyperbolic_D = swt(D_all,5,'haar');
    Hyperbolic_V = swt(V_all,5,'haar');
    
    
    
    
    if thr == 1
    
    T=sqrt(2*J*log10(2));
    
        else if thr == 2

                %T = sqrt(2*log(n*n)); 
                T=3;
            end
    end
    
    T=theta*T;
    
    
    
    
    
    Hyperbolic_H = Hyperbolic_H .* (abs(Hyperbolic_H) > T);
    Hyperbolic_D = Hyperbolic_D .* (abs(Hyperbolic_D) > T);
    Hyperbolic_V = Hyperbolic_V .* (abs(Hyperbolic_V) > T);
    
  %   Hyperbolic_H = max(abs(Hyperbolic_H) - T, 0);
   %  Hyperbolic_D = max(abs(Hyperbolic_D) - T, 0);
   %  Hyperbolic_V = max(abs(Hyperbolic_V) - T, 0);
    
    I_Hyperbolic_V = iswt(Hyperbolic_V,'haar');
    I_Hyperbolic_H = iswt(Hyperbolic_H,'haar');
    I_Hyperbolic_D = iswt(Hyperbolic_D,'haar');
    
    
        
    %threshold 
    %T=sqrt(2*log2((J-zero_levels)*n*n*3));
    
    %T=sqrt(7*2*log(2))^2;
    
    %T=sqrt(2*log(J-zero_levels)*n*n*3);

    
    H_all = reshape(I_Hyperbolic_H,n,s,s,J);
    D_all = reshape(I_Hyperbolic_D,n,s,s,J);
    V_all = reshape(I_Hyperbolic_V,n,s,s,J);
    
    
    
    
    
    
    
    for i=1:n
       
        
    Fisz_normal=squeeze(Fisz_normal_all{i});
    A=squeeze(A_all(i,:,:,:));
    D=squeeze(D_all(i,:,:,:));
    H=squeeze(H_all(i,:,:,:));
    V=squeeze(V_all(i,:,:,:));
        
    POS_TEST=squeeze(POS_TEST_all{i});
    
    D(POS_TEST)=D(POS_TEST).*Fisz_normal;
    H(POS_TEST)=H(POS_TEST).*Fisz_normal;
    V(POS_TEST)=V(POS_TEST).*Fisz_normal;
    
    Fisz_est(i,:,:) = squeeze(iswt2(A,H,V,D,BseOnd));
    
    end
 
    
end
