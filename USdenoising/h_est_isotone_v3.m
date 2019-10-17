function h=h_est_isotone_v3(x,M,AV)

%  Estimates the variance function h.
%   
%    takes:
%    x - data
%    M - half-span in kernel smoothing
%    method - either "isotone" (plain isotone regression) or 
%              "iso.correction" (kernel smoothing and then isotone regression
%                                of the smooth)
%    
%    returns:
%    list containing the domain points u and the values of h for those points


% Younes, 02/02/2016
% v3, blind, 15/05/2016


%filtered_image = medfilt2(x,[K K]);


filtered_image = AV;
alpha_hat = filtered_image(:);

s=51200;
s=s/200;
  
  
[domain_points, order_domain_points] = sort(alpha_hat);
  
  
 eps_hat_2 = (x(:) - alpha_hat).^2;
  
 
 h_hat = ksr(domain_points, eps_hat_2(order_domain_points),2*M+1, s);
 
 
 

 h.dot = lsqisotonic([1:s], h_hat.f, ones(1,s));
  
 h.u = h_hat.x;
  
end
