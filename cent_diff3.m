function [ xp ] = cent_diff3( x, dt )
%CENT_DIFF3 Numerical approximation of the derivative of x via a
%three-point central difference which provides second order accuracy (i.e.,
%truncation errors on the order of dt^2). Derivatives of the first and last
%data points are estimated using front_diff and back_diff, respectively.
%
% ARGUMENTS
% x - Vector or matrix containing numerical data. If x is a matrix
% differentation operates along the vertical dimensions (i.e., each column
% is differentiated).
% 
% dt - Unit change in time between samples (assumed constant).
%
% RETURNS
% xp - First-order derivative of x with respect to time. If x is an m x n
% vector or matrix, xp will have a size of m x n.
%
%=========================================================================%

[ r, c ] = size(x); %determine the size of x

xp = zeros(r, c); %create empty matrix for the derivative values

%for each column
for i = 1:c
    
    %for each row
    for j = 1:r
        
        if j == 1 %if at the first data point take forward difference
            
            xp(j,i) = front_diff( x(j:j+1,i), dt);
            
        elseif j == r % if at the last data point take backward difference
            
            xp(j,i) = back_diff( x(j-1:j,i), dt);
            
        else
                
            xp(j,i) = (x(j+1,i) - x(j-1,i)) / (2*dt); %calculate three-point central difference
            
        end %end if
        
    end %end row
    
end %end column

end %end function

