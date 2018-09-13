function [ sigma ] = Nye2Tensor( s )
%Nye2Tensor     converts a stress vector in Nye notation to tensor notation
%Input:
%  <s>          stress vector in Nye notation [sx,sy,sz,(sxy,sxz,syz)]
%Output:
%  <sigma>      stress tensor [sx,sxy,sxz;sxy,sy,syz;sxz,syz,sz]
%CL

sigma = [s(1),s(4),s(5);
         s(4),s(2),s(6);
         s(5),s(6),s(3)];

end

