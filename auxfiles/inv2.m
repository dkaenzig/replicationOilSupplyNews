function  InvA  = inv2(A)
% Use pseudo-inverse
temp=eye(size(A,2));
InvA=A\temp;

end

