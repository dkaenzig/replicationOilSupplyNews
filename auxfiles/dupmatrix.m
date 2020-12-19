function D = dupmatrix(n)

% QIP.OPEN_SYSTEMS.DUPLICATION   Generated a duplication matrix.
% requires: nothing
% author: Marcus da Silva
%
%    M = DUPLICATION(D) Generates a duplication matrix acting on the result of 
%    an elimination of a vectorization of a D by D matrix. When M multiplies such a vector, 
%    the result is a larger vector consisting of all parts of the symmetric matrix, 
%    in vectorized form.
%
%    See also: vec, row, rowinv, reshape
%
%   Copyright (C) 2010   Marcus P da Silva http://github.com/marcusps
% 
%   License: Distributed under Apache License, Version 2.0
%            http://www.apache.org/licenses/LICENSE-2.0
%

%    Copyright 2012 Marcus P. da Silva
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

vechindex = @(r,c,n) n * r + c - r*(r+1)/2;

% duplication
D = zeros(n*n,n*(n+1)/2);
for j=1:size(D,1)
  rowo = floor((j-1)/n);
  colo = mod(j-1,n);
  if (rowo <= colo)
    D(j,vechindex(rowo,colo,n)+1) = 1;
  else
    D(j,vechindex(colo,rowo,n)+1) = 1;
  end
end



