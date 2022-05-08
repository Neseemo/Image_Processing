function varargout = show_dictionary(D)
% SHOW_DICTIONARY(D) plot the 2d atoms in the dictionary D. The atoms are
% column vector obtained by vectorizing square 2d atoms

% for i=1:(size(D,2)-1)
%     v = sum((bsxfun(@minus,D(:,i),D(:,i+1:end))).^2,1);
%     [~,b] = min(v);
%     temp = D(:,i+1);
%     D(:,i+1) = D(:,b+i);
%     D(:,b+i) = temp;
% end

bound = 2;

psz = sqrt(size(D,1));

if (psz~=round(psz))
    error('Atoms must be squared');
end

np = size(D,2);
M = ceil(sqrt(np));
N = ceil(np/M);

img = ones(M*psz+bound*(M-1),N*psz+bound*(N-1));

for i=1:np
    m = mod(i,M);
    if (m==0)
        m = M;
    end
    n = (i-m)/M+1;
    
    m = (m-1)*psz + bound*(m-1) + 1;
    n = (n-1)*psz + bound*(n-1) + 1;
    
    atom =  reshape(D(:,i),psz,psz);
    if (min(atom(:)) ~= max(atom(:)))
        atom = (atom-min(atom(:)))/(max(atom(:))-min(atom(:)));
    end
    img(m:m+psz-1,n:n+psz-1) = atom;
end

imagesc(img);
colormap(gray);

if (nargout == 1)
    varargout{1} = img;
end






