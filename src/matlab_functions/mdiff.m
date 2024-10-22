function c=mdiff(a,b);


% c=mdiff(a,b) forms matrix c(n,m)=a(n)-b(m)
if(nargin==1)
b=a;
end
c=repmat(a(:),1,length(b(:)))-repmat(b(:).',length(a(:)),1);