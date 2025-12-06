%%input Data muliti-columns are OK, each column represents a time series
%%input method =1 use Toeplitz matrix; =2 use Trajectory matrix

% Tested on 9.13.0.2049777 (R2022b)
% Last modified by ypf0323@163.com, 12/18/2024
function Covariance=GrandMatrix(Data,M,Method)
   [N,L]=size(Data);
   N1=N-M+1;
   
   if(Method==1)
       for l1=1:L
           for l2=1:L
               for j1=1:M
                   for j2=1:M
                       tmin=max(1,1+j1-j2);
                       tmax=min(N,N+j1-j2);
                       N12=tmax-tmin+1;
                       dj=j2-j1;
                       xl1=Data(tmin:tmax,l1);
                       xl2=Data(tmin+dj:tmax+dj,l2);
                       T((l1-1)*M+j1,(l2-1)*M+j2)=sum(xl1.*xl2)/N12;
                   end
               end
           end
       end
       Covariance=T;
   end
   
   if(Method==2)
       for col=1:L
           for row=1:N1
               d0(row,:)=Data(row:row+M-1,col);
           end
           if(col==1)
               d=d0;
           else
               d=[d d0];
           end
       end
       Covariance=1/N1*d'*d;
   end
   
end