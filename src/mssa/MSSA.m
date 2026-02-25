%%input: time, for different channels with same time period and interval
%data,each column represents a time series,multi-columns are OK
%num:how many eofs or pcs to be stored
% method: 1£ºWhen the number of data is not large£¬2£ºWhen the number of data
% is large
% num: calculate the number of data
%
% Tested on 9.13.0.2049777 (R2022b)
% Last modified by ypf0323@163.com, 12/18/2024
% Last modified by zhongtian.ma@connect.polyu.hk, 06/12/2025

function [evalues,st_eofs,st_pcs,RC]=MSSA(time,data,M,Method,N)                                                                                                
tic
[H,L]=size(data);                          
%%%Calculate grand covariance matrix 1=Toeplize Matrix 2=Trajectory Matrix
T=GrandMatrix(data,M,Method);%Ð­·½²î¾ØÕó

% Eigenvalues method 1: QZ algorithm
% [eigenvectors,eigenvalues]=eig(T);
% evalues=diag(eigenvalues);                                         
% [evalues sort_index]=sort(evalues,'descend');
% evalues=evalues/sum(evalues);
% evectors=eigenvectors(:,sort_index);
% st_eofs=evectors(:,1:num);

% Eigenvalues method 2: SVD (this is faster)
[U,S,V] = svd(T);
evalues=diag(S);   
[evalues, sort_index]=sort(evalues,'descend');
evalues=evalues/sum(evalues);

% keep the num according to its evalues
% if num<0 && num>=-1
%     
%     evalues_sum=zeros(numel(evalues),1);
%     for i=1:numel(evalues)
%         evalues_sum(i)=sum(evalues(1:i));
%     end
% 
%     kk=find(evalues_sum>=-num);
% 
%     num=kk(1);
% end

evectors=U(:,sort_index);    
st_eofs=evectors(:,1:N);

st_pcs=zeros(H-M+1,N);
for t=1:H-M+1
    for i=1:L
        X1=data(t:t+M-1,i);
        st_pcs(t,:)=st_pcs(t,:)+X1'*st_eofs((i-1)*M+1:i*M,:);
    end
end

RC=zeros(H,N,L);
for l=1:L
    for k=1:N
        for t=1:M-1
            index1=t:-1:1;
            index2=1:t;
            RC(t,k,l)=1/t*st_pcs(index1,k)'*st_eofs((l-1)*M+index2,k);
        end
        for t=M:H-M+1
            index1=t:-1:t-M+1;
            index2=1:M;
            RC(t,k,l)=1/M*st_pcs(index1,k)'*st_eofs((l-1)*M+index2,k);
        end
        for t=H-M+2:H
            index1=H-M+1:-1:t-M+1;
            index2=t-H+M:M;
            RC(t,k,l)=1/(H-t+1)*st_pcs(index1,k)'*st_eofs((l-1)*M+index2,k);
        end
    end
end
timeused=toc;
fprintf(1,'M-SSA (SVD) Computational time    = %f seconds\n', timeused);
end