Lx=8;
Ly=8;

N_sites=Lx*Ly;
N_up=32;
N_dwn=32;

H_0 = zeros(N_sites,N_sites);

H_0(1,2)=-1;
H_0(1,Lx+1)=-1;
H_0(1,Lx+2)=-1;

H_0(N_sites,N_sites-1)=-1;
H_0(N_sites,N_sites-Lx)=-1;
H_0(N_sites,N_sites-Lx-1)=-1;


for i = 2:N_sites-1

	if (i == Lx)
		H_0(i,i-1)=-1;
                H_0(i,i+Lx)=-1;
        elseif (i == (Lx-1)*Ly + 1)
		H_0(i,i+1)=-1;
		H_0(i,i-Lx)=-1;
	elseif (i < Lx)
		H_0(i,i-1)=-1;
		H_0(i,i+1)=-1;
		H_0(i,i+Lx)=-1;
		H_0(i,i+Lx+1)=-1;
	elseif (i > (Lx-1)*Ly+1)
		H_0(i,i-1)=-1;
                H_0(i,i+1)=-1;
                H_0(i,i-Lx)=-1;
                H_0(i,i-Lx-1)=-1;
	elseif (mod(i,Lx) == 0)
		H_0(i,i-1)=-1;
                H_0(i,i-Lx-1)=-1;
		H_0(i,i-Lx)=-1;
                H_0(i,i+Lx)=-1;
	elseif (mod(i,Lx) == 1)
		H_0(i,i+1)=-1;
                H_0(i,i-Lx)=-1;
                H_0(i,i+Lx)=-1;
                H_0(i,i+Lx+1)=-1;

	else
		H_0(i,i-Lx-1)=-1;
		H_0(i,i-Lx)=-1;
		H_0(i,i-1)=-1;
		H_0(i,i+1)=-1;
		H_0(i,i+Lx)=-1;
		H_0(i,i+Lx+1)=-1;

	endif

endfor

[V,lambda] = eig(H_0);

psi=horzcat(V(:,1:N_up),V(:,1:N_dwn));
psi_up=horzcat(V(:,1:N_up));
psi_dwn=horzcat(V(:,1:N_dwn));

nij=psi*psi';
nij_up=psi_up*psi_up';
nij_dwn=psi_dwn*psi_dwn';

dnidnj=zeros(N_sites,N_sites);

for i=1:N_sites-1
	for j=i+1:N_sites
        dnidnj(i,j)=(nij_up(i,i)+nij_dwn(i,i))*(nij_up(j,j)+nij_dwn(j,j));
	dnidnj(i,j)=dnidnj(i,j)-nij_up(i,j)*nij_up(j,i)-nij_dwn(i,j)*nij_dwn(j,i);
	dnidnj(i,j)=dnidnj(i,j)-(nij(i,i)+nij(j,j))+1;
	d = dnidnj(i,j);
	save -ascii -append file2 d
	endfor
endfor



for i=1:N_up
	for j=1:N_sites
		p=psi(j,i);
		save -ascii -append psi.in p
	endfor
endfor
