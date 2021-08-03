	Lx=10;
	Ly=9;

	N_sites=106;
	N_up=53;
	N_dwn=53;

	my=(Ly+1)/2;

	dx=0.5;
	dy=sqrt(3.0)/2.0;

	sites=zeros(N_sites,2);
	k=1;
	for i=0:my-1
        	x=-i*dx;
        	y=-i*dy;
        	for j=1:Lx+i
                	sites(k,1)=x+j-1;
                	sites(k,2)=y;
                	k=k+1;
        	endfor
	endfor

	lmax=Lx+my-1;
	l=1;
	for i=my:Ly-1
        	x=x+dx;
        	y=-i*dy;
        	for j=1:lmax-l
                	sites(k,1)=x+j-1;
                	sites(k,2)=y;
                	k=k+1;
        	endfor
        	l=l+1;
	endfor


	H_0 = zeros(N_sites,N_sites);
	
	for i=1:N_sites
		for j=1:N_sites
			r2=(sites(i,1)-sites(j,1))**2+(sites(i,2)-sites(j,2))**2;
			if (r2 <= 1.5)
				if (r2 >= 0.5)
					H_0(i,j)=-1;
				endif
			endif
		endfor
	endfor
	
	[V,lambda] = eig(H_0);

	psi=horzcat(V(:,1:N_up),V(:,1:N_dwn));

	for i=1:N_up
        	for j=1:N_sites
                	p=psi(j,i);
                	save -ascii -append psi_hex.in p
        	endfor
	endfor

	n=diag(psi*psi');
	n

