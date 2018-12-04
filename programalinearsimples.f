
      program main
c
c	programa com elementos lineares simples 
c


	implicit none
      common ne,np,l,lec,imp,npi,nt
	integer nx,lec,imp
	real*8 hdif
	parameter (nx=500, hdif = 1.0E-6)
	integer nop(nx,3),kode(nx),locc(nx),idup(nx),ne,np,l,npi,nt
	real*8 x(nx),y(nx),g(nx,nx),fi(nx),dfi(nx),cx(nx)
	real*8 cy(nx),sol(nx),h(nx,nx),dsolx(nx),dsoly(nx),ccy(nx)
	real*8 a(nx,nx),b(nx,nx),cb(nx),xx(nx),yy(nx),d,ccx(nx)
	real*8 tempo,calcula_tempo
	character sai*20,tempo1*10,tempo2*10
	lec=5
	imp=6

      call input(cx,cy,x,y,kode,fi,nop,locc,Idup,nx,sai,ccx,ccy)
c	write(imp,*)'vou entrar na fmat'
c
c		CALCULO DO TEMPO DE PROCESSAMENTO
c
      call date_and_time(TIME=tempo1)

      call fmat(x,y,g,h,fi,dfi,kode,nop,Idup,a,b,cb,xx,yy,nx,ccx,ccy,
	*hdif)
c	write(imp,*)'ja sai da fmat'
      nt=np+npi 
      call slnpd(h,dfi,d,np,nx)
c	write(imp,*)'ja sai da slnpd'
      call inter(fi,dfi,kode,nx)
c	write(imp,*)'ja sai da inter'
c
c		CALCULO DO TEMPO DE PROCESSAMENTO
c
c      call date_and_time(TIME=tempo2)      
c      tempo=calcula_tempo(tempo1,tempo2)

	call pontos_internos(g,h,nx,fi,dfi,np,npi,hdif)

      call output (x,y,fi,dfi,cx,cy,sol,dsolx,dsoly,nx,sai,tempo)
      stop
      end

	real*8 function calcula_tempo(tempo1,tempo2)
	implicit none
	real*8 temp1(4),temp2(4)
	character tempo1*10,tempo2*10,aux*2,aux2*3

	aux=tempo1(1:2)
	read(aux,*) temp1(1)
	aux=tempo1(3:4)
	read(aux,*) temp1(2)
	aux=tempo1(5:6)
	read(aux,*) temp1(3)
	aux2=tempo1(8:10)
	read(aux2,*) temp1(4)
	aux=tempo2(1:2)
	read(aux,*) temp2(1)
	aux=tempo2(3:4)
	read(aux,*) temp2(2)
	aux=tempo2(5:6)
	read(aux,*) temp2(3)
	aux2=tempo2(8:10)
	read(aux2,*) temp2(4)
	if (temp1(1).gt.temp2(1)) then
	  temp2(1)=temp2(1)+24
	end if
	temp1(1)=temp1(1)*3600+temp1(2)*60+temp1(3)+temp1(4)/1000
	temp2(1)=temp2(1)*3600+temp2(2)*60+temp2(3)+temp2(4)/1000
	calcula_tempo=temp2(1)-temp1(1)
	end

      subroutine input(cx,cy,x,y,kode,fi,nop,locc,Idup,nx,sai,ccx,ccy)

	implicit none
	common ne,np,l,lec,imp,npi,nt
	integer nx,lec,imp,ne,np,npi,l,nt,i,locc(nx),kode(nx)
	integer idup(nx),nop(nx,3)
	real*8 cx(nx),cy(nx),x(nx),y(nx),fi(nx),ccx(nx),ccy(nx)
	character title*18,arqent*20,arqout*20,sai*20


c	open(unit=1,file='entrada.txt',status='old')
c	read(1,7) ARQENT
c	read(1,7) ARQOUT
c	read(1,7) sai
c	close(unit=1,status='delete')
c7	format (a20)

	write(*,'(A\)')' INFORME NOME DO ARQUIVO DE ENTRADA-->'
      read(*,'(A20)')ARQENT
      open(LEC,FILE=ARQENT)
      write(*,'(A\)')' INFORME NOME DO ARQUIVO DE SAIDA --->'
      read(*,'(A20)')ARQOUT
      open(IMP,FILE=ARQOUT)
c

      write(imp,100)
100   format(' ',120('*'))

      read(lec,150) title
150   format(18a4)
      write(imp,250) title
250   format(25x,18a4)

	read(lec,200) ne
	read(lec,200) np
      read(lec,200) npi
200   format(i5)
      nt=np+npi
      write(imp,300) ne,np,npi
300   format(//'data'//2x,'number of boundary elements=',
     1i3//,2x,'number of functional nodes=',i3//,2x,
     2'number of internal points where the function is calculated='
     3,i3,//)

      write(imp,500)
500   format (//2x,'coordinates of the extreme points of the boundary el
     1ements',//4x,'point',10x,'x',18x,'y')
	nt=np+npi
      do 750 i=1,nt
      read(lec,600)i, x(i),y(i)
600   format(i3,4f12.4)
      write(imp,700) i,x(i),y(i)
700   format(5x,i3,4(5x,e14.7))
750	continue

      write(imp,800)
800   format(//2x,'boundary conditions'//5x,'node',6x,'code',
     15x,'prescribed value')
      do 20 i=1,np
      read(lec,900) kode(i),Idup(i),fi(i)
900   format(I5,i6,f12.4)
      write(imp,950) i,kode(i),Idup(i),fi(i)
950   format (5x,i3,8x,I1,8x,I3,8x,e14.7)
 20   continue


	write(imp,880)
880   format(//,2x,'incidence of geometric nodes'//,5x,' first node',
     16x,'second node',/)
      do 261 i=1,ne
      read(lec,920) nop(i,1),nop(i,2)
920   format(i5,i5)
    	write(imp,960) i,nop(i,1),nop(i,2)
960   format(5x,i3,8x,i3,5x,i3)
261   continue

      return
      end
c
c
	subroutine slnpd(a,b,d,n,nx)
c      
      implicit none
	integer n1,n,k,k1,j,l,nx,i
	real*8 c,a(nx,nx),b(nx),d

      n1=n-1
      do 100 k=1,n1
      k1=k+1
      c=a(k,k)
      if (dabs(c)-0.000001)1,1,3
  1   do 7 j=k1,n
      if (dabs(a(j,k))-0.000001)7,7,5
  5   do 6 l=k,n
      c=a(k,l)
      a(k,l)=a(j,l)
  6   a(j,l)=c
      c=b(k)
      b(k)=b(j)
      b(j)=c
      c=a(k,k)
      go to 3
  7   continue
  8   write(6,2) k
  2   format('**** singularity in row',i5)
      d=0.
      go to 300
      
  3   c=a(k,k)
      do 4 j=k1,n
  4   a(k,j)=a(k,j)/c
      b(k)=b(k)/c
      do 10 i=k1,n
      c=a(i,k)
      do 9 j=k1,n
  9   a(i,j)=a(i,j)-c*a(k,j)
 10   b(i)=b(i)-c*b(k)
100   continue

      if (dabs(a(n,n))-0.000001)8,8,101
  
101   b(n)=b(n)/a(n,n)

      do 200 l=1,n1
      k=n-l
      k1=k+1
      do 200 j=k1,n
200   b(k)=b(k)-a(k,j)*b(j)

      d=1.
      do 250 i=1,n
250   d=d*a(i,i)
300   return
      end
c
c    

      subroutine fmat(x,y,g,h,fi,dfi,kode,nop,Idup,a,b,cb,xx,yy,nx,
     *ccx,ccy,hdif)
      implicit none
	common ne,np,l,lec,imp,npi,nt
	integer nx,nop(nx,3),idup(nx),l1,l3,i,ne,l2,l4,np,imp
	integer npi,nt,j,ng,l,lec,kode(nx),con1,con2,k,n,j11,j12,j21,j22
	integer count
	real*8 a11,a12,a21,a22,hm(nx,nx),v(2*nx,2*nx)
	real*8 xd(nx),x(nx),yd(nx),y(nx),xx(nx),yy(nx)
	real*8 g(nx,nx),h(nx,nx),h1,h2,g1,g2,dx1,dx2,ccx(nx),ccy(nx)
	real*8 dy1,dy2,ex1,ex2,ey1,ey2,ta,ra,c(nx,nx),fdo(nx),ge
	real*8 p(nx,nx),cd(nx),d,r,tok(nx,nx),alfa(nx,nx),teta
	real*8 gg,oneta(nx),ch,dfi(nx),cb(nx),fi(nx),q(nx,nx)
	real*8 a(nx,nx),b(nx,nx),cc(nx,nx),dc(nx),ag,k1,k2,pi
	real*8 dist,ala,am(nx,nx),f(nx,nx),a1,c1,a2,c2,FF(NX,NX)
	REAL*8 AX,AY,COMP,SENO(NX),COSE(NX),cci(nx),GA(NX,NX)
	parameter (pi=3.141592)
	real*8 hdif, xaux(nx),yaux(nx)

c	organização dos vetores xaux e yaux, que são os vetores com acréscimo de hdif

	do 199 i=1,npi
	xaux(i) = x(np+i) + hdif
	yaux(i) = y(np+i)
199	continue
	do 299 i=(npi+1),2*npi
	xaux(i) = x(np+i-npi)
	yaux(i) = y(np+i-npi) + hdif
299	continue
	do 399 i=1,2*npi
	x(nt+i) = xaux(i)
	y(nt+i) = yaux(i)
399	continue


c	 afastamento das coordenadas no caso de condições de Dirichlet nos nós duplos!
c      esta alteração vai afetar as matrizes H e G! As coordenadas x são alteradas
c      apenas no caso de dirichlet nos dois nós que são duplos!

	dist=0.13 !porcentagem do afastamento do elemento da aresta
 	kode(0)=1
	do 654 i=1,ne
	l1=nop(i,1)
	l3=idup(l1)
	if (l3.eq.0) go to 321
	ala=kode(l1)+kode(l3)
 	if(ala.ne.0) go to 321
	l2=nop(i,2)
	x(l1)=x(l1)+dist*(x(l2)-x(l1))
	y(l1)=y(l1)+dist*(y(l2)-y(l1))
321     continue
      l2=nop(i,2)
	l4=idup(l2)
	if (l4.eq.0) go to 543
	ala=kode(l2)+kode(l4)
	if (ala.ne.0) go to 543
	x(l2)=x(l2)-dist*(x(l2)-x(l1))
	y(l2)=y(l2)-dist*(y(l2)-y(l1))
543     continue
654     continue


c	 Usaremos o número total de pontos nodais, sendo a soma de pontos do contorno com os pontos
c	 internos, para calcular as matrizes H e G.
c
c	geração das matrizes h e g

	count = nt+2*npi

      do 10 i=1,count
      do 10 j=1,count
      g(i,j)=0
 10   h(i,j)=0
  	DO 12 I=1,count
	DO 18 J=1,NE
	L1=NOP(J,1)
	L2=NOP(J,2)
  	IF(I.EQ.L1)GO TO 20
	IF(I.EQ.L2)GO TO 22
      call inte(x(i),y(i),x(L1),y(L1),x(L2),y(L2),h1,h2,g1,g2,dx1,
     *dx2,dy1,dy2,ex1,ex2,ey1,ey2) 
      h(i,L2)=h(i,L2)+h2
      h(i,L1)=h(i,L1)+h1
      g(i,L2)=g(i,L2)+g2
      g(i,L1)=g(i,L1)+g1
	GO TO 26
  20  continue
      call inlo(x(L1),y(L1),x(L2),y(L2),g1,g2)
      g(i,L1)=g(i,L1)+g1
	h(i,L1)=h(i,L1)+0.0
	g(i,L2)=g(i,L2)+g2
	h(i,L2)=h(i,L2)+0.0
	go to 26
  22  continue
      call inlo(x(L1),y(L1),x(L2),y(L2),g1,g2)
	g(i,L1)=g(i,L1)+g2
 	h(i,L1)=h(i,L1)+0.0
	g(i,L2)=g(i,L2)+g1
	h(i,L2)=h(i,L2)+0.0
  26  continue
 18   continue
 12   continue
c  
c	 construção dos termos da diagonal de h
c 
 	DO 29 I=1,Nt
  	ta=0.0
	DO 28 J=1,Nt
	ta=ta+h(i,j) 
 28   continue
      h(i,i)=-ta
 29   continue

c	****************************************************************************
c	Neste seção serão realizadas as impressões das matrizes no arquivo de saída
c
c	impressão de matrizes
c
c	write(imp,*)'matriz g'
c	do 542 i=1,nt
c	write(imp,541)(g(i,j),j=1,nt)
c541   format(10f9.4)
c542   continue
c
c	write(imp,*)'matriz h'
c      do 545 i=1,nt
c	write(imp,546)(h(i,j),j=1,nt)
c546   format(10f9.4)
c545   continue
c
c	Fim da impressão das matrizes no arquivo de saída
c	**************************************************

      do 152 j=1,np
      if (kode(j)) 140,140,150
140   continue
      do 151 i=1,np
      ch=g(i,j)
      g(i,j)=-h(i,j)
      h(i,j)=-ch
151   continue
150   continue
152   continue

      do 161 i=1,np
      dfi(i)=0.
      do 160 j=1,np
      dfi(i)=dfi(i)+g(i,j)*fi(j)
160   continue
161   continue

c
c      write(imp,*)'montagem de fdi'
c	do 96 i=1,np
c   	write(imp,*)dfi(i) 
c96    continue

      return
      end

      subroutine inte(xp,yp,x1,y1,x2,y2,h1,h2,g1,g2,dx1,dx2,dy1,dy2,
 	*ex1,ex2,ey1,ey2)  
      implicit none  
	common ne,np,l,lec,imp,npi,nt
	integer i,ne,np,l,lec,imp,npi,nt
	real*8 gi(20),ome(20),ax,x1,x2,bx,ay,y1,y2,by,comp,h1,h2
	real*8 g1,g2,dx1,dx2,dy1,dy2,ex1,ex2,ey1,ey2,pi2,xco(20)
	real*8 yco(20),ra,xp,yp,aux,h,g,cjac,gx,gy,hx,hy,sen,cos


	 GI(1)=0.993128599185094
       GI(2)=-GI(1)
       GI(3)=0.963971927277913
       GI(4)=-GI(3)
       GI(5)=0.912234428251325
       GI(6)=-GI(5)
       GI(7)=0.839116971822218
       GI(8)=-GI(7)
       GI(9)=0.746331906460150
       GI(10)=-GI(9)
       GI(11)=0.636053680726515
       GI(12)=-GI(11)
       GI(13)=0.510867001950827
       GI(14)=-GI(13)
       GI(15)=0.373706088715419
       GI(16)=-GI(15)
       GI(17)=0.227785851141645
       GI(18)=-GI(17)
       GI(19)=0.076526521133497
       GI(20)=-GI(19)
*
*      pesos de Gauss
*
       OME(1)=0.017614007139152
       OME(2)=OME(1)
       OME(3)=0.040601429800386
       OME(4)=OME(3)
       OME(5)=0.062672048334109
       OME(6)=OME(5)
       OME(7)=0.083276741576704
       OME(8)=OME(7)
       OME(9)=0.101930119817240
       OME(10)=OME(9)
       OME(11)=0.118194531961518
       OME(12)=OME(11)
       OME(13)=0.131688638449176
       OME(14)=OME(13)
       OME(15)=0.142096109318382
       OME(16)=OME(15)
       OME(17)=0.149172986472603
       OME(18)=OME(17)
       OME(19)=0.152753387130725
       OME(20)=OME(19)


      ax=(x2-x1)/2
      bx=(x2+x1)/2
      ay=(y2-y1)/2
      by=(y1+y2)/2
	comp=2*dsqrt(ax**2+ay**2)
	cjac=comp/2
      SEN=(X1-X2)/COMP
      COS=(Y2-Y1)/COMP
		     
      h1=0
      h2=0
      g1=0
      g2=0
	dx1=0
	dx2=0
	dy1=0
	dy2=0
	ex1=0
	ex2=0
	ey1=0
	ey2=0
	pi2=6.28318
      do 40 i=1,20
      xco(i)=ax*gi(i)+bx
      yco(i)=ay*gi(i)+by
      ra=dsqrt((xco(i)-xp)**2+(yco(i)-yp)**2)
	aux=(xco(i)-xp)*cos+(yco(i)-yp)*sen
      h=-aux*ome(i)*cjac/(pi2*ra**2)
      g=dlog(1/ra)*ome(i)*cjac/pi2
      gx=ome(I)*cjac*(XCO(i)-XP)/(pi2*ra**2)
	gy=ome(I)*cjac*(YCO(i)-YP)/(pi2*ra**2)
	hx=OME(I)*cjac*(2*AUX*(XP-XCO(i))+COS*ra*ra)/(pi2*ra**4)
	hy=OME(I)*cjac*(2*AUX*(YP-YCO(i))+SEN*ra*ra)/(pi2*ra**4)
	h1=h1-(gi(i)-1)*h/2
      h2=h2+(gi(i)+1)*h/2
      g1=g1-(gi(i)-1)*g/2
      g2=g2+(gi(i)+1)*g/2	
 	dx1=dx1-(gi(i)-1)*gx/2
      dx2=dx2+(gi(i)+1)*gx/2
      dy1=dy1-(gi(i)-1)*gy/2
      dy2=dy2+(gi(i)+1)*gy/2
      ex1=ex1-(gi(i)-1)*hx/2
      ex2=ex2+(gi(i)+1)*hx/2
      ey1=ey1-(gi(i)-1)*hy/2
 40   ey2=ey2+(gi(i)+1)*hy/2
      return
      end

      subroutine inlo(x1,y1,x2,y2,g1,g2)
	implicit none
	real*8 sep,g1,x1,x2,y1,y2,g2

      sep=dsqrt((x2-x1)**2+(y2-y1)**2)
	g1=sep*(1.5-dlog(sep))/(2*6.28318)
      g2=sep*(0.5-dlog(sep))/(2*6.28318)
      continue
      return
      end
c
c
c
      subroutine inter(fi,dfi,kode,nx)
	implicit none
	common ne,np,l,lec,imp,npi,nt
	integer nx,kode(nx),i,np,k,l,m,locc(nx),l1,l2,nop(nx,3)
	integer j,ne,nt,lec,imp,npi
	real*8 ch,fi(nx),dfi(nx),sol(nx),dsolx(nx),dsoly(nx)
	real*8 cx(nx),cy(nx),x(nx),y(nx),a1,a2,b1,b2,dx1,dx2
	real*8 dy1,dy2,ex1,ex2,ey1,ey2,drumx,drumy,zog1,zog2
	real*8 hcsi(nx,nx),xx(nx),yy(nx),rex(nx,nx)
	real*8 doma(nx,nx),domb(nx,nx),gcsi(nx,nx),c1
	real*8 a(nx,nx),b(nx,nx),domt(nx),cb(nx)


      do 20 i=1,np
c	write(imp,*)'kode',kode(i)
       if (kode(i)) 20,20,10
 10   ch=fi(i)
      fi(i)=dfi(i)
      dfi(i)=ch
 20   continue

      return
      end
c
c
c
	 SUBROUTINE MINV(A,N)
	implicit none
	integer imp,k,n,j,i,ji,ki
	real*8 lo(500),mo(500),biga,a(500,500),hold

C     
C      INVERTE A MATRIZ ATRAVES DO METODO PADRAO DE GAUSS-JORDAN
C      O DETERMINANTE TAMBEM E CALCULADO. UM DETERMINANTE NULO INDICA QUE A 
C      MATRIZ E SINGULAR.
C
C      A - MATRIZ DE ENTRADA, DESTRUIDA NA OPERACAO E SUBSTITUIDA PELA
C          INVERSA RESULTANTE
C      N - ORDEM DA MATRIZ A
C      LO - VETOR DE TRABALHO
C      MO - VETOR DE TRABALHO
c     
      IMP=6
C
C
C      PROCURA DO MAIOR ELEMENTO
C
       DO 100 K=1,N
         LO(K)=K
         MO(K)=K
         BIGA=A(K,K)
         DO 10 J=K,N
         DO 10 I=K,N
           IF(dabs(BIGA).GE.dabs(A(I,J))) GO TO 10
           BIGA=A(I,J)
           LO(K)=I
           MO(K)=J
10       CONTINUE
C
C      MUDANCA DE LINHAS
C
         J=LO(K)
         IF(J.LE.K) GO TO 30
         DO 20 I=1,N
           HOLD=-A(K,I)
c           JI=KI-K+J
           A(K,I)=A(J,I)
           A(J,I)=HOLD
20       CONTINUE
C
C      MUDANCA DE COLUNAS
C
30       I=MO(K)
         IF(I.LE.K) GO TO 50
         DO 40 J=1,N
           HOLD=-A(J,K)
           A(J,K)=A(J,I)
           A(J,I)=HOLD
40       CONTINUE
C
C      DIVIDE COLUNA POR PIVOT NEGATIVO ( VALOR DOS ELEMENTOS DO PIVOT ESTAO
C      CONTABILIZADOS EM BIGA)
C
50       IF((dabs(BIGA)).GT.0)  GO TO 60
           WRITE (imp,160) K
           STOP
60       DO 70 I=1,N
           IF(I.EQ.K) GO TO 70
           A(I,K)=A(I,K)/(-BIGA)
70       CONTINUE
C
C      REDUCAO DA MATRIZ
C
         DO 80 I=1,N
           HOLD=A(I,K)
           DO 80 J=1,N
             IF(I.EQ.K) GO TO 80
             IF(J.EQ.K) GO TO 80
             A(I,J)=HOLD*A(K,J)+A(I,J)
80       CONTINUE
C
C      DIVIDE A LINHA PELO PIVOT
C
         DO 90 J=1,N
           IF(J.EQ.K) GO TO 90
           A(K,J)=A(K,J)/BIGA
90       CONTINUE
C
C      SUBSTITUI O PIVOT PELO RECIPROCO
C
         A(K,K)=1.0/BIGA
100    CONTINUE
C
C      MUDANCA FINAL DE LINHAS E COLUNAS
C
       K=N-1
110    CONTINUE
       I=LO(K)
       IF(I.LE.K) GO TO 130
       DO 120 J=1,N
         HOLD=A(J,K)
         A(J,K)=-A(J,I)
         A(J,I)=HOLD
120    CONTINUE
130    J=MO(K)
       IF(J.LE.K) GO TO 150
       DO 140 I=1,N
         HOLD=A(K,I)
         A(K,I)=-A(J,I)
         A(J,I)=HOLD
140    CONTINUE
150    K=K-1
       IF(K.GT.0)GO TO 110
       RETURN
160    FORMAT(/,20X,'DURANTE A OPERACAO DE INVERSAO DA MATRIZ',/,20X,
     1'DETERMINANTE NULO FOI DETECTADO.  O PROCESSAMENTO FOI',/,20X,
     1'INTERROMPIDO  O NUMERO DA LINHA OU COLUNA FOI :',I4,/,20x,
     1'dentro da rotina MINV - inversao de matriz')
       END

      subroutine output(x,y,fi,dfi,cx,cy,sol,dsolx,dsoly,nx,sai,tempo)
	implicit none
      common ne,np,l,lec,imp,npi,nt
	common dx1,dx2,dy1,dy2,ex1,ex2,ey1,ey2
 	integer	i,nt,k,l,ne,np,lec,imp,npi,nx,ini,fim
	real*8 x(nx),y(nx),fi(nx),dfi(nx),cx(nx),cy(nx)
	real*8 sol(nx),dsolx(nx),dsoly(nx),ti(nx)
	real*8 total,sn,cs,ym(nx),xm(nx),dx1,dx2
	real*8 dy1,dy2,ex1,ex2,ey1,ey2,aux,total2,tempo
	character sai*20

	write(imp,77) tempo
 77	format(/'Tempo de Processamento: ',f8.3,' segundos'/)
      write(imp,100)
100   format (' ',120('*'),//1x,'results'//2x,'boundary nodes'//16x,
     1'x',23x,'y',19x,'potential',10x,'potential derivative'/)
      do 10 i=1,nt
 10   write(imp,200) x(i),y(i),fi(i),dfi(i)
200   format (4(10x,e14.7))

c
c	     APENDICE PARA CALCULO DO ERRO COMETIDO


      ini=2
	if (ne.eq.32) then
	  ini=28
        fim=33
	else if (ne.eq.80) then
	  ini=64
	  fim=75
	else if (ne.eq.120) then
	  ini=94
	  fim=110
	else
	  ini=124
	  fim=145
	end if
	aux=0
	total=0.0
      DO 611 i=ini,fim


c	     RESPOSTA ANALITICA CONSTANTE
c	  ti=x(i)*(1-x(i)/2)
c	  ti=(1.0-x(i))
c	      RESPOSTA ANALITICA LINEAR
c	  ti=(x(i)/2)*(1-((x(i)**2)/3))
c	  ti=0.5-x(i)*2/3
c	     RESPOSTA ANALITICA SENOIDAL
c	  ti(i)=dsin(x(i))-x(i)*cos(1.0)
c       ti(i)=dcos(x(i))-cos(1.0)
c		      CALCULO DO ERRO
c	write(*,*)i,x(i),ti(i)
c	total=total+abs((abs(ti(i))-abs(dfi(i)))/ti(i))*100
c	aux=aux+1
	
c	WRITE(imp,581)dfi(i),ti(i)
c581   format (2X,f10.4,2x,f10.4)
611   CONTINUE

c      total=total/aux
c	write(imp,*)total


      return
      end


	!Subrotina para calcular o valor dos pontos internos

	subroutine pontos_internos(g,h,nx,fi,dfi,np,npi,hdif)
	
	implicit none

	integer np,npi,nx
	integer i,j,ini,fim
	
      real*8 g(nx,nx),h(nx,nx),fi(nx),dfi(nx)
	real*8 dfix(nx),dfiy(nx)
	real*8 somafi,somadfi
	real*8 hdif
	integer count

	ini = np+1
	fim = np+3*npi

c	Cálculo dos pontos internos

	do 980 i=ini,fim
	somafi  = 0.0
	somadfi = 0.0
	do 990 j=1,np
	somafi  = somafi + h(i,j) * fi(j)
	somadfi = somadfi + g(i,j) * dfi(j)
990	continue
	fi(i)  = somadfi - somafi
980	continue

c	Cálculo da derivada para os pontos internos

	count = 1

	do 970 i=(np+npi+1),(np+npi*2)
	dfix(count) = (fi(i)     - fi(i-npi))/hdif
	dfiy(count) = (fi(i+npi) - fi(i-npi))/hdif
	count = count + 1
970   continue 		

	count = 1

	do 123 i=ini,fim
	dfi(i) = sqrt((dfix(count)*dfix(count))+(dfiy(count)*dfiy(count))) 
	count = count + 1
123	continue

	return
	end