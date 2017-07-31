SUBROUTINE svdcmp_sp(a,w,v)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
  IMPLICIT NONE
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
  REAL(SP), DIMENSION(:), INTENT(OUT) :: w
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
  INTEGER(I4B) :: i,its,j,k,l,m,n,nm
  REAL(SP) :: anorm,c,f,g,h,s,scale,x,y,z
  REAL(SP), DIMENSION(SIZE(a,1)) :: tempm
  REAL(SP), DIMENSION(SIZE(a,2)) :: rv1,tempn

  INTERFACE
     FUNCTION pythag(a,b)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: pythag
     END FUNCTION pythag
  END INTERFACE


	m=SIZE(a,1)
	n=assert_eq(SIZE(a,2),SIZE(v,1),SIZE(v,2),SIZE(w),'svdcmp_sp')
	g=0.0
	scale=0.0
	DO i=1,n
		l=i+1
		rv1(i)=scale*g
		g=0.0
		scale=0.0
		IF (i <= m) THEN
			scale=SUM(ABS(a(i:m,i)))
			IF (scale /= 0.0) THEN
				a(i:m,i)=a(i:m,i)/scale
				s=DOT_PRODUCT(a(i:m,i),a(i:m,i))
				f=a(i,i)
				g=-SIGN(SQRT(s),f)
				h=f*g-s
				a(i,i)=f-g
				tempn(l:n)=MATMUL(a(i:m,i),a(i:m,l:n))/h
				a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
				a(i:m,i)=scale*a(i:m,i)
			END IF
		END IF
		w(i)=scale*g
		g=0.0
		scale=0.0
		IF ((i <= m) .AND. (i /= n)) THEN
			scale=SUM(ABS(a(i,l:n)))
			IF (scale /= 0.0) THEN
				a(i,l:n)=a(i,l:n)/scale
				s=DOT_PRODUCT(a(i,l:n),a(i,l:n))
				f=a(i,l)
				g=-SIGN(SQRT(s),f)
				h=f*g-s
				a(i,l)=f-g
				rv1(l:n)=a(i,l:n)/h
				tempm(l:m)=MATMUL(a(l:m,l:n),a(i,l:n))
				a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
				a(i,l:n)=scale*a(i,l:n)
			END IF
		END IF
	END DO
	anorm=MAXVAL(ABS(w)+ABS(rv1))
	DO i=n,1,-1
		IF (i < n) THEN
			IF (g /= 0.0) THEN
				v(l:n,i)=(a(i,l:n)/a(i,l))/g
				tempn(l:n)=MATMUL(a(i,l:n),v(l:n,l:n))
				v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
			END IF
			v(i,l:n)=0.0
			v(l:n,i)=0.0
		END IF
		v(i,i)=1.0
		g=rv1(i)
		l=i
	END DO
	DO i=MIN(m,n),1,-1
		l=i+1
		g=w(i)
		a(i,l:n)=0.0
		IF (g /= 0.0) THEN
			g=1.0_sp/g
			tempn(l:n)=(MATMUL(a(l:m,i),a(l:m,l:n))/a(i,i))*g
			a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
			a(i:m,i)=a(i:m,i)*g
		ELSE
			a(i:m,i)=0.0
		END IF
		a(i,i)=a(i,i)+1.0_sp
	END DO
	DO k=n,1,-1
		DO its=1,30
			DO l=k,1,-1
				nm=l-1
				IF ((ABS(rv1(l))+anorm) == anorm) EXIT
				IF ((ABS(w(nm))+anorm) == anorm) THEN
					c=0.0
					s=1.0
					DO i=l,k
						f=s*rv1(i)
						rv1(i)=c*rv1(i)
						IF ((ABS(f)+anorm) == anorm) EXIT
						g=w(i)
						h=pythag(f,g)
						w(i)=h
						h=1.0_sp/h
						c= (g*h)
						s=-(f*h)
						tempm(1:m)=a(1:m,nm)
						a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
						a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
					END DO
					EXIT
				END IF
			END DO
			z=w(k)
			IF (l == k) THEN
				IF (z < 0.0) THEN
					w(k)=-z
					v(1:n,k)=-v(1:n,k)
				END IF
				EXIT
			END IF
			IF (its == 30) CALL nrerror('svdcmp_sp: no convergence in svdcmp')
			x=w(l)
			nm=k-1
			y=w(nm)
			g=rv1(nm)
			h=rv1(k)
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
			g=pythag(f,1.0_sp)
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
			c=1.0
			s=1.0
			DO j=l,nm
				i=j+1
				g=rv1(i)
				y=w(i)
				h=s*g
				g=c*g
				z=pythag(f,h)
				rv1(j)=z
				c=f/z
				s=h/z
				f= (x*c)+(g*s)
				g=-(x*s)+(g*c)
				h=y*s
				y=y*c
				tempn(1:n)=v(1:n,j)
				v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
				v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
				z=pythag(f,h)
				w(j)=z
				IF (z /= 0.0) THEN
					z=1.0_sp/z
					c=f*z
					s=h*z
				END IF
				f= (c*g)+(s*y)
				x=-(s*g)+(c*y)
				tempm(1:m)=a(1:m,j)
				a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
				a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
			END DO
			rv1(l)=0.0
			rv1(k)=f
			w(k)=x
		END DO
	END DO
	END SUBROUTINE svdcmp_sp

	SUBROUTINE svbksb_sp(u,w,v,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(SP), DIMENSION(SIZE(x)) :: tmp
	mdum=assert_eq(SIZE(u,1),SIZE(b),'svbksb_sp: mdum')
	ndum=assert_eq((/SIZE(u,2),SIZE(v,1),SIZE(v,2),SIZE(w),SIZE(x)/),&
		'svbksb_sp: ndum')
	WHERE (w /= 0.0)
		tmp=MATMUL(b,u)/w
	ELSEWHERE
		tmp=0.0
	END WHERE
	x=MATMUL(v,tmp)
	END SUBROUTINE svbksb_sp

	FUNCTION pythag(a,b)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: pythag
	REAL(SP) :: absa,absb
	absa=ABS(a)
	absb=ABS(b)
	IF (absa > absb) THEN
		pythag=absa*SQRT(1.0_sp+(absb/absa)**2)
	ELSE
		IF (absb == 0.0) THEN
			pythag=0.0
		ELSE
			pythag=absb*SQRT(1.0_sp+(absa/absb)**2)
		END IF
	END IF
	END FUNCTION pythag
