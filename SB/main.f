      program fft2d_schroedinger_nonliner
      implicit none 

      integer i,j,k,l,lmax
      integer nx,ny,hstep,nco,na

c      parameter (hstep=131072)
      parameter (hstep = 65536) ! 時間ステップ数：2の累乗で与えること！
      parameter (lmax = 3 )       ! 繰り返し数
                                  ! 計算の総時間は t*hstep*lmaxとなる

      parameter (nx = 1024)        ! x方向要素数：2の累乗で与えること！
      parameter (ny = 1024)        ! y方向要素数：2の累乗で与えること！

      real*8 pi
      complex*16 imu
      real*8 t,time,tw               ! t： 時間刻み幅，time: 時間

      ! 共振器パラメータ
      real*8 delta                  
      real*8 dscrad
      real*8 back

      real*8 nin
      real*8 nout
      real*8 neff

      !ブロッホパラメータ
      real*8 kpa
      real*8 gpp
      real*8 gpl
      real*8 winf
      real*8 delta0
      real*8 pstrg

      ! 座標パラメータ
      real*8 xmax,ymax
      real*8 x(0:nx-1),y(0:ny-1),kx(0:nx-1),ky(0:ny-1)

      ! 吸収パラメータ
      real*8 rad
      real*8 r2
      real*8 rho
      real*8 cw
      real*8 absorbs(0:nx-1,0:ny-1)

      !窓関数パラメータ
      real*8 sigma
      real*8 dist
      real*8 cut(0:nx-1,0:ny-1)
 
      !屈折率
      real*8 xde
      real*8 refs(0:nx-1,0:ny-1)
      !共振器形状
      integer conreg(0:nx-1,0:ny-1)

      real*8 asum,echsum

      !電場、分極と分布反転
      complex*16 echld(0:nx-1,0:ny-1)     ! 電場
      real*8 echld_t(0:nx-1,0:ny-1)
      complex*16 p(0:nx-1,0:ny-1)         ! 分極(双極子モーメント) 
      real*8 ww(0:nx-1,0:ny-1)            ! 分布反転
      real*8 w0(0:nx-1,0:ny-1)

      complex*16 ech
      complex*16 oldp(0:nx-1,0:ny-1)

      ! A作用素とB作用素
      complex*16 aa(0:nx-1,0:ny-1),bb(0:nx-1,0:ny-1)

      ! 光スペクトル　パラメータ
      real*8 tau(0:hstep-1),kt(0:hstep-1)
      complex*16 psref(0:hstep-1)          ! 光スペクトル用モニター

      ! その他
      real*8 ax,ay
      real*8 ra,rb,raa,rbb
      real*8 kata,Btop,Cdesu

      pi  = 4.00d0*datan(1.00d0)  ! pi
      imu = dcmplx(0.00d0,1.00d0) ! imaginary number  

      na   = nx/512               ! 共振器や計算空間の大きさ調整
      nco  = 0                    ! グラフ出力番号の初期値
     
c###########################################################
c##   基本パラメータ設定
c###########################################################

      t=1d-1     ! 時間刻み幅 
      time=0d0   ! 初期時刻

c------------------------------
c     共振器パラメータ
c------------------------------
      delta = 1d-10            ! delta = 1d-10 (スタジアム)
c      delta = 1.d0               ! delta = 1     (円のとき) 
      nin=2.0d0                    ! 共振器内の屈折率
      nout=1.0d0                   ! 共振器外の屈折率
      neff=(nout*nout)/(nin*nin)
      back= 4.0d-3                 ! 共振器内の電場：吸収係数

      dscrad = 49d0/4d0  ! 共振器の半径 

c------------------------------
c     Blochパラメータ
c------------------------------
      kpa=0.5d0                  ! 結合係数κ
      gpp=6.0d-2                 ! 分極の緩和率 γ
      gpl=3.0d-2                 ! 分布反転の緩和率 γ//
      winf=3.0d-3                ! ポンピングパワー
      delta0 = 0d0           ! 2準位原子の共鳴周波数

      pstrg=2d0*pi*kpa/(nin*nin) ! 電場・分極結合係数
      
c############################################################
c##      xy座標と波数などの設定
c############################################################
      xmax  = 20d0*pi*dble(na) ! 変更の必要なし
      ymax  = 20d0*pi*dble(na) ! 変更の必要なし
      
      call mkcord2dim(x,y,kx,ky,xmax,ymax,nx,ny)
      call tcord(tau,kt,hstep,t)
c###########################################################
c##   吸収境界と屈折率分布の設定
c###########################################################
c----------------------------------------
c　吸収境界
c----------------------------------------
! --- 吸収境界パラメータ   ------         
      r2 = 50d0*dble(na)  ! 変更の必要なし
      rho=0.6d0           ! 変更の必要なし
      cw=9d0              ! 変更の必要なし
!------------------------------
      call absorber(absorbs,x,y,r2,rho,cw,nx,ny)

c----------------------------------------
c　屈折率分布
c----------------------------------------
      ! 境界の滑らかさを表すパラメータ
      xde=dscrad*0.05d0/4d0   !変更の必要なし。

      call refset(refs,x,y,delta,dscrad,xde,neff,conreg,nx,ny)
      
c###########################################################
c##   窓関数の設定
c###########################################################

! --- 窓関数パラメータ   ------         
      dist  = 62.00d0*dble(na)  ! 変更の必要なし
      sigma = 3.00d0            ! 変更の必要なし
!------------------------------
      
      do j=0,ny-1
         do i=0,nx-1
            rad=dsqrt(x(i)*x(i)+y(j)*y(j))
            if (rad.lt.dist) then
               cut(i,j)=1d0-exp(-(rad-dist)*(rad-dist)/sigma)
            else
               cut(i,j)=0d0
            endif
         enddo
      enddo

!      do i=0,nx-1
!         write(11,*)x(i),cut(i,0)
!         write(12,*)x(i),refs(i,0)
!         write(13,*)x(i),absorbs(i,0)
!         do j=0,ny-1
!            if ( conreg(i,j) .eq. 1 ) then
!               write(14,*)x(i),y(j)
!            endif      
!         enddo
!      enddo
!      close(14)
      
c###########################################################
c##   初期条件の設定　（必要に応じて変更可)
c###########################################################
c--- 電場の初期状態 -----
      do i=0,nx-1
         do j=0,ny-1
            rad=(x(i)-5d0)*(x(i)-5d0)+(y(j)-1d0)*(y(j)-1d0)
            echld(i,j)=1d-4*dexp(-rad/25d0 )
            echld(i,j)=cut(i,j)*echld(i,j)
            echld_t(i,j)=0d0

            !call setmix(ech,x(i),y(j),dscrad)
            !echld(i,j) = ech
         enddo
      enddo
      call plot2d(echld,nco,na,nx,ny,ww,echld_t,tw)
      call sphere(echld,x,y,xmax,ymax,nx,ny,nco)
      !if ( nx.eq.ny) then
      !   stop
      !endif

c--- 分極pと分布反転ww の初期状態 -----
      do j=0,ny-1
         do i=0,nx-1
            p(i,j)=dcmplx(0d0,0d0)
            w0(i,j)=winf*conreg(i,j)
            ww(i,j)=winf*conreg(i,j)
         enddo
      enddo

c###########################################################
c##   時間発展演算子 A 作用素、B作用素の設定 (以下,変更必要なし）
c###########################################################      
      do j=0,ny-1
         do i=0,nx-1
            aa(i,j) = exp(imu*t*refs(i,j)/4d0
     &           -t*((back*dble(conreg(i,j))+absorbs(i,j))/2d0))
            bb(i,j) = exp(-(kx(i)*kx(i)+ky(j)*ky(j))
     &           *imu*t/2d0)
         enddo
      enddo

c###########################################################
c##   integration  ==sympletic integrater                 ##
c###########################################################

c     call load(echld,p,ww,time,nco,nx,ny)
      
      do l=1,lmax

         do k=1,hstep

            do j=0,ny-1
               do i=0,nx-1
                  oldp(i,j)=p(i,j)
               enddo
            enddo
*     ===========================
*     時間発展
*     nlop: Bloch eqs. 
*     prop: Propagation of Echld
*     ===========================
            call nlop(echld,p,ww,t,nx,ny,kpa,gpp,gpl,winf,delta0,w0)
            call prop(echld,aa,bb,nx,ny)
            do j=0,ny-1
               do i=0,nx-1
                  echld(i,j)=echld(i,j)+pstrg*oldp(i,j)*t*conreg(i,j)
                  echld(i,j)=cut(i,j)*echld(i,j)
                  echld_t(i,j)=echld_t(i,j)
     +                  +dble(abs(echld(i,j)*echld(i,j)))*t
               enddo
            enddo
            time=time+t
            tw=tw+t

            !光スペクトル用のモニター
            psref(k-1)=echld(57,110)
            
            if ( mod(k,10) .eq. 0 ) then
               echsum=0d0
               do j=0,ny-1
                  do i=0,nx-1
                     asum = dble(abs(echld(i,j)))
                     echsum=echsum
     +                    + asum*asum*conreg(i,j)
                  enddo
               enddo
!               write(0,*)time,echsum
               write(93,*)time,echsum
               
            endif

         enddo
         nco=nco+1
         call plot2d(echld,nco,na,nx,ny,ww,echld_t,tw) ! 電場強度分布
         call power(hstep,kt,psref,nco)  ! 光スペクトル
         call sphere(echld,x,y,xmax,ymax,nx,ny,nco)
         call memory(echld,p,ww,time,nco,nx,ny)

         tw=0d0
         do j=0,ny-1
            do i=0,nx-1
               echld_t(i,j) = 0d0
            enddo
         enddo

      enddo

      stop
      end



