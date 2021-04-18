module allocate_vars

use global

implicit none

contains

subroutine global_allocations

    implicit none

    integer :: alloc_err

    allocate(T    (Ny,Nx), stat=alloc_err)
    allocate(Tptrb(Ny,Nx), stat=alloc_err)
    allocate(uy   (Ny,Nx), stat=alloc_err)
    allocate(phi  (Ny,Nx), stat=alloc_err)
    allocate(ux   (Ny,Nx), stat=alloc_err)

    allocate(nlT(Ny,Nx), stat=alloc_err)
    allocate(nlphi(Ny,Nx), stat=alloc_err)

    allocate(phi1(Ny,Nx), phi2(Ny,Nx), stat=alloc_err)
    allocate(V1(Ny,Nx), V2(Ny,Nx), stat=alloc_err)

    allocate(dyv1_B(Nx), dyv1_T(Nx), stat=alloc_err)
    allocate(dyv2_B(Nx), dyv2_T(Nx), stat=alloc_err)

    allocate(tT (Nx), stat=alloc_err)
    allocate(tux(Nx), stat=alloc_err)
    allocate(tuy(Nx), stat=alloc_err)
    allocate(tnlT(Nx), stat=alloc_err)
    allocate(tnlphi(Nx), stat=alloc_err)
    allocate(tphi(Nx), stat=alloc_err)

    allocate(xp(Nx), stat=alloc_err)
    allocate(yp(Ny), stat=alloc_err)
    allocate(zp(Nz), stat=alloc_err)

    allocate(kx(Nx), stat=alloc_err)
    allocate(kz(Nz), stat=alloc_err)
    allocate(kx_modes(Nx), stat=alloc_err)

    allocate(dynu(Ny-1), stat=alloc_err)
    allocate(g1(Ny), g2(Ny), g3(Ny), stat=alloc_err)
    allocate(h1(Ny), h2(Ny), h3(Ny), stat=alloc_err)

    allocate(tmp_phi(Ny), tmp_T(Ny), tmp_uy(Ny), stat=alloc_err)
    allocate(tmp_phi1(Ny), tmp_uy1(Ny), stat=alloc_err)
    allocate(tmp_K_phi(Ny), tmp_K_T(Ny), stat=alloc_err)
    allocate(phii(Ny,Nx), Ti(Ny,Nx), stat=alloc_err)
    allocate(uyi(Ny,Nx), uxi(Ny,Nx), stat=alloc_err)
    allocate(K1_phi(Ny,Nx), K1_T(Ny,Nx), stat=alloc_err)
    allocate(K2_phi(Ny,Nx), K2_T(Ny,Nx), stat=alloc_err)
    allocate(K3_phi(Ny,Nx), K3_T(Ny,Nx), stat=alloc_err)
    allocate(K1hat_phi(Ny,Nx), K1hat_T(Ny,Nx), stat=alloc_err)
    allocate(K2hat_phi(Ny,Nx), K2hat_T(Ny,Nx), stat=alloc_err)
    allocate(K3hat_phi(Ny,Nx), K3hat_T(Ny,Nx), stat=alloc_err)
    allocate(K4hat_phi(Ny,Nx), K4hat_T(Ny,Nx), stat=alloc_err)

    if (alloc_err /= 0) then
       write(*,*) "ERROR:  Global allocations failed."
       stop
    end if

    T        = (0.0_dp, 0.0_dp)
    Tptrb    = (0.0_dp, 0.0_dp)
    ux       = (0.0_dp, 0.0_dp)
    phi      = (0.0_dp, 0.0_dp)
    uy       = (0.0_dp, 0.0_dp)
    nlT      = (0.0_dp, 0.0_dp)
    nlphi    = (0.0_dp, 0.0_dp)

    tT     = (0.0_dp, 0.0_dp)
    tux    = (0.0_dp, 0.0_dp)
    tuy    = (0.0_dp, 0.0_dp)
    tnlT   = (0.0_dp, 0.0_dp)
    tnlphi = (0.0_dp, 0.0_dp)
    tphi   = (0.0_dp, 0.0_dp)

    phi1 = 0.0_dp
    phi2 = 0.0_dp
    V1   = 0.0_dp
    V2   = 0.0_dp
    dyv1_T = 0.0_dp
    dyv2_T = 0.0_dp
    dyv1_B = 0.0_dp
    dyv2_B = 0.0_dp

    phii       = (0.0_dp, 0.0_dp)
    Ti         = (0.0_dp, 0.0_dp)
    uyi        = (0.0_dp, 0.0_dp)
    uxi        = (0.0_dp, 0.0_dp)
    K1_phi     = (0.0_dp, 0.0_dp)
    K2_phi     = (0.0_dp, 0.0_dp)
    K3_phi     = (0.0_dp, 0.0_dp)
    K1_T       = (0.0_dp, 0.0_dp)
    K2_T       = (0.0_dp, 0.0_dp)
    K3_T       = (0.0_dp, 0.0_dp)
    K1hat_phi  = (0.0_dp, 0.0_dp)
    K2hat_phi  = (0.0_dp, 0.0_dp)
    K3hat_phi  = (0.0_dp, 0.0_dp)
    K4hat_phi  = (0.0_dp, 0.0_dp)
    K1hat_T    = (0.0_dp, 0.0_dp)
    K2hat_T    = (0.0_dp, 0.0_dp)
    K3hat_T    = (0.0_dp, 0.0_dp)
    K4hat_T    = (0.0_dp, 0.0_dp)
    tmp_phi    = (0.0_dp, 0.0_dp)
    tmp_T      = (0.0_dp, 0.0_dp)
    tmp_uy     = (0.0_dp, 0.0_dp)
    tmp_K_phi  = (0.0_dp, 0.0_dp)
    tmp_K_T    = (0.0_dp, 0.0_dp)
    tmp_phi1   = (0.0_dp, 0.0_dp)
    tmp_uy1    = (0.0_dp, 0.0_dp)

    kx       = 0.0_dp
    kx_modes = 0.0_dp
    kz       = 0.0_dp

end subroutine global_allocations

subroutine global_deallocations

    implicit none

    deallocate(T  )
    deallocate(ux  )
    deallocate(uy  )
    deallocate(phi)

    deallocate(xp, yp, zp)
    deallocate(kx, kz)

end subroutine global_deallocations

end module allocate_vars
