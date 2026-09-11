module imyann_nuppn_c_api
   use iso_c_binding, only: c_char, c_double, c_int, c_null_char
   use array_sizes, only: i282dim, i325dim, iAtdim, iCfdim, nre, nsp
   use alpha_decays, only: alpha_decays_init
   use bader_deuflhard, only: bader_deuflhard_init
   use constants, only: HALF, ONE, ZERO
   use evaluate_rates, only: evaluate_all_rates
   use frame_knobs, only: istart, readframeinput
   use fuller, only: fuller_init
   use jac_rhs, only: calculate_dxdt, dxdt, makemap
   use jbj16, only: jbj_init
   use kadonis, only: kadonis_init
   use netgen, only: netgen_init
   use nkk04, only: nkk_init
   use nse_swj, only: compute_nse, nse_init
   use nuc_data, only: considerisotope, considerreaction, ispe, niso, nuc_data_init, zis
   use other_nuc, only: other_nuc_init
   use physics, only: rnetw2007
   use physics_knobs, only: detailed_balance, index_reaclib, ininet, &
      readphysicsinput, rho_nw_ini, t9_nw_ini, use_cache, ye_nw_ini, yps_nw_ini
   use rates, only: k1, k2, k3, k4, k5, k6, k7, k8, rates_init, v
   use reaclib, only: reaclib_create_masks, reaclib_init, reaclib_preprocessor
   use reaction_info, only: bind_energy_diff, i_ba, i_bm, i_bn, i_bp, i_ec, &
      i_ve, i_vp, i_vsa, i_vsn, i_vsp, ilabb, reaction_info_init
   use reverse, only: reverse_init
   use screening, only: screen_init
   use solver, only: integrate_network, solver_init
   use solver_diagnostics, only: nsubt
   use solver_knobs, only: readsolverinput
   use utils, only: calculate_ye, r8
   use vital, only: vital_init
   use neutrinos, only: neutrino, neutrinos_init
   implicit none
   private

   public :: nuppn_init, nuppn_num_species, nuppn_get_species, &
      nuppn_integrate, nuppn_integrate_to_time, nuppn_last_substeps, &
      nuppn_history_size, nuppn_get_history_step, &
      nuppn_compute_nse, nuppn_get_dxdt, nuppn_evaluate_reactions, &
      nuppn_num_reactions, nuppn_get_reaction

   logical :: initialized = .false.
   logical :: nse_ready = .false.
   integer(c_int) :: nactive = 0
   integer(c_int) :: nactive_reactions = 0
   integer(c_int) :: nvar_saved = 0, nvrel_saved = 0
   integer(c_int) :: active_species(nsp) = 0
   integer(c_int) :: active_reactions(nre) = 0
   integer(c_int) :: ppn_to_active(nsp) = -1
   integer(c_int) :: last_substeps = 0
   integer(c_int) :: last_history_count = 0
   integer(c_int) :: last_history_capacity = 0
   real(r8) :: reaction_flow(nre) = ZERO
   real(r8) :: an(nsp), zn(nsp), qi(nre)
   real(c_double), allocatable :: last_history_time(:)
   real(c_double), allocatable :: last_history_rho(:)
   real(c_double), allocatable :: last_history_temp(:)
   real(c_double), allocatable :: last_history_dt(:)
   real(c_double), allocatable :: last_history_dedt(:)
   real(c_double), allocatable :: last_history_xnuc(:,:)
   integer(c_int), allocatable :: last_history_substeps(:)
   type(neutrino) :: nu_saved

   common / cnetw / an, zn

contains

   subroutine reset_history(capacity)
      integer(c_int), intent(in) :: capacity

      if (allocated(last_history_time)) deallocate(last_history_time)
      if (allocated(last_history_rho)) deallocate(last_history_rho)
      if (allocated(last_history_temp)) deallocate(last_history_temp)
      if (allocated(last_history_dt)) deallocate(last_history_dt)
      if (allocated(last_history_dedt)) deallocate(last_history_dedt)
      if (allocated(last_history_substeps)) deallocate(last_history_substeps)
      if (allocated(last_history_xnuc)) deallocate(last_history_xnuc)

      last_history_count = 0_c_int
      last_history_capacity = max(capacity, 0_c_int)
      if (last_history_capacity <= 0_c_int .or. nactive <= 0_c_int) return

      allocate(last_history_time(last_history_capacity))
      allocate(last_history_rho(last_history_capacity))
      allocate(last_history_temp(last_history_capacity))
      allocate(last_history_dt(last_history_capacity))
      allocate(last_history_dedt(last_history_capacity))
      allocate(last_history_substeps(last_history_capacity))
      allocate(last_history_xnuc(nactive, last_history_capacity))
   end subroutine reset_history

   subroutine append_history(time, rho, temp_k, dt, dedt, substeps, yps)
      real(r8), intent(in) :: time, rho, temp_k, dt, dedt, yps(nsp)
      integer, intent(in) :: substeps
      integer :: i, row

      if (.not. allocated(last_history_xnuc)) return
      if (last_history_count >= last_history_capacity) return

      row = int(last_history_count) + 1
      last_history_count = last_history_count + 1_c_int
      last_history_time(row) = time
      last_history_rho(row) = rho
      last_history_temp(row) = temp_k
      last_history_dt(row) = dt
      last_history_dedt(row) = dedt
      last_history_substeps(row) = int(substeps, c_int)
      do i = 1, nactive
         last_history_xnuc(i, row) = yps(active_species(i))
      end do
   end subroutine append_history

   integer(c_int) function nuppn_init() bind(C)
      real(r8) :: t9, rho, dt_yr, t_max_yr, dt_max_yr, dt_factor
      integer :: i

      if (initialized) then
         nuppn_init = 0_c_int
         return
      end if

      nuppn_init = 1_c_int
      an = ZERO
      zn = ZERO
      qi = ZERO

      call reaction_info_init()
      call vital_init()
      call fuller_init()
      call other_nuc_init()
      call readframeinput(t9, rho, dt_yr, t_max_yr, dt_max_yr, dt_factor)
      call readsolverinput()
      call readphysicsinput()

      if (allocated(niso)) deallocate(niso)
      if (index_reaclib == 3) then
         allocate(niso(0:i325dim, 0:iCfdim, 2))
      else
         allocate(niso(0:i282dim, 0:iAtdim, 2))
      end if

      istart = 0
      call screen_init()
      call rates_init()
      call reaclib_init()
      call netgen_init()
      call alpha_decays_init()
      call neutrinos_init(nu_saved)
      call jbj_init()
      call nkk_init()
      call bader_deuflhard_init()
      call kadonis_init()

      if (ininet == 4) then
         nuppn_init = 4_c_int
         return
      end if
      call rnetw2007(ye_nw_ini, qi, an, zn, nvar_saved, nvrel_saved, &
         rho_nw_ini, t9_nw_ini, yps_nw_ini)

      istart = 1
      call nuc_data_init()
      call solver_init()
      call makemap()
      if (detailed_balance) call reverse_init()
      call reaclib_create_masks()
      if (use_cache) call reaclib_preprocessor()

      nactive = 0_c_int
      ppn_to_active = -1_c_int
      do i = 1, nsp
         if (considerisotope(i)) then
            nactive = nactive + 1_c_int
            active_species(nactive) = i
            ppn_to_active(i) = nactive - 1_c_int
         end if
      end do

      nactive_reactions = 0_c_int
      do i = 1, nre
         if (considerreaction(i)) then
            nactive_reactions = nactive_reactions + 1_c_int
            active_reactions(nactive_reactions) = i
         end if
      end do

      initialized = .true.
      nuppn_init = 0_c_int
   end function nuppn_init

   integer(c_int) function nuppn_num_species() bind(C)
      nuppn_num_species = nactive
   end function nuppn_num_species

   integer(c_int) function nuppn_get_species(index, name, name_len, z, a) bind(C)
      integer(c_int), value :: index, name_len
      character(kind=c_char) :: name(*)
      integer(c_int), intent(out) :: z, a
      integer :: ppn_index

      if (.not. initialized .or. index < 0_c_int .or. index >= nactive) then
         z = 0_c_int
         a = 0_c_int
         call copy_c_string('', name, name_len)
         nuppn_get_species = 1_c_int
         return
      end if

      ppn_index = active_species(index + 1_c_int)
      z = int(zn(ppn_index), c_int)
      a = int(an(ppn_index), c_int)
      call copy_c_string(zis(ppn_index), name, name_len)
      nuppn_get_species = 0_c_int
   end function nuppn_get_species

   integer(c_int) function nuppn_integrate(rho0, temp0_k, rho1, temp1_k, &
      xnuc, count, dt, dedt) bind(C)
      real(c_double), value :: rho0, temp0_k, rho1, temp1_k, dt
      integer(c_int), value :: count
      real(c_double) :: xnuc(*)
      real(c_double), intent(out) :: dedt
      real(r8) :: yps(nsp), ye, rho, t9, t9_0, t9_1
      integer :: ierr, i

      dedt = 0.0_c_double
      last_substeps = 0_c_int
      nuppn_integrate = unpack_xnuc(xnuc, count, yps)
      if (nuppn_integrate /= 0_c_int) return

      rho = HALF * (rho0 + rho1)
      t9_0 = temp0_k / 1.0e9_r8
      t9_1 = temp1_k / 1.0e9_r8
      t9 = HALF * (t9_0 + t9_1)
      call calculate_ye(yps, an, zn, ye, considerisotope)
      call evaluate_all_rates(ye, nvar_saved, nvrel_saved, rho, t9, yps, nu_saved)
      dedt = energy_flux(yps)
      call integrate_network(nvar_saved, yps, t9_0, t9_1, rho0, rho1, ye, dt, &
         nvrel_saved, nu_saved, ierr, ZERO)
      if (ierr /= 0) then
         nuppn_integrate = int(ierr, c_int)
         return
      end if

      last_substeps = int(nsubt, c_int)
      do i = 1, nactive
         xnuc(i) = yps(active_species(i))
      end do
      nuppn_integrate = 0_c_int
   end function nuppn_integrate

   integer(c_int) function nuppn_integrate_to_time(rho, temp_k, xnuc, count, &
      final_time, initial_dt, max_dt, dt_factor, max_steps, dedt) bind(C)
      real(c_double), value :: rho, temp_k, final_time, initial_dt, max_dt, dt_factor
      integer(c_int), value :: count, max_steps
      real(c_double) :: xnuc(*)
      real(c_double), intent(out) :: dedt
      real(r8) :: yps(nsp), ye, t9, time, step_dt, dt
      integer :: ierr, i, step, total_substeps

      dedt = 0.0_c_double
      last_substeps = 0_c_int
      nuppn_integrate_to_time = unpack_xnuc(xnuc, count, yps)
      if (nuppn_integrate_to_time /= 0_c_int) return

      t9 = temp_k / 1.0e9_r8
      time = ZERO
      step_dt = min(initial_dt, max_dt)
      total_substeps = 0
      call reset_history(max_steps + 1_c_int)
      call append_history(time, rho, temp_k, step_dt, ZERO, 0, yps)
      do step = 1, max_steps
         if (time >= final_time) exit
         dt = min(step_dt, final_time - time)
         if (dt <= ZERO) exit

         call calculate_ye(yps, an, zn, ye, considerisotope)
         call evaluate_all_rates(ye, nvar_saved, nvrel_saved, rho, t9, yps, nu_saved)
         dedt = energy_flux(yps)
         call integrate_network(nvar_saved, yps, t9, t9, rho, rho, ye, dt, &
            nvrel_saved, nu_saved, ierr, time)
         if (ierr /= 0) then
            last_substeps = int(total_substeps, c_int)
            nuppn_integrate_to_time = int(ierr, c_int)
            return
         end if

         total_substeps = total_substeps + nsubt
         time = time + dt
         call append_history(time, rho, temp_k, dt, dedt, nsubt, yps)
         step_dt = min(dt * dt_factor, max_dt)
      end do

      if (time < final_time) then
         last_substeps = int(total_substeps, c_int)
         nuppn_integrate_to_time = 3_c_int
         return
      end if

      last_substeps = int(total_substeps, c_int)
      do i = 1, nactive
         xnuc(i) = yps(active_species(i))
      end do
      nuppn_integrate_to_time = 0_c_int
   end function nuppn_integrate_to_time

   integer(c_int) function nuppn_last_substeps() bind(C)
      nuppn_last_substeps = last_substeps
   end function nuppn_last_substeps

   integer(c_int) function nuppn_history_size() bind(C)
      nuppn_history_size = last_history_count
   end function nuppn_history_size

   integer(c_int) function nuppn_get_history_step(index, time, rho, temp_k, &
      dt, dedt, substeps, xnuc, count) bind(C)
      integer(c_int), value :: index, count
      real(c_double), intent(out) :: time, rho, temp_k, dt, dedt
      integer(c_int), intent(out) :: substeps
      real(c_double) :: xnuc(*)
      integer :: i, row

      time = ZERO
      rho = ZERO
      temp_k = ZERO
      dt = ZERO
      dedt = ZERO
      substeps = 0_c_int
      if (.not. allocated(last_history_xnuc) .or. index < 0_c_int .or. &
         index >= last_history_count .or. count /= nactive) then
         nuppn_get_history_step = 1_c_int
         return
      end if

      row = int(index) + 1
      time = last_history_time(row)
      rho = last_history_rho(row)
      temp_k = last_history_temp(row)
      dt = last_history_dt(row)
      dedt = last_history_dedt(row)
      substeps = last_history_substeps(row)
      do i = 1, nactive
         xnuc(i) = last_history_xnuc(i, row)
      end do
      nuppn_get_history_step = 0_c_int
   end function nuppn_get_history_step

   integer(c_int) function nuppn_compute_nse(rho, temp_k, ye, xnuc, count) bind(C)
      real(c_double), value :: rho, temp_k, ye
      integer(c_int), value :: count
      real(c_double) :: xnuc(*)
      real(r8) :: yps(nsp), mu_p, mu_n
      integer :: ierr, iter, i

      if (.not. initialized) then
         nuppn_compute_nse = 1_c_int
         return
      end if
      if (count /= nactive) then
         nuppn_compute_nse = 2_c_int
         return
      end if

      if (.not. nse_ready) then
         call nse_init()
         nse_ready = .true.
      end if

      call compute_nse(temp_k / 1.0e9_r8, rho, ye, yps, mu_p, mu_n, iter, ierr)
      if (ierr /= 0) then
         nuppn_compute_nse = int(ierr, c_int)
         return
      end if

      do i = 1, nactive
         xnuc(i) = yps(active_species(i))
      end do
      nuppn_compute_nse = 0_c_int
   end function nuppn_compute_nse

   integer(c_int) function nuppn_get_dxdt(rho, temp_k, xnuc, count, dxdt_out) bind(C)
      real(c_double), value :: rho, temp_k
      integer(c_int), value :: count
      real(c_double), intent(in) :: xnuc(*)
      real(c_double) :: dxdt_out(*)
      real(r8) :: yps(nsp), ye, dt_est, t9
      integer :: i

      nuppn_get_dxdt = unpack_xnuc(xnuc, count, yps)
      if (nuppn_get_dxdt /= 0_c_int) return

      t9 = temp_k / 1.0e9_r8
      call calculate_ye(yps, an, zn, ye, considerisotope)
      call evaluate_all_rates(ye, nvar_saved, nvrel_saved, rho, t9, yps, nu_saved)
      call calculate_dxdt(yps, dt_est)
      do i = 1, nactive
         dxdt_out(i) = dxdt(active_species(i))
      end do
      nuppn_get_dxdt = 0_c_int
   end function nuppn_get_dxdt

   integer(c_int) function nuppn_evaluate_reactions(rho, temp_k, xnuc, count) &
      bind(C)
      real(c_double), value :: rho, temp_k
      integer(c_int), value :: count
      real(c_double), intent(in) :: xnuc(*)
      real(r8) :: yps(nsp), ye, t9

      nuppn_evaluate_reactions = unpack_xnuc(xnuc, count, yps)
      if (nuppn_evaluate_reactions /= 0_c_int) return

      t9 = temp_k / 1.0e9_r8
      call calculate_ye(yps, an, zn, ye, considerisotope)
      call evaluate_all_rates(ye, nvar_saved, nvrel_saved, rho, t9, yps, nu_saved)
      call compute_reaction_flows(yps)
      nuppn_evaluate_reactions = 0_c_int
   end function nuppn_evaluate_reactions

   integer(c_int) function nuppn_num_reactions() bind(C)
      nuppn_num_reactions = nactive_reactions
   end function nuppn_num_reactions

   integer(c_int) function nuppn_get_reaction(index, in1, in1_count, in2, &
      in2_count, out1, out1_count, out2, out2_count, rate, flow, q_value, &
      weak) bind(C)
      integer(c_int), value :: index
      integer(c_int), intent(out) :: in1, in1_count, in2, in2_count
      integer(c_int), intent(out) :: out1, out1_count, out2, out2_count, weak
      real(c_double), intent(out) :: rate, flow, q_value
      integer :: reaction_index

      in1 = -1_c_int
      in2 = -1_c_int
      out1 = -1_c_int
      out2 = -1_c_int
      in1_count = 0_c_int
      in2_count = 0_c_int
      out1_count = 0_c_int
      out2_count = 0_c_int
      rate = 0.0_c_double
      flow = 0.0_c_double
      q_value = 0.0_c_double
      weak = 0_c_int

      if (.not. initialized .or. index < 0_c_int .or. &
         index >= nactive_reactions) then
         nuppn_get_reaction = 1_c_int
         return
      end if

      reaction_index = active_reactions(index + 1_c_int)
      in1 = active_index(k1(reaction_index))
      in2 = active_index(k3(reaction_index))
      out1 = active_index(k7(reaction_index))
      out2 = active_index(k5(reaction_index))
      in1_count = int(k2(reaction_index), c_int)
      in2_count = int(k4(reaction_index), c_int)
      out1_count = int(k8(reaction_index), c_int)
      out2_count = int(k6(reaction_index), c_int)
      rate = real(v(reaction_index), c_double)
      flow = real(reaction_flow(reaction_index), c_double)
      q_value = real(bind_energy_diff(reaction_index), c_double)
      if (is_weak_reaction(reaction_index)) weak = 1_c_int
      nuppn_get_reaction = 0_c_int
   end function nuppn_get_reaction

   integer(c_int) function unpack_xnuc(xnuc, count, yps)
      integer(c_int), value :: count
      real(c_double), intent(in) :: xnuc(*)
      real(r8), intent(out) :: yps(nsp)
      integer :: i

      if (.not. initialized) then
         unpack_xnuc = 1_c_int
         return
      end if
      if (count /= nactive) then
         unpack_xnuc = 2_c_int
         return
      end if

      yps = ZERO
      do i = 1, nactive
         yps(active_species(i)) = xnuc(i)
      end do
      unpack_xnuc = 0_c_int
   end function unpack_xnuc

   subroutine compute_reaction_flows(yps)
      real(r8), intent(in) :: yps(nsp)
      real(r8) :: xmol(nsp)
      integer :: gamma, i

      xmol = ZERO
      reaction_flow = ZERO
      where (considerisotope .and. an /= ZERO)
         xmol = yps / an
      end where
      gamma = ispe('OOOOO')
      if (gamma > 0) xmol(gamma) = ONE

      do i = 1, nre
         if (.not. considerreaction(i)) cycle
         if (k1(i) < 1 .or. k1(i) > nsp) cycle
         if (k3(i) < 1 .or. k3(i) > nsp) cycle
         reaction_flow(i) = input_factor(k2(i)) * v(i) * &
            xmol(k1(i)) ** k2(i) * xmol(k3(i))
      end do
   end subroutine compute_reaction_flows

   real(r8) function energy_flux(yps)
      real(r8), intent(in) :: yps(nsp)
      integer :: i

      call compute_reaction_flows(yps)
      energy_flux = ZERO
      do i = 1, nre
         if (.not. considerreaction(i)) cycle
         energy_flux = energy_flux + reaction_flow(i) * bind_energy_diff(i)
      end do
   end function energy_flux

   real(r8) function input_factor(count)
      integer, intent(in) :: count
      select case (count)
      case (0, 1)
         input_factor = ONE
      case (2)
         input_factor = HALF
      case (3)
         input_factor = ONE / 6.0_r8
      case (4)
         input_factor = ONE / 24.0_r8
      case default
         input_factor = ZERO
      end select
   end function input_factor

   integer(c_int) function active_index(ppn_index)
      integer, intent(in) :: ppn_index

      if (ppn_index < 1 .or. ppn_index > nsp) then
         active_index = -1_c_int
      else
         active_index = ppn_to_active(ppn_index)
      end if
   end function active_index

   logical function is_weak_reaction(reaction_index)
      integer, intent(in) :: reaction_index

      select case (ilabb(reaction_index))
      case (i_bm, i_ec, i_bn, i_bp, i_ba, i_ve, i_vp, i_vsp, i_vsn, i_vsa)
         is_weak_reaction = .true.
      case default
         is_weak_reaction = .false.
      end select
   end function is_weak_reaction

   subroutine copy_c_string(text, out, out_len)
      character(len=*), intent(in) :: text
      character(kind=c_char) :: out(*)
      integer(c_int), value :: out_len
      integer :: i, n

      if (out_len <= 0_c_int) return
      do i = 1, out_len
         out(i) = c_null_char
      end do
      n = min(len_trim(text), int(out_len) - 1)
      do i = 1, n
         out(i) = text(i:i)
      end do
   end subroutine copy_c_string

end module imyann_nuppn_c_api
