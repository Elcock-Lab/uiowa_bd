
      subroutine update_counters
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        

      nstep=nstep+1                                                 
      necount=necount+1
      ntcount=ntcount+1
      nmcount=nmcount+1
      tot_time=tot_time+time_step                                       
      if(num_umb_cnt.eq.num_umb_stp) num_umb_cnt=0
      if(num_lst_cnt.eq.num_lst_stp) num_lst_cnt=0
      if(num_fmd_cnt.eq.num_fmd_stp) num_fmd_cnt=0
      if(num_hyd_cnt.eq.num_hyd_stp) num_hyd_cnt=0
      if(num_bal_cnt.eq.num_bal_stp) num_bal_cnt=0
      if(num_ene_scrn_cnt.eq.1000) num_ene_scrn_cnt=0                 
      num_ene_scrn_cnt=num_ene_scrn_cnt+1                           
      num_umb_cnt=num_umb_cnt+1                          
      num_lst_cnt=num_lst_cnt+1                          
      num_fmd_cnt=num_fmd_cnt+1                          
      num_hyd_cnt=num_hyd_cnt+1                          
      num_bal_cnt=num_bal_cnt+1                          

      return
      end

