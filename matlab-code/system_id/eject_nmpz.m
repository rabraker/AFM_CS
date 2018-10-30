function sys_eject = eject_nmpz(sys)
  
  z = zero(sys);
  z_nmp = z(abs(z) > 1);
  
  g_eject = zpk([], z_nmp, 1, sys.Ts);
  g_eject = g_eject/dcgain(g_eject);
  
  sys_eject = minreal(g_eject*sys);
  
  
end