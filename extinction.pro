;Model after Hayes and Latham 1975
;for mean values see Green (1992), cfa-www.harvard.edu/icq/ICQExtinct.html


pro extinction, alt, wavelength_Kieffer, A_prime

  altitude = alt/1000.d0

  ;lambda = 0.51d0	;wavelength in um
  lambda = wavelength_Kieffer/1000.d0

; ozone extinction, value from Schaefer 1992: http://www.cfa.harvard.edu/icq/ICQExtinct.html

  A_oz = 0.016d0


; Rayleigh Scattering

  c1 = 0.23465d0 + (107.6d0/(146.0d0-(1./lambda)^2.d0)) + (0.93161d0/(41.0d0-(1.d0/lambda)^2.d0))
  A_ray = 0.0094977d0/lambda^4.d0*c1*c1*exp(-altitude/7.996d0)


; Aerosol scattering

  A_aer = 0.05d0*lambda^(-1.3d0)*exp(-altitude/1.5d0)

  A_prime = A_ray + A_aer + A_oz

;Burki et al. 1995, p7: for La Silla
  ;A_ray = 0.00685d0*lambda^(-4.05d0) ;nearly identical to model
  ;A_aer = 0.0116d0*lambda^(-1.39d0)	;factor 1.42 higher

return
end