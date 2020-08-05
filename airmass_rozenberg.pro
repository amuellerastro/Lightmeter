; compute airmass from Rozenberg 1966, see Krisciunas, Schaefer (1991): A Model of the Brightness of Moonlight
; limited to 40 air masses

function airmass_rozenberg, zd

  AM = 1.d0/(cos(zd*!dtor) + (0.025d0*exp(-11.0d0*cos(zd*!dtor))))

return, AM

end