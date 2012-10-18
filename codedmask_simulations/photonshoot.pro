;******************************************************************************
;creates photons on the source shoots them over the mask and holds
;placement data on the detector for the ones who passes the mask
;------------------------------------------------------------------------------
pro photonshoot,nofphot,maskhit,dethit,image,pixenergy,aperture, $
                osource=source,odetector=detector,omask=mask,maxdelta=deltamax,$
                plot=plot,intsource=intsource
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;Yigit Dallilar 14.10.2012
;INPUT
;nofphot     : number of photons to be simulated
;OUTPUT
;maskhit     : mask hit positions
;dethit      : detector hit positions
;detector    : detector structure
;OPTIONAL
;maxdelta    : maximum delta range on the denominators
;osource     : photon source can be added from outside
;odetector   : detector can be added from outside
;omask       : mask can be added from outside
;KEYWORD
;plot        : realise some plots
;intsource   : tries intelligent source who looks to the center of the mask
;------------------------------------------------------------------------------

;variable initialisation
  if not keyword_set(source) then $
     source = create_struct('radius',100,'theta',0.3*!pi,'phi',0.75*!pi, $
                            'pos',[5.5,-5.5,5.5],'energy',1)
  if not keyword_set(detector) then $
     detector=create_struct('z',-20,'pixsize',1.,'pixenergy',dblarr(200,200))
  tmp = sqrt(n_elements(detector.pixenergy)) 
  detector=create_struct(detector,'pos',dblarr(tmp,tmp,2))
  if not keyword_set(mask) then getmask,[2,2],73,mask,pixsize=1.
  if not keyword_set(deltamax) then deltamax=0.

  radius=[1,1,sqrt(randomu(systime(1),nofphot+1))]*(source.radius)
  angle=[!pi*0.2,-!pi*0.2,randomu(systime(1)+1,nofphot+1)*2*!pi]
  delta=(randomu(systime(1)+2,nofphot,3)-0.5)*deltamax
  xp=radius[*]*cos(angle[*])
  yp=radius[*]*sin(angle[*])
  refpos=dblarr(2,3)
  photenergy=dblarr(nofphot)+source.energy
  pos=dblarr(nofphot,3)
  maskhit=dblarr(nofphot,2)
  maskpix=lonarr(nofphot,2)
  dethit=dblarr(nofphot,2)
  denom=dblarr(3)
  tarray = (findgen(tmp)-(tmp/2.-0.5))*detector.pixsize
  for i=0,tmp-1 do detector.pos[*,i,0]=tarray
  for i=0,tmp-1 do detector.pos[i,*,1]=tarray

  if keyword_set(intsource) then begin
     intsource=create_struct('radius',100000,'theta',!pi*0.2,'phi',-!pi*0.5)
     source.pos[0]=intsource.radius*sin(intsource.theta)*cos(intsource.phi)
     source.pos[1]=intsource.radius*sin(intsource.theta)*sin(intsource.phi)
     source.pos[2]=intsource.radius*cos(intsource.theta)
     source.theta=!pi*0.5-intsource.theta
     source.phi=!pi+intsource.phi
     denom[0]=source.pos[0]
     denom[1]=source.pos[1]
     denom[2]=source.pos[2]
  endif else begin
;calculation of denominators with two reference positions in the circle 
     for i=0,1 do begin
        refpos[i,0] = +xp[i]*sin(source.phi)+ $
                      yp[i]*sin(source.theta)*cos(source.phi)+source.pos[0]
        refpos[i,1] = -xp[i]*cos(source.phi)+ $
                      yp[i]*sin(source.theta)*sin(source.phi)+source.pos[1]
        refpos[i,2] = yp[i]*cos(source.theta)+source.pos[2]
     endfor
     denom[2]=(refpos[1,0]-source.pos[0])*(refpos[0,1]-source.pos[1]) $
              -(refpos[1,1]-source.pos[1])*(refpos[0,0]-source.pos[0])      
     denom[0]=(refpos[1,1]-source.pos[1])*(refpos[0,2]-source.pos[2]) $
              -(refpos[1,2]-source.pos[2])*(refpos[0,1]-source.pos[1])
     denom[1]=(refpos[1,2]-source.pos[2])*(refpos[0,0]-source.pos[0]) $
              -(refpos[1,0]-source.pos[0])*(refpos[0,2]-source.pos[2])
  endelse

;cartesian coordinates and denominator calculation
  for i=0,nofphot-1 do begin
     pos[i,0] = +xp[i+2]*sin(source.phi)+ $
                yp[i+2]*sin(source.theta)*cos(source.phi)+source.pos[0]
     pos[i,1] = -xp[i+2]*cos(source.phi)+ $
                yp[i+2]*sin(source.theta)*sin(source.phi)+source.pos[1]
     pos[i,2] = yp[i+2]*cos(source.theta)+source.pos[2]
  endfor

;mask hit positions and pixels are written z=0
  for i=0,nofphot-1 do begin
     maskhit[i,0]= $
        -((denom[0]*(1+delta[i,0]))/(denom[2]*(1+delta[i,2])))*pos[i,2]+pos[i,0]
     maskhit[i,1]= $
        -((denom[1]*(1+delta[i,1]))/(denom[2]*(1+delta[i,2])))*pos[i,2]+pos[i,1]
     maskpix[i,0]=where(maskhit[i,0] le mask.pos[*,0,0]+mask.pixsize*0.5 and $
                maskhit[i,0] ge mask.pos[*,0,0]-mask.pixsize*0.5)
     maskpix[i,1]=where(maskhit[i,1] le mask.pos[0,*,1]+mask.pixsize*0.5 and $
                maskhit[i,1] ge mask.pos[0,*,1]-mask.pixsize*0.5)               
  endfor

;detector hit positions and energy distributed to pixels are calculated
  for i=0,nofphot-1 do begin
     dethit[i,0]=(((denom[0]*(1+delta[i,0]))/(denom[2]*(1+delta[i,2])))* $
                  (-pos[i,2]+detector.z)+pos[i,0])* $
                 mask.apert[maskpix[i,0],maskpix[i,1]]
     dethit[i,1]=(((denom[1]*(1+delta[i,1]))/(denom[2]*(1+delta[i,2])))* $
                  (-pos[i,2]+detector.z)+pos[i,1])* $
                 mask.apert[maskpix[i,0],maskpix[i,1]]
     if mask.apert[maskpix[i,0],maskpix[i,1]] gt 0 then begin
        ndx=where(dethit[i,0] le detector.pos[*,0,0]+detector.pixsize*0.5 and $
                  dethit[i,0] ge detector.pos[*,0,0]-detector.pixsize*0.5)
        ndy=where(dethit[i,1] le detector.pos[0,*,1]+detector.pixsize*0.5 and $
                  dethit[i,1] ge detector.pos[0,*,1]-detector.pixsize*0.5)
        detector.pixenergy[ndx,ndy]=detector.pixenergy[ndx,ndy]+photenergy[i]
     endif
  endfor

  ;output value 
  pixenergy=detector.pixenergy
  aperture=mask.apert
  image=convol(pixenergy,aperture)

  ;plotting image
  if keyword_set(plot) then begin
     window,0,xsize=1200,ysize=600
     !p.multi=[0,4,2]
     contour,(image/max(image))^1,nlevel=100,/fill,xr=[75,125],yr=[75,125]
     contour,(image/max(image))^4,nlevel=100,/fill,xr=[75,125],yr=[75,125]
     contour,(image/max(image))^7,nlevel=100,/fill,xr=[75,125],yr=[75,125]
     surface,(image/max(image))^1
     contour,detector.pixenergy,nlevel=100,/fill
     range=[-tmp*detector.pixsize*0.5,tmp*detector.pixsize*0.5]
     plot,dethit[*,0],dethit[*,1],psym=3,xr=range,yr=range
     plot,maskhit[*,0],maskhit[*,1],psym=3,xr=range,yr=range,color=1000
     getmask,[2,2],73,/plot,pixsize=mask.pixsize
     plot,radius*cos(angle),radius*sin(angle),psym=3
  endif

end
;******************************************************************************
