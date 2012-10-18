;******************************************************************************
; this a num x num mura aperture function generator
;------------------------------------------------------------------------------
pro muragen,num,mura,pos,pixsize=sizepix,plotmura=plotmura
;------------------------------------------------------------------------------
;Yigit Dallilar 11.11.2012
;INPUTS 
;num     : for num x num array
;OUTPUT
;mura    : mura num x num aperture function ( 1 for open , 0 for closed)
;pos     : position of every pixels of mura
;OPTIONAL OUTPUT
;pixsize : one pixel size for created mura
;------------------------------------------------------------------------------

;variable initialisation  
  if keyword_set(plotmura) then pmura=1 else pmura=0
  if not keyword_set(sizepix) then sizepix=1
  tmp=''
  mura = bytarr(num,num)
  pos = dblarr(num,num,2)
  cntrl=0              ;0 if there is no num existed in the file,1 for found

;writes postion data
  for i=0,num-1 do begin 
     pos[*,i,0]=(findgen(num)-(num-1)/2)*sizepix
     pos[i,*,1]=((num-1)/2-findgen(num))*sizepix
  endfor

;reads file
  openr,1,'quadresidu.txt'
  while ~ eof(1) do begin
     readf,1,tmp
     tmp1 = strsplit(tmp,' ',/extract)
     if long(tmp1(0)) eq num then begin
        cntrl=1
        break
     endif
  endwhile
  close,1

;checks if num found then starts to create mura
  if cntrl eq 1 then begin
     arr = byte(tmp1(1))  ;arr has quadresidu data
     mura[0,1:num-1]=1
     for i=0,num-2 do begin
        for j=0,num-2 do begin
           if arr[i] eq arr[j] then mura[i+1,j+1] = 1 $
           else mura[i+1,j+1] = 0
        endfor
     endfor
;draws mura
     if pmura eq 1 then begin 
        plotmura,mura,pos,sizepix
     endif
        
  endif

end
;*******************************************************************************
