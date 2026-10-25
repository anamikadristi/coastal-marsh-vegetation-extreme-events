# Collapse CRMS hourly hydrographic rows to one row per station and day:
# salinity sum and count, water elevation relative to the marsh surface sum and count,
# and hours with water above the marsh surface (flooded hours).
# Input columns (CIMS export): 1 Station ID, 2 Date, 5 Sensor Environment, 11 Adjusted Salinity (ppt), 15 Adjusted Water Elevation to Marsh (ft)
BEGIN{FS=","; OFS=","; print "station,date,sal_sum,sal_n,wlm_sum,wlm_n,flood_h"}
NR==1{next}
$1 !~ /-H[0-9][0-9]$/ {next}
$5 != "Surface Water" {next}
{
  k=$1 FS $2
  if(k!=prev && prev!=""){ flush() }
  prev=k
  if($11!=""){ss+=$11; sn++}
  if($15!=""){ws+=$15; wn++; if($15>0) fh++}
}
function flush(){ print prev, ss, sn, ws, wn, fh; ss=sn=ws=wn=fh=0 }
END{ if(prev!="") flush() }
