rm_onlyDigits=function(digCnt,TotCnt,x)
{
  rtString=vector()
  k=1
  for(i in 1:length(digCnt))
  {
    if((TotCnt[i]) != (digCnt[i]))
    {
      rtString[k]=x[i]
      k=k+1
    }
    
  }
  return(rtString)
}

