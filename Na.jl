function Na(Vs)
  am = Vs == -40 ? 1 : 0.1*(Vs+40)/(1-exp(-(Vs-40)/10))
  bm = 4*exp(-(Vs+65)/18)
  ah = 0.07*exp(-(Vs+65)/20)
  bh = 1/(1+exp(-(Vs+35)/10))

  minf = am/(am+bm)
  taum = 1/(am+bm)
  hinf = ah/(ah+bh)
  tauh = 1/(ah+bh)

  return minf, taum, hinf, tauh 
end 

