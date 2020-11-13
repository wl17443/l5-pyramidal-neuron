
using Plots

# Units 
uF = 1e-6
mV = 1e-3
ms = 1e-3
uA = 1e-6
mS = 1e-3

# A few constants...
Cm = 100*uF # Membrane Conductance 
p = 0.5   # Somatic area/Total area 
dt = 0.5*ms 

beta_m = -1.2*mV 
gamma_m = 18*mV
minf(Vs) = 0.5 * ( 1 + tanh((Vs - beta_m) / gamma_m))

gc = 1*mS 
Ids(Vd, Vs) = gc * (Vd - Vs)

Ena = 55*mV 
gna = 45*mS
Ina(Vs) = gna * minf(Vs) * (Vs - Ena)

beta_w = 0*mV
gamma_w = 10*mV
winf(Vs) = 0.5 * (1 + tanh((Vs - beta_w) / gamma_w))

tau_w(Vs) = 1/cosh((Vs-beta_w) / (2 * gamma_w))

phi_w = 0.15
fw(Vs, w) = phi_w * (winf(Vs) - w) / tau_w(Vs)

gk = 18*mS 
Ek = -80*mV
Ik(Vs, w) = gk * w * (Vs - Ek)

gsl = 0.1 * mS
Esl = -65*mV 
Isl(Vs) = gsl * (Vs - Esl) 

function fvs(Vs, Vd, Is, w)
  return (Is/p + Ids(Vd, Vs)/p - Ina(Vs) - Ik(Vs, w) - Isl(Vs)) / Cm 
end 

ninf(Vd) = 1 / (1 + exp(-(Vd + 9) / 0.5))

tau_n = 15*ms 
fn(Vd, n) = (ninf(Vd) - n) / tau_n

hinf(Vd) = 1 / (1 + exp((Vd + 21) / 0.5))

tau_h = 80*ms 
fh(Vd, h) = (hinf(Vd) - h) / tau_h

gca = 0.8*mS
Eca = 140*mV
Ica(Vd, n, h) = gca * n * h * (Vd - Eca)

gdl = 0.1*mS
Edl = -65*mV 
Idl(Vd) = gdl * (Vd - Edl)

function fvd(Vd, Vs, Id, n, h)
  return (Id/(1-p) - Ids(Vd, Vs)/(1-p) - Ica(Vd, n, h) - Idl(Vd)) / Cm
end 

Vrests = -60*mV 
Vrestd = -70*mV 
Vth    = -47*mV

function simulate(T=5, dt=1*ms)
  voltage_trace_s = zeros(Int(T/dt))
  voltage_trace_d = zeros(Int(T/dt))

  voltage_trace_s[1] = -64.8*mV
  voltage_trace_d[1] = -64.8594*mV

  steady_input_current_s = 34*uA
  steady_input_current_d = 0*uA

  input_current_s = zeros(Int(T/dt))
  input_current_d = zeros(Int(T/dt))
  for i in 1:Int(T/dt)
    if i >= 100 && i < 200
      input_current_s[i] = steady_input_current_s
    end
  end 

  w = 0.049
  n = 0.9650
  h = 1

  for timestep in 2:Int(T/dt)
    voltage_s = voltage_trace_s[timestep-1]
    voltage_d = voltage_trace_d[timestep-1]
    println("v_s: $voltage_s, v_d: $voltage_d, w: $w, n: $n, h: $h")
    
    if voltage_trace_s[timestep-1] >= Vth
      w = w + fw(voltage_trace_s[timestep-1], w) * dt
      voltage_trace_s[timestep-1] = 0
      voltage_trace_s[timestep] = Vrests
    else 
      w = w + fw(voltage_trace_s[timestep-1], w) * dt
      voltage_trace_s[timestep] = voltage_trace_s[timestep-1] + fvs(voltage_trace_s[timestep-1], voltage_trace_d[timestep-1], input_current_s[timestep-1], w) * dt
    end 

    if voltage_trace_d[timestep-1] >= Vth
      n = n + fn(voltage_trace_d[timestep-1], n) * dt
      h = h + fh(voltage_trace_d[timestep-1], h) * dt
      voltage_trace_d[timestep-1] = 0
      voltage_trace_d[timestep] = Vrestd
    else 
      n = n + fn(voltage_trace_d[timestep-1], n) * dt
      h = h + fh(voltage_trace_d[timestep-1], h) * dt 
      voltage_trace_d[timestep] = voltage_trace_d[timestep-1] + fvd(voltage_trace_d[timestep-1], voltage_trace_s[timestep-1], input_current_d[timestep-1], n, h) * dt
    end
  end 
  return (voltage_trace_d, voltage_trace_s)
end 

voltage_d, voltage_s = simulate()
gr()
x = range(1, stop=Int(5/(1*ms)))
y = voltage_s
z = voltage_d
plt = plot(x, y)
# plt2 = plot!(plt, x, z)
gui(plt)
