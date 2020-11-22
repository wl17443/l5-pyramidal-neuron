
using Plots

# Units 
uF = 1e-6
mV = 1
ms = 1e-3
uA = 1e-6
mS = 1e-3
uS = 1e-6
um = 1e-6
MOhm = 1e6
nM = 1e-9
zF = 1e-21
mM = 1e-3

T = 5
dt = 0.5e-3

# A few constants...
# Cm = 2e-6 # Membrane Conductance 
Cm = Dict("s" => 0.26e-9, "d" => 0.12e-9)
Rm = Dict("s" => 50e6, "d" => 43e6)
Rt = 65e6
Tadj = 2.3^((34-21)/10)
Vsh = 8 # Shift in the dendritic channels kinetics to account fo rht eshift in the dendritic resting potential 

Vs_0 = -65.0 
Vd_0 = -55.0 


actvar = Dict("mNa"  => 0.0, 
              "hNa"  => 0.0, 
              "hKdr" => 0.0,
              "mNap" => 0.0,
              "hNap" => 0.0,
              "Ca2+" => 76.6e-6,
              "mCaL" => 76.6e-6, 
              "mKs"  => 0.0,
              "hKs"  => 0.0,
              "mh"   => 0.0,
              "mM"   => 0.0)

# Ionic Channels (Somatic)

# Sodium 
amNa(Vs) = Vs == -40 ? 1 : (Vs+40)/(10*(1-exp(-(Vs-40)/10)))
bmNa(Vs) = 4*exp(-(Vs+65)/18)
ahNa(Vs) = 0.07*exp(-(Vs+65)/20)
bhNa(Vs) = 1/(1+exp(-(Vs+35)/10))

minfNa(Vs) = amNa(Vs)/(amNa(Vs)+bmNa(Vs))
taumNa(Vs) = 1/(amNa(Vs)+bmNa(Vs))
hinfNa(Vs) = ahNa(Vs)/(ahNa(Vs)+bhNa(Vs))
tauhNa(Vs) = 1/(ahNa(Vs)+bhNa(Vs))

gbar_Na = 18e-6; E_Na = 50
INa(Vs) = gbar_Na*(actvar["mNa"]^3)*actvar["hNa"]*(Vs-E_Na)

# Potassium 
ahKdr(Vs) = Vs == -55 ? 0.1 : 0.01*(Vs+55)/(1-exp(-(Vs+55)/10))
bhKdr(Vs) = 0.125*exp(-(Vs+65)/80)

hinfKdr(Vs) = ahKdr(Vs)/(ahKdr(Vs)+bhKdr(Vs))
tauhKdr(Vs) = 1e-3/(ahKdr(Vs)+bhKdr(Vs))


gbar_Kdr = 5e-6; E_Kdr = -85 
IKdr(Vs) = gbar_Kdr*(actvar["hKdr"]^4)*(Vs-E_Kdr)

# Leak Current 
E_sl = -31.5
Isl(Vs) = (Vs-E_sl)/Rm["s"]

function fvs(Vs, Vd, Is)
  return (-INa(Vs)-IKdr(Vs)-Isl(Vs)+(Vd-Vs)/Rt+Is)/Cm["s"]
end 

amNap(Vd) = abs(Vd+38.0) < 2.22e-16 ? 0.182*6/exp(-(Vd+38)/6) : 0.182*(Vd+38)/(1-exp(-(Vd+38)/6))
bmNap(Vd) = abs(Vd+38.0) < 2.22e-16 ? 0.124*6/exp((Vd+38)/6) : -0.124*(Vd+38)/(1-exp((Vd+38)/6));
minfNap(Vd) = 1/(1+exp(-(Vd+52.6)/4.6))
# minfNap(Vd) = amNap(Vd)/(amNap(Vd)+bmNap(Vd))
taumNap(Vd) = 1e-3*6/(Tadj*(amNap(Vd)+bmNap(Vd)))

ahNap(Vd) = abs(Vd+17) < 2.22e-16 ? 2.88e-6*4.63/exp((Vd+17)/4.63) : -2.88e-6*(Vd+17)/(1-exp((Vd+17)/4.63))
bhNap(Vd) = abs(Vd+64.4) < 2.22e-16 ? 6.94e-6*2.63/exp(-(Vd+64.4)/2.63) : 6.94e-6*(Vd+64.4)/(1-exp(-(Vd+64.4)/2.63))
hinfNap(Vd) = 1/(1+exp((Vd+48.8)/10))
# hinfNap(Vd) = ahNap(Vd)/(ahNap(Vd)+bhNap(Vd))
tauhNap(Vd) = 1e-3/(Tadj*(ahNap(Vd)+bhNap(Vd)))

gbar_Nap = 0.022e-6; E_Nap = 50 
INap(Vd) = gbar_Nap*(actvar["mNap"]^3)*actvar["hNap"]*(Vd-E_Nap)

gbar_Ca = 3.85e-6 ; R = 8.314e3; z_Ca = 2; F = 9.648e4; temp = 310.15; tau_R = 80e-3; Ca_o = 2; 
Ca_const = (R*temp)/(z_Ca*F)
E_Ca = Ca_const*log(Ca_o/actvar["Ca2+"])

amCaL(Vd) = 1.6/(exp(-0.072*(Vd-5))+1)
bmCaL(Vd) = abs(Vd+8.69) < 2.22e-16 ? 0.02*5.36/exp((Vd+8.69)/5.36) : 0.02*(Vd+8.69)/(exp((Vd+8.69)/5.36)-1)

minfCaL(Vd) = amCaL(Vd)/(amCaL(Vd)+bmCaL(Vd))
taumCaL(Vd) = 1e-3/(amCaL(Vd)+bmCaL(Vd))

ICaL(Vd) = gbar_Ca*(actvar["mCaL"]^2)*(Vd-E_Ca)
Ad = 9302.3e-8; d = 0.1; gamma = 0.02; factorCa = 1e4*gamma/(Ad*d); Cai_inf = 76.6e-6
Cam_0 = minfCaL(Vd_0)
ICa_0 = gbar_Ca*(Cam_0^2)*(Vd_0-E_Ca)

minfKs(Vd) = 1/(1+exp(-(Vd+11)/12))
taumKs(Vd) = Vd < -50 ? 1e-3*(1.25+175.03*exp(0.026*(Vd+10)))/Tadj : 1e-3*(1.25+13*exp(-0.026*(Vd+10)))/Tadj
hinfKs(Vd) = 1/(1+exp((Vd+64)/11))
tauhKs(Vd) = 1e-3*(360+(1010+24*(Vd+65))*exp(-((Vd+85)/48)^2))/Tadj

gbar_Ks = 28e-6; E_Ks = -85
IKs(Vd) = gbar_Ks*(actvar["mKs"]^2)*actvar["hKs"]*(Vd-E_Ks)

amh(Vd) = isnan((6.43*(Vd+154))/(exp((Vd+154)/11.9)-1)) ? (6.43)/((1/11.9)*exp((Vd+154)/11.9)) : 6.43*(Vd+154)/(exp((Vd+154)/11.9)-1)
bmh(Vd) = 193*exp(Vd/33.1)

minfh(Vd) = amh(Vd)/(amh(Vd)+bmh(Vd))
taumh(Vd) = 1/(amh(Vd)+bmh(Vd))

gbar_h = 0.865e-6; E_Ih = -45
Ih(Vd) = gbar_h*actvar["mh"]*(Vd-E_Ih)

amM(Vd) = 0.0033*exp(0.1*(Vd+35))
bmM(Vd) = 0.0033*exp(-0.1*(Vd+35))

minfM(Vd) = amM(Vd)/(amM(Vd)+bmM(Vd))
taumM(Vd) = 1e-3/(Tadj*(amM(Vd)+bmM(Vd)))

gbar_IM = 1e-6
IM(Vd) = gbar_IM*actvar["mM"]*(Vd-E_Ks)

E_dl = -48.1
Idl(Vd) = (Vd-E_dl)/Rm["d"]

function fvd(Vd, Vs, Id)
  return (-INap(Vd)-ICaL(Vd)-Ih(Vd)-IM(Vd)-IKs(Vd)-Idl(Vd)+(Vs-Vd)/Rt+Id)/Cm["d"]
end 

fmNa(Vs) = (minfNa(Vs)-actvar["mNa"])/taumNa(Vs)
fhNa(Vs) = (hinfNa(Vs)-actvar["hNa"])/tauhNa(Vs)
fhKdr(Vs) = (hinfKdr(Vs)-actvar["hKdr"])/tauhKdr(Vs)

fmNap(Vd) = (minfNap(Vd)-actvar["mNap"])/taumNap(Vd)
fhNap(Vd) = (hinfNap(Vd)-actvar["hNap"])/tauhNap(Vd)
fCa2Plus(Vd) = -factorCa*((ICaL(Vd)-ICa_0)/(z_Ca*F))-(actvar["Ca2+"]-Cai_inf)/tau_R
fmCaL(Vd) = (minfCaL(Vd)-actvar["mCaL"])/taumCaL(Vd)
fmKs(Vd) = (minfKs(Vd)-actvar["mKs"])/taumKs(Vd)
fhKs(Vd) = (hinfKs(Vd)-actvar["hKs"])/tauhKs(Vd)
fmh(Vd) = (minfh(Vd)-actvar["mh"])/taumh(Vd)
fmM(Vd) = (minfM(Vd)-actvar["mM"])/taumM(Vd)

function simulate()
  voltage_trace_s = zeros(Int(ceil(T/dt)))
  voltage_trace_d = zeros(Int(ceil(T/dt)))

  voltage_trace_s[1] = Vs_0
  voltage_trace_d[1] = Vd_0

  steady_input_current_s = 20e-6
  steady_input_current_d = 0e-6

  input_current_s = zeros(Int(ceil(T/dt)))
  input_current_d = zeros(Int(ceil(T/dt)))

  actvar["mNa"] = minfNa(Vs_0); actvar["hNa"] = hinfNa(Vs_0); actvar["hKdr"] = hinfKdr(Vs_0); actvar["mNap"] = minfNap(Vd_0-Vsh);
  actvar["hNap"] = hinfNap(Vd_0-Vsh); actvar["mKs"] = minfKs(Vd_0-Vsh); actvar["hKs"] = hinfKs(Vd_0-Vsh); actvar["mh"] = minfh(Vd_0-Vsh);
  actvar["mM"] = minfM(Vd_0-Vsh); 

  for i in 1:Int(ceil(T/dt))
    if i >= 1000 && i < 3000
      input_current_s[i] = steady_input_current_s
    end
  end 

  for i in 1:Int(ceil(T/dt))
    if i >= 3000 && i < 7000
      input_current_s[i] = steady_input_current_s
    end
  end 

  for timestep in 2:Int(ceil(T/dt))
    voltage_s = voltage_trace_s[timestep-1]
    voltage_d = voltage_trace_d[timestep-1]
    println("v_s: $voltage_s, v_d: $voltage_d")
    
    # Update activation variables 
    actvar["mNa"] = actvar["mNa"]+fmNa(voltage_trace_s[timestep-1])*dt
    actvar["hNa"] = actvar["hNa"]+fhNa(voltage_trace_s[timestep-1])*dt
    actvar["hKdr"] = actvar["hKdr"] + fhKdr(voltage_trace_s[timestep-1])*dt
    voltage_trace_s[timestep] = voltage_trace_s[timestep-1] + fvs(voltage_trace_s[timestep-1], voltage_trace_d[timestep-1], input_current_s[timestep-1]) * dt
    println("actvar: $actvar")

    actvar["mNap"] = actvar["mNap"]+fmNap(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["hNap"] = actvar["hNap"]+fhNap(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["mCaL"] = actvar["mCaL"]+fmCaL(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["Ca2+"] = actvar["Ca2+"]+fCa2Plus(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["mh"] = actvar["mh"]+fmh(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["mM"] = actvar["mM"]+fmM(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["mKs"] = actvar["mKs"]+fmKs(voltage_trace_d[timestep-1]-Vsh)*dt
    actvar["hKs"] = actvar["hKs"]+fhKs(voltage_trace_d[timestep-1]-Vsh)*dt
    voltage_trace_d[timestep] = voltage_trace_d[timestep-1]+fvd(voltage_trace_d[timestep-1], voltage_trace_s[timestep-1], input_current_d[timestep-1]) * dt
  end 
  return (voltage_trace_d, voltage_trace_s)
end 

voltage_d, voltage_s = simulate()
gr()
x = range(1, stop=Int(ceil(T/dt)))
y = voltage_s
z = voltage_d
plt = plot(x, y, label="Somatic Voltage")
plt2 = plot(x, z, label="Dendritic Voltage")

plot(plt, plt2, layout=(2,1))

