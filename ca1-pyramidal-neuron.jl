# Units
mV = 1^-3 
mS = 1^-3
uF = 1^-6

# Neuron parameters of CA1 two-compartment model 
neuronParameters = Dict( "VNa"   => 120.0*mV, 
                         "VCa"   => 140.0*mV, 
                         "VK"    => -15.0*mV,
                         "Vleak" => 0.0*mV,
                         "VExc"  => 60.0*mV,
                         "gc"    => 1.5*mS,
                         "p"     => 0.5,
                         "Cm"    => 3.0*uF)

ionicConductances = Dict( "gleakS" => 0.1*mS,
                          "gNaS"   => 30.0*mS, 
                          "gKDRS"  => 17.0*mS,
                          "gKAHPS" => 0.8*mS,
                          "gKCS"   => 15.0*mS,
                          "gleakD" => 0.1*mS,
                          "gNaD"   => 0.0*mS,
                          "gCaD"   => 5.0*mS,
                          "gKDRD"  => 0.0*mS,
                          "gKAHPD" => 0.8*mS,
                          "gKCD"   => 5.0*mS)

calcium_conc = Dict( "Ca2+S" => 0,
                     "Ca2+D" => 0)

                
function Ileak(compartment, voltage)
  return ionicConductances["gleak$compartment"]*(voltage-neuronParameters["Vleak"])
end 

function INa(compartment, voltage, h)
  return ionicConductances["gNa$compartment"]*minf^2*h*(voltage-neuronParameters["VNa"])
end 

function IKDR(compartment, voltage, n)
  return ionicConductances["gKDR$compartment"]*n*(voltage-neuronParameters["VK"])
end 

function ICa(compartment, voltage, s)
  return ionicConductances["gCa$compartment"]*s^2*(voltage-neuronParameters["VCa"])
end 

function IKAHP(compartment, voltage, q)
  return ionicConductances["gKAHP$compartment"]*q*(voltage-neuronParameters["VK"])
end 

function IKC(compartment, calcium_conce, voltage, c)
  return ionicConductances["gKC$compartment"]*c*min(1, calcium_conce/250)*(voltage-neuronParameters["VK"])
end 

# Functions with respect to time 
function dCa_dt(calcium_current, calcium_conce, phi, beta)
  return -phi*calcium_current-beta*calcium_conce
end 

# Rate constants 
function alpha(constant, voltage)
  if constant == 'm'
    return 0.32*(13.1-voltage)/(exp^((13.1-voltage)/4)-1)
  elseif constant == 'h'
    return 0.128*exp^((17-voltage)/18)
  elseif constant == 'n'
    return (0.016*(35.1-voltage))/(exp^((35.1-voltage)/5)-1)
  elseif constant == 's'
    return 1.6/(1+exp^(-0.072*(voltage-65)))
  elseif constant == 'c'
    if voltage > 50
      return 2*exp^((6.5-v)/27)
    else 
      return exp^((voltage-10)/11-(voltage-6.5)/27)/18.975
    end 
  end 
end 

function beta(constant, voltage)