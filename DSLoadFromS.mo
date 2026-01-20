
model DSLoadFromS
  "Draws complex power S (pu) from a PwPin using S = V*conj(I), Re/Im only"

  // Inputs from FMU/master
  Modelica.Blocks.Interfaces.RealInput S_re "Real part of S (pu)";
  Modelica.Blocks.Interfaces.RealInput S_im "Imag part of S (pu)";

  // Electrical connector to the bus
  OpenIPSL.Interfaces.PwPin p;

protected 
  Real P;
  Real Q;
  Real Vr;
  Real Vi;
  Real V2;
  Real Ir;
  Real Ii;

equation
  // Complex power S = P + j Q
  P = S_re;
  Q = S_im;

  // Local bus voltage from pin
  Vr = p.vr;
  Vi = p.vi;

  // |V|^2 (small epsilon avoids division by zero)
  V2 = Vr*Vr + Vi*Vi + 1e-6;

  // I = conj(S / V) in real arithmetic:
  // Re(S/V) = (P*Vr + Q*Vi)/|V|^2
  // Im(S/V) = (Q*Vr - P*Vi)/|V|^2
  // I = conj(S/V) => Ir = Re(S/V), Ii = -Im(S/V)
  Ir = (P*Vr + Q*Vi) / V2;
  Ii = (P*Vi - Q*Vr) / V2;

  // Currents into the component
  p.ir = Ir;
  p.ii = Ii;

end DSLoadFromS;
