model Two_Areas_TS_FMU
  "TS FMU: Sds + h_in inputs, Vpcc + Im + h_out outputs (bus7 is PCC)"

  import Modelica.Blocks.Interfaces.*;

// ============ FMU INPUTS ============
  RealInput Sds_re "Real part of S_DS (pu)";
  RealInput Sds_im "Imag part of S_DS (pu)";

  RealInput hin_re "Incoming history term h_in real part";
  RealInput hin_im "Incoming history term h_in imag part";

// ============ FMU OUTPUTS ============
  RealOutput Vpcc_re "TS PCC voltage real part (pu)";
  RealOutput Vpcc_im "TS PCC voltage imag part (pu)";

  RealOutput Im_re "Interface current Im (TS -> DS) real part (pu)";
  RealOutput Im_im "Interface current Im (TS -> DS) imag part (pu)";

  RealOutput hout_re "Outgoing history term h_out real part";
  RealOutput hout_im "Outgoing history term h_out imag part";

// ============ EXTEND THE TS MODEL ============
  // This pulls in Two_Areas_PSSE_AVR, including bus7, lines, etc.
  extends OpenIPSL.Examples.TwoAreas.Two_Areas_PSSE_AVR;

  // ============ DS LOAD INJECTION ============
  DSLoadFromS dsLoad
    annotation (Placement(transformation(origin = {-108, -44}, extent = {{-20, -20}, {20, 20}})));

// ============ FICTITIOUS LINE PARAMETERS ============
  // Yc is the characteristic admittance of the fictitious line (pu)
  // These are PARAMETERS: you will set them from Python based on Zc_pu.
  parameter Real Yc_re = 0.0 "Real part of Yc (pu)";
  parameter Real Yc_im = 0.0 "Imag part of Yc (pu)";

protected 
  Real Vk_re;
  Real Vk_im;
  Real Im_loc_re;
  Real Im_loc_im;

equation
// ---- connect DS load to bus7 and feed Sds ----
  connect(dsLoad.p, bus7.p);
  dsLoad.S_re = Sds_re;
  dsLoad.S_im = Sds_im;
// ---- PCC voltage from bus7 pin ----
  Vk_re = bus7.p.vr;
  Vk_im = bus7.p.vi;

  Vpcc_re = Vk_re;
  Vpcc_im = Vk_im;
// ---- Interface current Im (TS -> DS) ----
// dsLoad.p.ir/ii are currents INTO the DSLoad component.
// For co-simulation we interpret Im as current from TS to DS,
// which is the same direction as p.ir (into the load).
  Im_loc_re = dsLoad.p.ir;
  Im_loc_im = dsLoad.p.ii;

  Im_re = Im_loc_re;
  Im_im = Im_loc_im;
// ---- History update: H_out = Yc * Vk + Im ----
// This matches the definition used in the initialization script:
//   Hk = Yc * Vk + Ik
// All quantities are current-like (pu):
  hout_re = (Yc_re*Vk_re - Yc_im*Vk_im) + Im_loc_re;
  hout_im = (Yc_re*Vk_im + Yc_im*Vk_re) + Im_loc_im;

end Two_Areas_TS_FMU;
