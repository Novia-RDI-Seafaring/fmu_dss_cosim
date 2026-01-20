model IEEE_14_Buses_FMU
  /**
     * TS_FMU: Wrapper around example1_1 (Transmission) using 'extends' so B2.p is at top level.
     * PCC = Bus2.
     */
  extends OpenIPSL.Examples.IEEE14.IEEE_14_Buses();
  import PwPinFictitiousPort;
  // Fictitious line port at the PCC (Bus2)
  parameter Real Yc_re = 0 "Real part of line admittance (pu)";
  parameter Real Yc_im = 0 "Imag part of line admittance (pu)";
  PwPinFictitiousPort pcc(Yc_re = Yc_re, Yc_im = Yc_im);
  // I/O seen by the FMU master
  input Real Hin_re;
  input Real Hin_im;
  output Real V_re;
  output Real V_im;
  output Real I_re;
  output Real I_im;
equation
// Now 'B2.p' is directly in scope thanks to 'extends'
  connect(pcc.p, B7.p);
// pass-through I/O
  pcc.Hin_re = Hin_re;
  pcc.Hin_im = Hin_im;
  V_re = pcc.V_re;
  V_im = pcc.V_im;
  I_re = pcc.I_re;
  I_im = pcc.I_im;
  annotation(
    experiment(StartTime = 0, StopTime = 5, Tolerance = 1e-6, Interval = 0.001),
    uses(OpenIPSL(version = "3.1.0-dev")));
end IEEE_14_Buses_FMU;
